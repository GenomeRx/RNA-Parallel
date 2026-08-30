## The parallel machinery.
##
## Backends, chunk arithmetic, and what happens when a worker dies or a backend misbehaves.
## These do not test the statistics; they test that the dispatch layer cannot corrupt a
## result or swallow a failure.

## Which multiprocessing frameworks are compatible, asserted rather than claimed.
## Every backend must reproduce the serial result bit for bit. All parallelism
## funnels through one dispatch point, so this is the test that keeps that true.

test_that("the backend list is what the docs promise", {
  expect_setequal(combat_backends(),
                  c("mclapply", "future", "BiocParallel", "foreach", "serial"))
})

test_that("an unknown backend is refused, not silently ignored", {
  d <- make_counts(20, G = 60, n_per_batch = c(4, 4))
  expect_error(ComBat_seq_parallel(d$counts, d$batch, parallel_backend = "sparklyr"),
               "should be one of")
})

test_that("every installed backend gives a bit-identical result", {
  vendor <- backend_fn()
  d <- make_counts(21, G = 300, n_per_batch = c(7, 6))
  ref <- quietly(vendor(d$counts, d$batch, group = NULL))

  needs <- c(mclapply = "parallel", future = "future.apply",
             BiocParallel = "BiocParallel", foreach = "doParallel", serial = "base")

  exercised <- 0L
  for (be in combat_backends()) {
    pkg <- needs[[be]]
    # skip() aborts the whole test, so one absent package used to drop every backend after
    # it from this comparison. next drops only that backend, and the counter below refuses a
    # run in which nothing was exercised.
    if (pkg != "base" && !requireNamespace(pkg, quietly = TRUE)) next
    exercised <- exercised + 1L
    if (be == "future") {
      old <- future::plan(future::multisession, workers = 2)
      on.exit(future::plan(old), add = TRUE)
    }
    got <- quietly(ComBat_seq_parallel(d$counts, d$batch, group = NULL,
                                       workers = 2L, chunks = 4L,
                                       parallel_backend = be))
    expect_identical(got, ref, info = paste("backend:", be))
  }
  # a loop that exercised nothing would otherwise report green
  expect_gte(exercised, 2L)
})

test_that("chunk order is preserved, which is what identical() depends on", {
  # A backend returning results out of order would rbind the genes scrambled. The
  # matrix would still be the right shape, so only a value check catches it.
  d <- make_counts(22, G = 240, n_per_batch = c(6, 6))
  ref <- quietly(ComBat_seq_parallel(d$counts, d$batch, group = NULL,
                                     parallel_backend = "serial", workers = 1L))
  for (be in c("mclapply", "BiocParallel")) {
    if (be == "BiocParallel" && !requireNamespace("BiocParallel", quietly = TRUE)) next
    got <- quietly(ComBat_seq_parallel(d$counts, d$batch, group = NULL,
                                       workers = 4L, chunks = 8L,
                                       parallel_backend = be))
    expect_identical(rownames(got), rownames(d$counts))
    expect_identical(got, ref, info = paste("backend:", be))
  }
})

test_that("combat.fork = FALSE forces serial on every backend", {
  d <- make_counts(23, G = 150, n_per_batch = c(5, 5))
  ref <- quietly(ComBat_seq_parallel(d$counts, d$batch, group = NULL,
                                     parallel_backend = "serial", workers = 1L))
  old <- getOption("combat.fork")
  options(combat.fork = FALSE)
  on.exit(options(combat.fork = old), add = TRUE)
  for (be in combat_backends()) {
    got <- quietly(ComBat_seq_parallel(d$counts, d$batch, group = NULL,
                                       workers = 4L, parallel_backend = be))
    expect_identical(got, ref, info = paste("backend:", be))
  }
})

test_that("the future backend warns instead of silently running serially", {
  skip_if_not_installed("future.apply")
  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)
  idx <- rnaparallel:::combat_row_chunks(20L, chunks = 4L)
  expect_warning(
    rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), workers = 2L,
                                               parallel_backend = "future"),
    "resolves in one process")
})

test_that("a custom executor function is accepted and gives identical results", {
  skip_on_os("windows")   # these executors call mclapply, which cannot fork on Windows

  # The extension point. Anything with an lapply shape plugs in, which is how the
  # many frameworks were verified this way without the package
  # growing a branch per framework.
  vendor <- backend_fn()
  d <- make_counts(24, G = 250, n_per_batch = c(6, 6))
  ref <- quietly(vendor(d$counts, d$batch, group = NULL))

  # preschedule=TRUE, which the named mclapply backend does not use
  presched <- function(idx, f, workers)
    parallel::mclapply(idx, f, mc.cores = min(workers, length(idx), test_max_cores()), mc.preschedule = TRUE)
  expect_identical(
    quietly(ComBat_seq_parallel(d$counts, d$batch, group = NULL,
                                workers = 2L, parallel_backend = presched)), ref)

  # a plain serial closure, proving no framework is required at all
  expect_identical(
    quietly(ComBat_seq_parallel(d$counts, d$batch, group = NULL, workers = 2L,
                                parallel_backend = function(idx, f, workers) lapply(idx, f))), ref)
})

test_that("a custom executor returning the wrong shape is caught, not silently used", {
  idx <- rnaparallel:::combat_row_chunks(20L, chunks = 4L)
  f <- function(i) sum(i)

  # drops a chunk: would rbind a short matrix downstream
  expect_error(
    rnaparallel:::combat_parallel_lapply(
      idx, f, workers = 2L,
      parallel_backend = function(idx, f, workers) lapply(idx[-1], f)),
    "must return a list of length 4")

  # returns a vector rather than a list
  expect_error(
    rnaparallel:::combat_parallel_lapply(
      idx, f, workers = 2L,
      parallel_backend = function(idx, f, workers) unlist(lapply(idx, f))),
    "must return a list")
})

test_that("combat.fork = FALSE overrides a custom executor too", {
  # the escape hatch must not be defeatable by supplying your own executor
  called <- FALSE
  spy <- function(idx, f, workers) { called <<- TRUE; lapply(idx, f) }
  old <- getOption("combat.fork")
  options(combat.fork = FALSE)
  on.exit(options(combat.fork = old), add = TRUE)
  idx <- rnaparallel:::combat_row_chunks(10L, chunks = 2L)
  rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), workers = 4L,
                                             parallel_backend = spy)
  expect_false(called)
})

test_that("a dispatch inside a worker does not open a second pool", {
  # ComBat-seq dispatches the tagwise loop across BATCHES and ships the vendor closure, whose
  # environment still carries the rebound estimateGLMTagwiseDisp; inside the worker that symbol
  # dispatched AGAIN over gene rows. Measured on Windows before the guard: workers = 2L gave 2
  # outer and 4 nested processes, which is workers + workers^2, or 272 at the 16-worker arm.
  #
  # The old guard was mc.allow.recursive = FALSE, an argument to parallel::mclapply, so it
  # covered the fork branch alone. Windows runs mclapply serially and reaches its workers only
  # through foreach/PSOCK, so the one platform that needed the guard never had it.
  #
  # A spy executor rather than real workers: it runs in-process, so this is deterministic and
  # costs nothing, and it counts DISPATCHES, which is the thing that must not happen twice.
  dispatches <- 0L
  spy <- function(idx, f, workers) { dispatches <<- dispatches + 1L; lapply(idx, f) }

  inner <- function(i) {
    rnaparallel:::combat_parallel_lapply(
      rnaparallel:::combat_row_chunks(10L, chunks = 2L), function(j) sum(j),
      workers = 4L, parallel_backend = spy, cells = Inf, min_cells = 0)
    sum(i)
  }
  out <- rnaparallel:::combat_parallel_lapply(
    rnaparallel:::combat_row_chunks(20L, chunks = 2L), inner,
    workers = 2L, parallel_backend = spy, cells = Inf, min_cells = 0)

  # exactly one: the outer dispatch. Both inner calls must have taken the serial path.
  expect_identical(dispatches, 1L)
  expect_length(out, 2L)
  # and the flag must not survive into the caller's session
  expect_identical(Sys.getenv("RNAPARALLEL_IN_WORKER"), "")
})

test_that("the nesting guard does not disturb a caller's own parallel loop", {
  # The documented pattern is to spend the worker budget INSIDE a loop body. Nothing marks a
  # caller's own workers, so a companion called from within their loop must still dispatch.
  dispatches <- 0L
  spy <- function(idx, f, workers) { dispatches <<- dispatches + 1L; lapply(idx, f) }
  # a caller's loop that is not ours, then our dispatch inside it
  invisible(lapply(1:2, function(k)
    rnaparallel:::combat_parallel_lapply(
      rnaparallel:::combat_row_chunks(10L, chunks = 2L), function(j) sum(j),
      workers = 2L, parallel_backend = spy, cells = Inf, min_cells = 0)))
  expect_identical(dispatches, 2L)
})

test_that("the two nesting guards answer differently for a caller's own fork", {
  # The spy tests above run a CUSTOM executor, which never reaches the mclapply branch, so
  # neither of them exercises mc.allow.recursive = FALSE. That guard is still passed, and on
  # the default backend it fires for ANY enclosing fork child, including one the caller made.
  # The env-var guard fires only for our own workers. The two therefore disagree about a
  # caller's own loop, and the README says which is which, so it is pinned here by process
  # count rather than left to prose.
  skip_on_os("windows")
  skip_if(!identical(Sys.getenv("NOT_CRAN"), "true"), "forking test")

  pids <- function() rnaparallel:::combat_parallel_lapply(
    rnaparallel:::combat_row_chunks(8L, chunks = 4L), function(j) Sys.getpid(),
    workers = 4L, parallel_backend = "mclapply", cells = Inf, min_cells = 0)

  # at top level the dispatch forks, so the chunks report several distinct PIDs
  expect_gt(length(unique(unlist(pids()))), 1L)

  # inside the caller's OWN fork child it degrades to serial: one PID, the child's
  inner <- parallel::mclapply(1:2, function(k) length(unique(unlist(pids()))), mc.cores = 2L)
  expect_identical(unlist(inner), c(1L, 1L))
})

test_that("a dispatch too small to be worth a fork runs serially instead", {
  # ComBat-seq dispatches once per batch, so a 100-batch design hands over thin
  # slices. Below the threshold the fork cost more than the work: measured 0.80x
  # against plain sva on 500 genes x 1000 samples x 100 batches.
  spy_calls <- 0L
  spy <- function(idx, f, workers) { spy_calls <<- spy_calls + 1L; lapply(idx, f) }
  idx <- rnaparallel:::combat_row_chunks(100L, chunks = 4L)
  f <- function(i) sum(i)

  old <- getOption("combat.min.cells")
  on.exit(options(combat.min.cells = old), add = TRUE)
  options(combat.min.cells = 2e4)

  rnaparallel:::combat_parallel_lapply(idx, f, workers = 4L,
                                             parallel_backend = spy, cells = 5000)
  expect_identical(spy_calls, 0L)

  rnaparallel:::combat_parallel_lapply(idx, f, workers = 4L,
                                             parallel_backend = spy, cells = 5e5)
  expect_identical(spy_calls, 1L)

  # the gate is opt-out, and 0 restores the old unconditional dispatch
  options(combat.min.cells = 0)
  rnaparallel:::combat_parallel_lapply(idx, f, workers = 4L,
                                             parallel_backend = spy, cells = 5000)
  expect_identical(spy_calls, 2L)
})

test_that("the size gate does not change the numbers, only who computes them", {
  skip_on_os("windows")   # these executors call mclapply, which cannot fork on Windows

  d <- make_counts(11, G = 300L, n_per_batch = c(12, 10, 11))
  old <- getOption("combat.min.cells")
  on.exit(options(combat.min.cells = old), add = TRUE)

  options(combat.min.cells = 0)
  hot <- quietly(ComBat_seq_parallel(d$counts, d$batch, group = NULL, workers = 2L))
  options(combat.min.cells = 1e9)
  gated <- quietly(ComBat_seq_parallel(d$counts, d$batch, group = NULL, workers = 2L))

  expect_identical(gated, hot)
})

## ---- the compatibility matrix ------------------------------------------------
## One backend on one code path proves very little. These cross mechanism against
## argument path, including the two paths that consume RNG.

# spans both mechanisms on purpose: fork shares pages, sockets ship a copy
matrix_backends <- function() {
  b <- list(serial = "serial", mclapply = "mclapply")
  b[["presched"]] <- function(idx, f, w)
    parallel::mclapply(idx, f, mc.cores = min(w, length(idx), test_max_cores()), mc.preschedule = TRUE)

  # Created and torn down inside the call, so no socket is ever open while another backend in
  # this matrix forks. Holding one cluster across the sweep instead was measurably worse: a
  # forked child inherits the parent's socket file descriptors, and the suite segfaulted on a
  # later send. The package's own pool guards this by refusing a handle it did not create
  # (R/helper_seq_parallel.R, combat_cluster); a hand-rolled cluster here has no such guard, so
  # it simply must not outlive the call.
  #
  # Standing a fresh cluster up repeatedly is what made a worker connection break about one run
  # in five, so the creation retries rather than failing the suite on a port that was still in
  # TIME_WAIT.
  b[["parLapply_PSOCK"]] <- function(idx, f, w) {
    n <- min(w, length(idx), test_max_cores())
    # Create and tear down inside the call: a cluster that outlives it has its socket file
    # descriptors inherited by the forking backends this same matrix runs, and the suite
    # segfaulted on a later send. Standing one up per dispatch instead breaks a connection
    # mid-parLapply now and then under the churn, so the whole create-use-destroy cycle
    # retries rather than the creation alone. Retrying is safe because a dispatch is a pure
    # function of its chunk: the same idx recomputed gives the same answer, and nothing here
    # accumulates.
    for (attempt in 1:4) {
      cl <- tryCatch({
        cc <- parallel::makeCluster(n, type = "PSOCK")
        parallel::clusterEvalQ(cc, suppressMessages(library(edgeR)))
        cc
      }, error = function(e) NULL)
      if (is.null(cl)) { Sys.sleep(0.25 * attempt); next }
      out <- tryCatch(parallel::parLapply(cl, idx, f), error = function(e) e)
      try(parallel::stopCluster(cl), silent = TRUE)
      if (!inherits(out, "error")) return(out)
      # Only a broken socket is retryable. An error raised BY the chunk must surface on the
      # first attempt: the next test in this file deliberately makes a worker throw and
      # asserts no backend swallows it, and a blanket retry turned that assertion into a
      # skip, which is the exact failure it exists to catch.
      if (!grepl("writing to connection|reading from connection|invalid connection",
                 conditionMessage(out))) {
        stop(out)
      }
      Sys.sleep(0.25 * attempt)
    }
    testthat::skip("PSOCK could not complete a dispatch in four attempts")
  }
  if (requireNamespace("BiocParallel", quietly = TRUE)) b[["BiocParallel"]] <- "BiocParallel"
  b
}

test_that("every backend agrees on every argument path, RNG paths included", {
  skip_on_os("windows")   # matrix_backends()'s "presched" entry calls mclapply directly
  vendor <- backend_fn()
  d <- make_counts(30, G = 200, n_per_batch = c(7, 7, 6), with_group = TRUE)
  n <- ncol(d$counts)
  covar <- cbind(cov1 = rep_len(c(0, 1), n), cov2 = rep_len(c(0, 0, 1, 1), n))

  paths <- list(
    list(nm = "default",     a = list(group = NULL),                     seed = FALSE),
    list(nm = "group",       a = list(group = d$group),                  seed = FALSE),
    list(nm = "covar_mod",   a = list(group = NULL, covar_mod = covar),  seed = FALSE),
    list(nm = "full_mod=F",  a = list(group = d$group, full_mod = FALSE), seed = FALSE),
    # shrink routes through monte_carlo_int_NB, which calls sample(). If a backend
    # perturbed the RNG stream this is the cell that would catch it.
    list(nm = "shrink",      a = list(group = d$group, shrink = TRUE, shrink.disp = TRUE),
         seed = TRUE),
    list(nm = "subset.n",    a = list(group = d$group, shrink = TRUE, gene.subset.n = 40),
         seed = TRUE))

  bes <- matrix_backends()
  for (p in paths) {
    if (p$seed) set.seed(4242)
    ref <- quietly(do.call(vendor, c(list(counts = d$counts, batch = d$batch), p$a)))
    for (bn in names(bes)) {
      if (p$seed) set.seed(4242)
      got <- quietly(do.call(ComBat_seq_parallel, c(
        list(counts = d$counts, batch = d$batch, workers = 3L, chunks = 4L,
             parallel_backend = bes[[bn]]), p$a)))
      expect_identical(got, ref, info = paste0(bn, " / ", p$nm))
    }
  }
})

test_that("no backend fails silently when a chunk errors", {
  # A partial or wrong result returned quietly is far worse than an error. Some
  # backends throw, some hand back an error object; either is fine, silence is not.
  idx <- rnaparallel:::combat_row_chunks(40L, chunks = 4L)
  # keyed on a row that lands in exactly one chunk under ANY split, contiguous or
  # interleaved. The old predicate was i[1] > 10, which silently stopped firing when
  # chunking became round-robin because every chunk then starts at row 1..4.
  boom <- function(i) if (37L %in% i) stop("worker boom") else sum(i)

  bes <- matrix_backends()
  for (bn in names(bes)) {
    be <- bes[[bn]]
    surfaced <- tryCatch({
      parts <- rnaparallel:::combat_parallel_lapply(idx, boom, workers = 3L,
                                                          parallel_backend = be)
      # came back rather than throwing: the check must reject it
      tryCatch({ rnaparallel:::combat_parallel_check(parts, "probe"); FALSE },
               error = function(e) TRUE)
    }, error = function(e) TRUE)
    expect_true(surfaced, info = paste("backend swallowed a worker error:", bn))
  }
})

test_that("an active future plan does not disturb the other backends", {
  # A caller's session may set a plan globally; a fork backend must not be affected by it.
  skip_if_not_installed("future.apply")
  vendor <- backend_fn()
  d <- make_counts(31, G = 180, n_per_batch = c(6, 6))
  ref <- quietly(vendor(d$counts, d$batch, group = NULL))

  old <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old), add = TRUE)
  for (bn in c("mclapply", "serial")) {
    expect_identical(
      quietly(ComBat_seq_parallel(d$counts, d$batch, group = NULL, workers = 3L,
                                  parallel_backend = bn)),
      ref, info = paste(bn, "under an active multisession plan"))
  }
})

test_that("the option sets the default backend", {
  old <- getOption("combat.backend")
  options(combat.backend = "serial")
  on.exit(options(combat.backend = old), add = TRUE)
  expect_identical(formals(ComBat_seq_parallel)$parallel_backend |> eval(), "serial")
})

test_that("clusters are cached, not rebuilt per dispatch", {
  skip_on_os("windows")   # asks for a FORK cluster explicitly

  # The bug this guards: the foreach backend used to call makeCluster() inside every
  # dispatch. One ComBat-seq run dispatches at least five times, typically 3 + 2*n_batch,
  # and a 4-worker PSOCK cluster
  # costs ~150 ms to build, which made that backend measure 25x slower than serial.
  skip_on_cran()
  rnaparallel:::combat_cluster_stop()
  a <- rnaparallel:::combat_cluster(2L, "FORK")
  b <- rnaparallel:::combat_cluster(2L, "FORK")
  expect_identical(a, b)                       # same object, not a rebuild
  expect_gte(rnaparallel::combat_cluster_stop(), 1L)
})

test_that("a stopped cluster is rebuilt rather than reused dead", {
  skip_on_cran()
  rnaparallel:::combat_cluster_stop()
  cl <- rnaparallel:::combat_cluster(2L, "FORK")
  parallel::stopCluster(cl)                    # kill it behind the cache's back
  fresh <- rnaparallel:::combat_cluster(2L, "FORK")
  expect_true(inherits(fresh, "cluster"))
  expect_length(parallel::clusterCall(fresh, function() 1L), 2L)
  rnaparallel::combat_cluster_stop()
})

test_that("a machine with no ps binary still builds and reuses a cluster", {
  skip_on_cran()
  skip_on_os("windows")

  # The identity probe shells out to ps. Making a missing identity fatal turned ps into an
  # undeclared hard dependency: on a distroless container or a sandbox that blocks subprocess
  # spawning, every foreach dispatch died where 0.4.4 had worked. Identity is a hardening
  # bonus over the PID, not a precondition for building a cluster.
  script <- tempfile(fileext = ".R")
  writeLines(c(
    'Sys.setenv(PATH = tempdir())',
    sprintf('suppressMessages(pkgload::load_all(%s, quiet = TRUE))',
            deparse(normalizePath(testthat::test_path("..", "..")))),
    'stopifnot(!nzchar(Sys.which("ps")))',
    'cl <- rnaparallel:::combat_cluster(2L, "FORK")',
    'stopifnot(length(cl) == 2L)',
    'stopifnot(all(is.na(rnaparallel:::.combat_clusters$FORK$worker_identities)))',
    'stopifnot(length(rnaparallel:::combat_cluster(2L, "FORK")) == 2L)',
    'idx <- rnaparallel:::combat_row_chunks(12L, chunks = 3L)',
    'o <- rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), workers = 2L,',
    '       parallel_backend = "foreach", cells = Inf, min_cells = 0)',
    'stopifnot(identical(unlist(o), vapply(idx, sum, integer(1))))',
    'rnaparallel::combat_cluster_stop()',
    'cat("NO_PS_OK\\n")'), script)
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)),
    stdout = TRUE, stderr = TRUE,
    env = c("NOT_CRAN=true",
            paste0("R_LIBS=", shQuote(paste(.libPaths(), collapse = .Platform$path.sep)))),
    timeout = 60))
  expect_true("NO_PS_OK" %in% out, info = paste(out, collapse = "\n"))
})

test_that("a retired entry that will not die is given up on rather than retried forever", {
  skip_on_cran()
  skip_on_os("windows")

  outsider <- parallel::mcparallel(Sys.sleep(30))
  on.exit({
    if (isTRUE(tools::pskill(outsider$pid, 0L))) tools::pskill(outsider$pid, tools::SIGKILL)
    try(parallel::mccollect(outsider, wait = FALSE), silent = TRUE)
  }, add = TRUE)

  cache <- rnaparallel:::.combat_clusters
  old_retired <- cache$retired
  on.exit(cache$retired <- old_retired, add = TRUE)

  # a live PID this session recorded but cannot identify: status "unknown", so it is signalled
  # and waited on. Before the attempt cap it survived every round, and combat_cluster() paid
  # the full timeout on every later dispatch, forever, saying nothing.
  ns <- asNamespace("rnaparallel")
  old_kill <- get("combat_wait_pids", envir = ns)
  assignInNamespace("combat_wait_pids", function(...) FALSE, ns = "rnaparallel")
  on.exit(assignInNamespace("combat_wait_pids", old_kill, ns = "rnaparallel"), add = TRUE)

  cache$retired <- list(list(pid = Sys.getpid(), pids = outsider$pid,
                             identities = NA_character_, tries = 0L))
  tries <- rnaparallel:::.combat_retire_tries
  for (i in seq_len(tries - 1L)) {
    suppressWarnings(rnaparallel:::combat_retired_reap(timeout = 0))
    expect_length(cache$retired, 1L)
  }
  expect_warning(rnaparallel:::combat_retired_reap(timeout = 0), "gave up killing")
  expect_null(cache$retired)
})

test_that("a forked child never reaps the workers its parent recorded", {
  skip_on_cran()
  skip_on_os("windows")

  cache <- rnaparallel:::.combat_clusters
  old_retired <- cache$retired
  on.exit(cache$retired <- old_retired, add = TRUE)

  outsider <- parallel::mcparallel(Sys.sleep(30))
  on.exit({
    if (isTRUE(tools::pskill(outsider$pid, 0L))) tools::pskill(outsider$pid, tools::SIGKILL)
    try(parallel::mccollect(outsider, wait = FALSE), silent = TRUE)
  }, add = TRUE)

  # stamped with a PID that is not this process, exactly as a child inheriting the parent's
  # list would see it. The child must forget the entry without signalling anything behind it.
  cache$retired <- list(list(pid = Sys.getpid() + 1L, pids = outsider$pid,
                             identities = NA_character_, tries = 0L))
  expect_identical(rnaparallel:::combat_retired_reap(timeout = 0), 0L)
  expect_true(isTRUE(tools::pskill(outsider$pid, 0L)))
  expect_null(cache$retired)
})

test_that("teardown reports success on the fallback path, not only through stopCluster", {
  skip_on_cran()
  skip_on_os("windows")

  cl <- parallel::makeCluster(2L, type = "FORK")
  pids <- unlist(parallel::clusterCall(cl, Sys.getpid), use.names = FALSE)
  on.exit(for (p in pids) if (isTRUE(tools::pskill(p, 0L))) tools::pskill(p, tools::SIGKILL),
          add = TRUE)

  # suspected_dead skips stopCluster entirely, which is the path every partially dead cluster
  # takes. It closed every node correctly and still reported FALSE on macOS and Linux.
  expect_true(rnaparallel:::combat_cluster_teardown(cl, suspected_dead = TRUE))
})

test_that("combat_cluster_stop counts clusters, not reaped leftovers", {
  skip_on_cran()
  skip_on_os("windows")

  rnaparallel::combat_cluster_stop()
  cache <- rnaparallel:::.combat_clusters
  old_retired <- cache$retired
  on.exit({ cache$retired <- old_retired; rnaparallel::combat_cluster_stop() }, add = TRUE)

  # a leftover entry for a PID that is already gone used to be counted as a stopped cluster
  cache$retired <- list(list(pid = Sys.getpid(), pids = 999999L,
                             identities = NA_character_, tries = 0L))
  expect_identical(rnaparallel::combat_cluster_stop(), 0L)

  cl <- rnaparallel:::combat_cluster(2L, "FORK")
  expect_length(cl, 2L)
  expect_identical(rnaparallel::combat_cluster_stop(), 1L)
})

test_that("retired cleanup never kills a reused PID", {
  skip_on_cran()
  skip_on_os("windows")

  outsider <- parallel::mcparallel(Sys.sleep(30))
  on.exit({
    if (isTRUE(tools::pskill(outsider$pid, 0L))) {
      tools::pskill(outsider$pid, tools::SIGKILL)
    }
    try(parallel::mccollect(outsider, wait = FALSE), silent = TRUE)
  }, add = TRUE)

  cache <- rnaparallel:::.combat_clusters
  old_retired <- cache$retired
  on.exit(cache$retired <- old_retired, add = TRUE)
  cache$retired <- list(list(
    pid = Sys.getpid(),
    pids = outsider$pid,
    identities = "not-the-outsider"
  ))

  expect_identical(rnaparallel:::combat_retired_reap(timeout = 0), 1L)
  # NOT pskill(pid, 0L): an mcparallel child that has been SIGKILLed but not collected is a
  # zombie, and signal 0 still succeeds against it. That assertion returns TRUE for exactly
  # the process this test exists to prove was left alone. A zombie's command line reads
  # <defunct>, so the package's own identity probe classifies it as foreign instead.
  now <- rnaparallel:::combat_pid_identities(outsider$pid)
  expect_identical(rnaparallel:::combat_pid_status(outsider$pid, now), "owned")
  expect_false(grepl("defunct", now))
  expect_null(cache$retired)
})

test_that("a partially dead cluster is torn down whole, leaving no worker behind", {
  skip_on_cran()
  skip_on_os("windows")

  script <- testthat::test_path("fixtures", "partial-dead-cluster.R")
  libs <- paste(.libPaths(), collapse = .Platform$path.sep)
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)),
    stdout = TRUE, stderr = TRUE,
    env = c("NOT_CRAN=true", paste0("R_LIBS=", shQuote(libs)),
            paste0("RNAPARALLEL_SRC=", shQuote(normalizePath(testthat::test_path("..", ".."))))),
    timeout = 30
  ))

  # The parent owns cleanup if the disposable R process dies before on.exit runs
  pid_lines <- grep("_PIDS=", out, value = TRUE)
  pids <- as.integer(unlist(strsplit(sub("^[^=]*=", "", pid_lines), ",", fixed = TRUE)))
  pids <- pids[!is.na(pids)]
  on.exit(for (pid in pids) {
    if (isTRUE(tools::pskill(pid, 0L))) tools::pskill(pid, tools::SIGKILL)
  }, add = TRUE)

  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  expect_identical(status, 0L, info = paste(out, collapse = "\n"))
  expect_true("PARTIAL_DEAD_CLUSTER_OK" %in% out, info = paste(out, collapse = "\n"))
})

test_that("a caller's own foreach backend returns chunks in dispatch order", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("doParallel")
  skip_if_not_installed("foreach")

  cl <- parallel::makeCluster(4L, type = "FORK")
  doParallel::registerDoParallel(cl)
  on.exit({
    foreach::registerDoSEQ()
    parallel::stopCluster(cl)
  }, add = TRUE)

  # A chunk count that is not a multiple of the backend width, because that is where any
  # regrouping shows up: reassembly is by tag, and a permutation that slipped through would
  # bind genes into another chunk's rows silently rather than fail.
  idx <- rnaparallel:::combat_row_chunks(21L, chunks = 7L)
  f <- function(i) cbind(i, i * 2L)
  expect_identical(
    rnaparallel:::combat_parallel_lapply(idx, f, workers = 2L, parallel_backend = "foreach",
                                         cells = Inf, min_cells = 0),
    lapply(idx, f))

  # a chunk that errors must arrive as a condition in ITS OWN slot, named by its own index
  parts <- rnaparallel:::combat_parallel_lapply(
    idx, function(i) {
      if (3L %in% i) stop("foreign foreach fixture boom")
      i
    }, workers = 2L, parallel_backend = "foreach", cells = Inf, min_cells = 0
  )
  boom <- which(vapply(idx, function(i) 3L %in% i, logical(1)))
  expect_true(inherits(parts[[boom]], "condition"))
  expect_false(any(vapply(parts[-boom], function(p) inherits(p, "condition"), logical(1))))
  expect_error(rnaparallel:::combat_parallel_check(parts, "foreign foreach", idx),
               "foreign foreach fixture boom")
})

test_that("a caller-registered foreach backend is not throttled to `workers`", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("doParallel")
  skip_if_not_installed("foreach")

  cl <- parallel::makeCluster(4L, type = "FORK")
  doParallel::registerDoParallel(cl)
  on.exit({
    foreach::registerDoSEQ()
    parallel::stopCluster(cl)
  }, add = TRUE)

  # 0.4.5 grouped chunks so a caller's backend could not exceed `workers`. It was aimed at a
  # nesting bug that never reproduced here, and it cost 2.08x on a wide registered backend,
  # held a whole group's results in each worker, and clamped a remote backend to this
  # machine's core count. A backend the caller registered is theirs; its width governs.
  idx <- rnaparallel:::combat_row_chunks(80L, chunks = 8L)
  out <- rnaparallel:::combat_parallel_lapply(
    idx, function(i) { Sys.sleep(0.05); Sys.getpid() }, workers = 2L,
    parallel_backend = "foreach", cells = Inf, min_cells = 0
  )
  expect_gt(length(unique(unlist(out))), 2L)
})

test_that("a caller's doRNG registration cannot move the master random stream", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("doRNG")
  skip_if_not_installed("doParallel")

  cl <- parallel::makeCluster(2L, type = "FORK")
  on.exit({ foreach::registerDoSEQ(); parallel::stopCluster(cl) }, add = TRUE)
  doParallel::registerDoParallel(cl)
  doRNG::registerDoRNG(42L)

  idx <- rnaparallel:::combat_row_chunks(12L, chunks = 4L)
  set.seed(7); expected <- sample(1e6L, 2L)
  # ComBat-seq's own sample() is on the serial side, so a stream the dispatch moved is a
  # different answer from the vendor under shrink = TRUE rather than a slower one
  set.seed(7)
  rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), workers = 2L,
                                       parallel_backend = "foreach",
                                       cells = Inf, min_cells = 0)
  expect_identical(sample(1e6L, 2L), expected)
})

test_that("the foreach cap does not group mclapply chunks", {
  skip_on_cran()
  skip_on_os("windows")

  idx <- rnaparallel:::combat_row_chunks(80L, chunks = 8L)
  out <- rnaparallel:::combat_parallel_lapply(
    idx, function(i) Sys.getpid(), workers = 2L, parallel_backend = "mclapply",
    cells = Inf, min_cells = 0
  )
  expect_identical(length(unique(unlist(out))), length(idx))
})

test_that("the mclapply backend works with one core", {
  skip_on_os("windows")

  idx <- rnaparallel:::combat_row_chunks(12L, chunks = 3L)
  out <- rnaparallel:::combat_parallel_lapply(
    idx, function(i) list(pid = Sys.getpid(), value = sum(i)), workers = 1L,
    parallel_backend = "mclapply", cells = Inf, min_cells = 0
  )

  expect_identical(vapply(out, `[[`, integer(1), "pid"), rep(Sys.getpid(), length(idx)))
  expect_identical(vapply(out, `[[`, integer(1), "value"), vapply(idx, sum, integer(1)))

  # Assert the PACKAGE's shortcut, not mclapply's. mc.cores = 1 runs in the master anyway, so
  # the two assertions above pass with the workers <= 1 shortcut deleted outright; only
  # watching whether mclapply is entered at all can tell the difference.
  calls <- 0L
  suppressMessages(trace("mclapply", where = asNamespace("parallel"), print = FALSE,
                         tracer = quote(NULL), exit = quote(NULL)))
  on.exit(suppressMessages(untrace("mclapply", where = asNamespace("parallel"))), add = TRUE)
  suppressMessages(trace("mclapply", where = asNamespace("parallel"), print = FALSE,
                         tracer = substitute(calls <<- calls + 1L)))
  rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), workers = 1L,
                                       parallel_backend = "mclapply",
                                       cells = Inf, min_cells = 0)
  expect_identical(calls, 0L)
})

test_that("the foreach backend still agrees after the caching change", {
  skip_if_not_installed("doParallel")
  vendor <- backend_fn()
  d <- make_counts(40, G = 200, n_per_batch = c(6, 6))
  ref <- quietly(vendor(d$counts, d$batch, group = NULL))
  got <- quietly(ComBat_seq_parallel(d$counts, d$batch, group = NULL,
                                     workers = 2L, parallel_backend = "foreach"))
  expect_identical(got, ref)
  rnaparallel::combat_cluster_stop()
})

test_that("a backend that reorders or duplicates chunks cannot scramble the result", {
  # Checking only the LENGTH of what a backend returns cannot tell a correct answer from
  # the same chunks in the wrong order, and binding a reordered list scrambles genes with
  # no error at all. Every job carries its chunk number so this is detectable.
  idx <- rnaparallel:::combat_row_chunks(20L, chunks = 4L)
  f <- function(i) sum(i)
  correct <- lapply(idx, f)

  # right length, wrong order: silently wrong before the tag existed
  reversed <- rnaparallel:::combat_parallel_lapply(
    idx, f, workers = 2L, parallel_backend = function(idx, f, w) rev(lapply(idx, f)))
  expect_identical(reversed, correct)

  # right length, one chunk delivered twice and another lost
  expect_error(
    rnaparallel:::combat_parallel_lapply(
      idx, f, workers = 2L,
      parallel_backend = function(idx, f, w) { o <- lapply(idx, f); o[[2]] <- o[[1]]; o }),
    "duplicated or out-of-range")
})

test_that("a chunk that comes back the wrong size is refused, not bound", {
  # the shape that introduced NA rows during the bind rather than failing
  idx <- rnaparallel:::combat_row_chunks(20L, chunks = 4L)
  parts <- lapply(idx, function(i) matrix(0, nrow = length(i), ncol = 2))
  expect_silent(rnaparallel:::combat_parallel_check(parts, "unit", idx))

  parts[[2]] <- parts[[2]][-1, , drop = FALSE]
  expect_error(rnaparallel:::combat_parallel_check(parts, "unit", idx),
               "came back with")
})

test_that("an uninstalled named backend is refused at any dispatch size", {
  # the dependency check used to sit below the size gate, so an unavailable framework
  # succeeded quietly on a small dispatch and errored on an otherwise identical large one
  skip_if(requireNamespace("future.apply", quietly = TRUE),
          "future.apply is installed, so this path cannot be exercised here")
  idx <- rnaparallel:::combat_row_chunks(20L, chunks = 4L)
  expect_error(
    rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), workers = 2L,
                                               parallel_backend = "future", cells = 10),
    "needs future.apply")
})
test_that("row chunks cover every gene exactly once", {
  for (ntag in c(1L, 2L, 7L, 100L)) {
    for (nch in c(1L, 2L, 3L, 8L, 1000L)) {
      idx <- rnaparallel:::combat_row_chunks(ntag, chunks = nch)
      expect_identical(sort(unlist(idx)), seq_len(ntag))
      expect_lte(length(idx), ntag)
    }
  }
})

test_that("chunks are clamped, so a huge value cannot fork once per gene", {
  # the bug this guards: chunks = nrow(counts) once forked 5000 processes, one per row
  expect_length(rnaparallel:::combat_row_chunks(10L, chunks = 100000L), 10L)
  expect_length(rnaparallel:::combat_row_chunks(10L, chunks = 0L), 1L)
  expect_length(rnaparallel:::combat_row_chunks(10L, chunks = -5L), 1L)
  expect_length(rnaparallel:::combat_row_chunks(10L, chunks = NA_integer_), 1L)
})

test_that("worker and chunk layouts do not change the answer", {
  d <- make_counts(11, G = 250, n_per_batch = c(6, 6))
  ref <- quietly(ComBat_seq_parallel(d$counts, d$batch, group = NULL, workers = 1L, chunks = 1L))
  for (layout in list(c(1, 8), c(2, 2), c(4, 1), c(4, 3), c(2, 64))) {
    got <- quietly(ComBat_seq_parallel(d$counts, d$batch, group = NULL,
                                  workers = layout[1], chunks = layout[2]))
    expect_identical(got, ref,
                     info = sprintf("workers=%s chunks=%s", layout[1], layout[2]))
  }
})

test_that("a single batch fails the same way the vendor fails", {
  # Not a companion defect. sva's ComBat_seq builds model.matrix(~ -1 + batch), which
  # cannot take a one-level factor, so the vendor errors here too. Being a drop-in
  # means reproducing the failure, not papering over it.
  vendor <- backend_fn()
  d <- make_counts(12, G = 60, n_per_batch = c(6, 6), adversarial = FALSE)
  one_batch <- rep("only", ncol(d$counts))

  ref <- tryCatch(quietly(vendor(d$counts, one_batch, group = NULL)),
                  error = function(e) conditionMessage(e))
  par <- tryCatch(quietly(ComBat_seq_parallel(d$counts, one_batch, group = NULL, workers = 2L)),
                  error = function(e) conditionMessage(e))
  expect_true(is.character(ref))
  expect_identical(par, ref)
})

test_that("workers must be a positive integer", {
  d <- make_counts(13, G = 50, n_per_batch = c(4, 4))
  expect_error(ComBat_seq_parallel(d$counts, d$batch, workers = 0L), "positive integer")
  expect_error(ComBat_seq_parallel(d$counts, d$batch, workers = NA_integer_), "positive integer")
})

test_that("interleaved chunks balance a matrix whose genes are ordered", {
  # The reason for round-robin: a filtered or sorted count matrix puts cheap genes
  # together and expensive ones together, so contiguous blocks hand one worker all
  # the work. Measured 313x imbalance on a matrix sorted by expression.
  ntag <- 400L
  cost <- seq_len(ntag)                      # stand-in for per-gene cost, monotone
  inter <- rnaparallel:::combat_row_chunks(ntag, chunks = 4L, interleave = TRUE)
  contig <- rnaparallel:::combat_row_chunks(ntag, chunks = 4L, interleave = FALSE)

  bal <- function(idx) { w <- vapply(idx, function(i) sum(cost[i]), 0); max(w) / min(w) }
  expect_lt(bal(inter), 1.05)                # essentially even
  expect_gt(bal(contig), 3)                  # badly skewed on ordered input
  # both still cover every row exactly once
  for (idx in list(inter, contig)) expect_identical(sort(unlist(idx)), seq_len(ntag))
})

test_that("interleaving does not permute the output", {
  # The correctness risk the reorder exists to remove: interleaved chunks come back
  # in chunk order, so without combat_row_order the genes would be silently scrambled.
  skip_if_not_installed("sva")
  d <- make_counts(50, G = 300, n_per_batch = c(6, 6))
  ref <- quietly(sva::ComBat_seq(d$counts, d$batch, group = NULL))
  for (ch in c(1L, 3L, 4L, 7L, 16L)) {
    got <- quietly(ComBat_seq_parallel(d$counts, d$batch, group = NULL,
                                       workers = 2L, chunks = ch))
    expect_identical(got, ref, info = paste("chunks =", ch))
    expect_identical(rownames(got), rownames(d$counts))
  }
})

test_that("combat_row_order is identity for a contiguous split", {
  idx <- rnaparallel:::combat_row_chunks(100L, chunks = 5L, interleave = FALSE)
  expect_identical(rnaparallel:::combat_row_order(idx), seq_len(100L))
})
## Exercising the real out-of-memory path needs a matrix large enough to be
## slow, machine-dependent and rude to CI. The condition it produces is a
## NULL element in the mclapply result, so the branch is driven directly instead.

test_that("a killed worker reports the death, not a conditionMessage failure", {
  # mclapply returns NULL for a child the kernel killed. The old code did
  # conditionMessage(attr(NULL, "condition")) and threw its own dispatch error,
  # burying the real OOM cause.
  parts <- list(matrix(1), NULL, matrix(3))
  expect_error(rnaparallel:::combat_parallel_check(parts, "unit"),
               "worker process died")
  expect_error(rnaparallel:::combat_parallel_check(parts, "unit"), "unit")
})

test_that("the killed-worker message never mentions conditionMessage dispatch", {
  parts <- list(NULL)
  msg <- tryCatch(rnaparallel:::combat_parallel_check(parts, "unit"),
                  error = function(e) conditionMessage(e))
  expect_false(grepl("applicable method", msg))
  expect_false(grepl("conditionMessage", msg))
  expect_match(msg, "worker process died")
})

test_that("a thrown error surfaces its own message", {
  err <- try(stop("boom in the worker"), silent = TRUE)
  parts <- list(matrix(1), err)
  msg <- tryCatch(rnaparallel:::combat_parallel_check(parts, "unit"),
                  error = function(e) conditionMessage(e))
  expect_match(msg, "boom in the worker")
})

test_that("both failure modes at once are both reported", {
  err <- try(stop("boom"), silent = TRUE)
  parts <- list(NULL, err)
  msg <- tryCatch(rnaparallel:::combat_parallel_check(parts, "unit"),
                  error = function(e) conditionMessage(e))
  expect_match(msg, "worker process died")
  expect_match(msg, "boom")
})

test_that("a clean result passes straight through unchanged", {
  parts <- list(matrix(1), matrix(2))
  expect_identical(rnaparallel:::combat_parallel_check(parts, "unit"), parts)
})

test_that("a killed worker survives the untagging layer instead of vanishing", {
  # The regression that mattered: `res[[i]] <- NULL` DELETES a list element rather than
  # storing NULL, so a dead chunk disappeared, every later chunk shifted up one slot, and
  # genes were bound into the wrong rows with combat_parallel_check reporting success.
  # The pre-existing test built `parts` by hand and called the check directly, so it never
  # went through the layer where the NULL was eaten.
  idx <- rnaparallel:::combat_row_chunks(100L, chunks = 4L)
  f <- function(i) matrix(0, nrow = length(i), ncol = 2)

  dead <- function(idx, f, w) { o <- lapply(idx, f); o[2] <- list(NULL); o }
  out <- rnaparallel:::combat_parallel_lapply(idx, f, workers = 2L,
                                                    parallel_backend = dead)
  expect_length(out, length(idx))
  expect_true(is.null(out[[2]]))
  expect_error(rnaparallel:::combat_parallel_check(out, "unit", idx),
               "worker process died")
})

test_that("a result count that does not match the chunk count is refused", {
  # `got != want` recycles when the lengths differ, so three results against four chunks
  # compared all-FALSE and passed with only a warning.
  idx <- rnaparallel:::combat_row_chunks(100L, chunks = 4L)
  parts <- lapply(idx, function(i) matrix(0, nrow = length(i), ncol = 2))
  expect_error(rnaparallel:::combat_parallel_check(parts[-1], "unit", idx),
               "result\\(s\\) for 4 chunk\\(s\\)")
})

test_that("a GLM chunk of the wrong length is refused after the bind", {
  # A list-shaped chunk cannot be row-checked before binding, so an over-long one produced
  # a correctly SHAPED result with half the genes carrying another gene's coefficient.
  skip_if_not_installed("edgeR")
  set.seed(51)
  y <- matrix(rnbinom(20 * 8, mu = 200, size = 4), 20, 8)
  # A continuous column, so this design does not reach the one-group kernel and the row split
  # is actually taken. cbind(1, rep(0:1, each = 4)) is one-group and now returns whole.
  design <- cbind(1, rep(0:1, each = 4), rnorm(8))
  off <- matrix(log(colSums(y)), 20, 8, byrow = TRUE)

  liar <- function(idx, f, w) {
    o <- lapply(idx, f)
    # executors receive the tagged wrapper, so the fit is at $value
    v <- o[[2]]$value
    v$coefficients <- rbind(v$coefficients, v$coefficients[1, , drop = FALSE])
    o[[2]]$value <- v
    o
  }
  expect_error(
    rnaparallel:::glmFit_rows_parallel(y, design = design, dispersion = 0.1,
                                             offset = off, workers = 2L, chunks = 4L,
                                             parallel_backend = liar),
    "bound to")
})

test_that("glmFit_rows_parallel refuses to split without an explicit offset", {
  # each worker would rebuild library sizes from its own genes: measured 1.8 of coefficient
  # movement, returned silently
  skip_if_not_installed("edgeR")
  set.seed(52)
  y <- matrix(rnbinom(20 * 8, mu = 200, size = 4), 20, 8)
  design <- cbind(1, rep(0:1, each = 4))
  expect_error(
    rnaparallel:::glmFit_rows_parallel(y, design = design, dispersion = 0.1,
                                             offset = NULL, workers = 2L, chunks = 4L,
                                             parallel_backend = "serial"),
    "explicit offset")
})
