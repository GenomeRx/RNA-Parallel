rp_os_children <- function() {
  kids <- tryCatch(ps::ps_children(ps::ps_handle()), error = function(e) list())
  vapply(kids, ps::ps_pid, integer(1))
}

rp_pid_state <- function(pid) {
  s <- suppressWarnings(system2("ps", c("-o", "stat=", "-p", pid), stdout = TRUE, stderr = FALSE))
  s <- trimws(s[nzchar(trimws(s))])
  if (length(s)) s[1L] else ""
}

rp_wait_gone <- function(pid, secs = 3) {
  deadline <- Sys.time() + secs
  while (nzchar(rp_pid_state(pid)) && Sys.time() < deadline) Sys.sleep(0.05)
  rp_pid_state(pid)
}

rp_subprocess_libs <- function() {
  paste0("R_LIBS=", shQuote(paste(.libPaths(), collapse = .Platform$path.sep)))
}


# ---- the progress reporter -------------------------------------------------------

test_that("the BiocParallel backend returns in a fresh session instead of hanging on the reporter", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("BiocParallel")
  skip_if_not_installed("sva")

  script <- tempfile(fileext = ".R")
  writeLines(c(
    rp_load_line(),
    "utils::assignInNamespace(\"rp_reporter_visible\", function() TRUE, \"rnaparallel\")",
    "options(combat.min.cells = 0, combat.min.disp.cells = 0, combat.min.batch.cells = 0,",
    "        combat.min.glm.cells = 0, combat.mem.guard = FALSE)",
    "set.seed(1)",
    "m <- matrix(rnbinom(400 * 12, mu = 50, size = 5), 400, 12)",
    "batch <- factor(rep(1:2, each = 6))",
    "invisible(capture.output(x <- ComBat_seq_parallel(m, batch, group = NULL, workers = 2L,",
    "                                                  parallel_backend = \"BiocParallel\")))",
    "invisible(capture.output(ref <- sva::ComBat_seq(m, batch, group = NULL)))",
    "cat(\"\\nBIOC_OK\", identical(x, ref), \"\\n\")"), script)
  t0 <- proc.time()[["elapsed"]]
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)),
    stdout = TRUE, stderr = TRUE, env = c("NOT_CRAN=true", rp_subprocess_libs()), timeout = 90))
  secs <- proc.time()[["elapsed"]] - t0
  expect_true(any(grepl("BIOC_OK TRUE", out, fixed = TRUE)),
              info = sprintf("%.1f s: %s", secs, paste(tail(out, 5), collapse = "\n")))
  expect_lt(secs, 60)
})

test_that("a BiocParallel dispatch does not collect the caller's own mcparallel job", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("BiocParallel")
  withr::local_options(combat.progress = FALSE)

  job <- parallel::mcparallel({ Sys.sleep(1); "user result" })
  idx <- rnaparallel:::combat_row_chunks(12L, chunks = 2L)
  out <- rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), workers = 2L,
                                              parallel_backend = "BiocParallel",
                                              cells = Inf, min_cells = 0)
  expect_identical(out, lapply(idx, sum))
  got <- parallel::mccollect(job, wait = TRUE)
  expect_identical(got[[1L]], "user result")
})

test_that("the progress reporter is invisible to parallel::children() and leaves no zombie", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("ps")
  withr::local_options(combat.progress = TRUE, combat.progress.dir = withr::local_tempdir())
  local_mocked_bindings(rp_reporter_visible = function() TRUE, .package = "rnaparallel")

  children <- utils::getFromNamespace("children", "parallel")
  spare <- vapply(children(), function(p) as.integer(p$pid), integer(1))
  on.exit(rnaparallel:::combat_reap(spare), add = TRUE)
  before_os <- rp_os_children()
  before_mc <- length(children())
  seen <- NULL
  exec <- function(idx, f, w) {
    seen <<- list(reporter = isTRUE(rnaparallel:::.rp_dispatch$reporter),
                  mc = length(children()), os = setdiff(rp_os_children(), before_os))
    lapply(idx, f)
  }
  idx <- rnaparallel:::combat_row_chunks(12L, chunks = 3L)
  out <- rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), workers = 2L,
                                              parallel_backend = exec, cells = Inf, min_cells = 0)
  expect_identical(out, lapply(idx, sum))
  expect_true(seen$reporter)
  expect_length(seen$os, 1L)
  expect_identical(seen$mc, before_mc)
  expect_identical(length(children()), before_mc)
  if (length(seen$os) == 1L) expect_identical(rp_wait_gone(seen$os), "")
})

test_that("a future multicore plan is not throttled by the reporter", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("future.apply")
  withr::local_options(combat.progress = TRUE, combat.progress.dir = withr::local_tempdir())
  local_mocked_bindings(rp_reporter_visible = function() TRUE, .package = "rnaparallel")
  old <- future::plan(future::multicore, workers = 2)
  on.exit(future::plan(old), add = TRUE)

  ws <- character()
  idx <- rnaparallel:::combat_row_chunks(12L, chunks = 2L)
  out <- withCallingHandlers(
    rnaparallel:::combat_parallel_lapply(idx, function(i) { Sys.sleep(0.3); sum(i) },
                                         workers = 2L, parallel_backend = "future",
                                         cells = Inf, min_cells = 0),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_identical(out, lapply(idx, function(i) sum(i)))
  expect_false(any(grepl("active multicore processes", ws)), info = paste(ws, collapse = "\n"))
})

test_that("combat.progress = FALSE starts no reporter", {
  skip_on_os("windows")
  skip_if_not_installed("ps")
  withr::local_options(combat.progress = FALSE, combat.progress.dir = withr::local_tempdir())
  local_mocked_bindings(rp_reporter_visible = function() TRUE, .package = "rnaparallel")

  before_os <- rp_os_children()
  seen <- NULL
  exec <- function(idx, f, w) {
    seen <<- list(reporter = isTRUE(rnaparallel:::.rp_dispatch$reporter),
                  os = setdiff(rp_os_children(), before_os))
    lapply(idx, f)
  }
  idx <- rnaparallel:::combat_row_chunks(12L, chunks = 3L)
  rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), workers = 2L,
                                       parallel_backend = exec, cells = Inf, min_cells = 0)
  expect_false(seen$reporter)
  expect_length(seen$os, 0L)
})

test_that("serial and one-worker dispatches fork no reporter", {
  skip_on_os("windows")
  skip_if_not_installed("ps")
  withr::local_options(combat.progress = TRUE, combat.progress.dir = withr::local_tempdir())
  local_mocked_bindings(rp_reporter_visible = function() TRUE, .package = "rnaparallel")

  idx <- rnaparallel:::combat_row_chunks(12L, chunks = 3L)
  for (case in list(list(be = "serial", w = 2L), list(be = "mclapply", w = 1L))) {
    before_os <- rp_os_children()
    seen <- integer()
    f <- function(i) { seen <<- c(seen, length(setdiff(rp_os_children(), before_os))); sum(i) }
    out <- rnaparallel:::combat_parallel_lapply(idx, f, workers = case$w,
                                                parallel_backend = case$be,
                                                cells = Inf, min_cells = 0)
    expect_identical(out, lapply(idx, sum), info = case$be)
    expect_identical(seen, rep(0L, length(idx)), info = case$be)
  }
})

test_that("mclapply inside a caller's own fork is counted serial and forks no reporter", {
  skip_on_os("windows")
  skip_if_not_installed("ps")
  withr::local_options(combat.progress = TRUE, combat.progress.dir = withr::local_tempdir())
  local_mocked_bindings(rp_reporter_visible = function() TRUE, .package = "rnaparallel")

  job <- parallel::mcparallel({
    rnaparallel:::rp_count_reset()
    base <- rp_os_children()
    kids <- integer()
    idx <- rnaparallel:::combat_row_chunks(12L, chunks = 3L)
    out <- rnaparallel:::combat_parallel_lapply(
      idx, function(i) { kids <<- c(kids, length(setdiff(rp_os_children(), base))); sum(i) },
      workers = 2L, parallel_backend = "mclapply", cells = Inf, min_cells = 0)
    list(ok = identical(out, lapply(idx, sum)), par = rnaparallel:::.rp_dispatch$par,
         ser = rnaparallel:::.rp_dispatch$ser, kids = kids)
  })
  res <- parallel::mccollect(job)[[1L]]
  expect_true(is.list(res), info = paste(format(res), collapse = " "))
  expect_true(res$ok)
  expect_identical(res$par, 0L)
  expect_identical(res$ser, 1L)
  expect_true(all(res$kids == 0L))
})

test_that("the reporter exits when its master is killed instead of pinning the master's pages", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("ps")

  d <- withr::local_tempdir()
  script <- file.path(d, "master.R")
  writeLines(c(
    rp_load_line(),
    "utils::assignInNamespace(\"rp_reporter_visible\", function() TRUE, \"rnaparallel\")",
    sprintf("d <- %s", deparse(d)),
    "dir.create(file.path(d, 'progress'))",
    "options(combat.progress = TRUE, combat.progress.dir = file.path(d, 'progress'))",
    "exec <- function(idx, f, w) {",
    "  Sys.sleep(1)",
    "  kids <- vapply(ps::ps_children(ps::ps_handle()), ps::ps_pid, integer(1))",
    "  writeLines(as.character(kids), file.path(d, 'kids'))",
    "  tools::pskill(Sys.getpid(), tools::SIGKILL)",
    "}",
    "rnaparallel:::combat_parallel_lapply(list(1:3, 4:6), function(i) sum(i), 2L, exec,",
    "                                     cells = Inf, min_cells = 0)"), script)
  log <- file.path(d, "master.log")
  suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)),
                           stdout = log, stderr = log,
                           env = c("NOT_CRAN=true", rp_subprocess_libs()), timeout = 60))
  kids <- if (file.exists(file.path(d, "kids"))) as.integer(readLines(file.path(d, "kids"))) else integer()
  on.exit(for (k in kids) try(tools::pskill(k, tools::SIGKILL), silent = TRUE), add = TRUE)
  expect_length(kids, 1L)
  for (k in kids) expect_identical(rp_wait_gone(k, 4), "", info = paste("reporter", k))
})


# ---- the orphan check ------------------------------------------------------------

test_that("a socket worker that cannot see the master's pid is not killed by the orphan check", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("ps")
  withr::local_options(combat.progress = FALSE)

  gone <- parallel::mcparallel(NULL)
  parallel::mccollect(gone)
  cl <- parallel::makePSOCKcluster(1L)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  rp_cluster_load(cl)
# A pid that names nothing stands in for a worker on another host or in another PID namespace.
  exec <- function(idx, f, w) {
    environment(f)$master_pid <- gone$pid
    parallel::clusterApply(cl, idx, f)
  }
  f <- function(i) sum(i)
  environment(f) <- baseenv()
  idx <- rnaparallel:::combat_row_chunks(12L, chunks = 3L)
  out <- tryCatch(rnaparallel:::combat_parallel_lapply(idx, f, workers = 2L, parallel_backend = exec,
                                                       cells = Inf, min_cells = 0),
                  error = function(e) conditionMessage(e))
  expect_identical(out, lapply(idx, sum))
})

test_that("a fork whose master died is stopped even while the master is an unreaped zombie", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("ps")

  d <- withr::local_tempdir()
  script <- file.path(d, "master.R")
  writeLines(c(
    rp_load_line(),
    "options(combat.progress = FALSE)",
    sprintf("d <- %s", deparse(d)),
    "writeLines(as.character(Sys.getpid()), file.path(d, 'master.pid'))",
    "exec <- function(idx, f, w) {",
    "  w <- parallel::mcparallel({ Sys.sleep(2); f(idx[[1L]]); writeLines('ran', file.path(d, 'marker')) },",
    "                            detached = TRUE)",
    "  writeLines(as.character(w$pid), file.path(d, 'worker.pid'))",
    "  tools::pskill(Sys.getpid(), tools::SIGKILL)",
    "}",
    "rnaparallel:::combat_parallel_lapply(list(1:3, 4:6), function(i) sum(i), 2L, exec,",
    "                                     cells = Inf, min_cells = 0)"), script)
# sh execs a sleep that never reaps, so the master stays a zombie after it exits.
  cmd <- sprintf("echo $$ > %s; NOT_CRAN=true %s %s --vanilla %s > %s 2>&1 & exec sleep 30",
                 shQuote(file.path(d, "sh.pid")), rp_subprocess_libs(),
                 shQuote(file.path(R.home("bin"), "Rscript")), shQuote(script),
                 shQuote(file.path(d, "master.log")))
  system2("sh", c("-c", shQuote(cmd)), wait = FALSE)
  rd <- function(nm) {
    p <- file.path(d, nm)
    if (file.exists(p)) suppressWarnings(as.integer(readLines(p, warn = FALSE)[1L])) else NA_integer_
  }
  deadline <- Sys.time() + 40
  while (is.na(rd("worker.pid")) && Sys.time() < deadline) Sys.sleep(0.1)
  on.exit({
    for (p in c(rd("worker.pid"), rd("sh.pid"))) if (!is.na(p)) try(tools::pskill(p, tools::SIGKILL), silent = TRUE)
  }, add = TRUE)
  skip_if(is.na(rd("worker.pid")), "the master never dispatched")
  Sys.sleep(0.5)
  skip_if_not(startsWith(rp_pid_state(rd("master.pid")), "Z"), "the master was reaped, not left a zombie")
  Sys.sleep(3)
  expect_false(file.exists(file.path(d, "marker")))
})


# ---- memory readers --------------------------------------------------------------

test_that("the memory readers return real numbers off Linux through ps", {
  skip_if(file.exists("/proc/meminfo"), "this machine has /proc")
  skip_on_os("windows")
  skip_if_not_installed("ps")
  rss <- rnaparallel:::rp_mem_rss()
  peak <- rnaparallel:::rp_mem_peak()
  total <- rnaparallel:::rp_mem_total()
  expect_true(is.numeric(rss) && !is.na(rss) && rss > 0)
  expect_true(is.numeric(peak) && !is.na(peak) && peak >= 0.9 * rss)
  expect_true(is.numeric(total) && !is.na(total) && total > rss)
})

test_that("rp_mem_rss reads VmRSS in kB, independent of page size", {
  f <- rnaparallel:::rp_mem_rss
  e <- new.env(parent = environment(f))
  e$file.exists <- function(x) identical(x, "/proc/self/status") || base::file.exists(x)
  e$readLines <- function(con, ...) {
    if (!identical(con, "/proc/self/status")) return(base::readLines(con, ...))
    c("Name:\tR", "VmPeak:\t  999999 kB", "VmRSS:\t  123456 kB", "RssAnon:\t  4 kB")
  }
  environment(f) <- e
  expect_identical(f(), 123456 * 1024)
})

test_that("rp_mem_rss agrees with ps's page-count reading on Linux", {
  skip_if_not(file.exists("/proc/self/status"), "no /proc on this machine")
  skip_if_not_installed("ps")
  rss <- rnaparallel:::rp_mem_rss()
  ref <- ps::ps_memory_info()[["rss"]]
  expect_true(abs(rss - ref) / ref < 0.2)
})


# ---- the caller's random stream --------------------------------------------------

test_that("a seeding custom executor cannot move the caller's random stream", {
  skip_if_not_installed("future.apply")
  skip_if_not_installed("sva")
  withr::local_options(combat.progress = FALSE)
  withr::local_preserve_seed()
  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)
  seeded <- function(idx, f, w) future.apply::future_lapply(idx, f, future.seed = TRUE)
  idx <- rnaparallel:::combat_row_chunks(12L, chunks = 4L)

  set.seed(7)
  s0 <- .Random.seed
  invisible(rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), 2L, seeded,
                                                 cells = Inf, min_cells = 0))
  expect_identical(.Random.seed, s0)

  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  invisible(rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), 2L, seeded,
                                                 cells = Inf, min_cells = 0))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))

  set.seed(3)
  G <- 300; n <- 12
  cts <- matrix(rnbinom(G * n, mu = 80, size = 4), G, n,
                dimnames = list(paste0("g", 1:G), paste0("s", 1:n)))
  batch <- rep(1:2, each = 6); grp <- rep(1:2, 6)
  set.seed(11)
  invisible(capture.output(ref <- sva::ComBat_seq(cts, batch, group = grp, shrink = TRUE,
                                                  gene.subset.n = 40)))
  after_ref <- sample(1e6L, 3L)
  set.seed(11)
  invisible(capture.output(got <- ComBat_seq_parallel(cts, batch, group = grp, shrink = TRUE,
                                                      gene.subset.n = 40, workers = 2L,
                                                      parallel_backend = seeded)))
  expect_identical(got, ref)
  expect_identical(sample(1e6L, 3L), after_ref)
})


# ---- drift gates -----------------------------------------------------------------

test_that("the ComBat-seq gate refuses a backend whose across-batch calls are no longer rebindable", {
  skip_if_not_installed("sva")
  expect_no_error(rnaparallel:::combat_backend(sva::ComBat_seq))
  edit <- function(from, to) {
    f <- sva::ComBat_seq
    txt <- deparse(body(f), width.cutoff = 500L)
    stopifnot(sum(grepl(from, txt, fixed = TRUE)) == 1L)
    body(f) <- str2lang(paste(sub(from, to, txt, fixed = TRUE), collapse = "\n"))
    f
  }
  common <- edit("disp_common <- sapply(", "disp_common <- base::sapply(")
  tagwise <- edit("genewise_disp_lst <- lapply(", "genewise_disp_lst <- base::lapply(")
  expect_error(rnaparallel:::combat_backend(common), "estimateGLMCommonDisp")
  expect_error(rnaparallel:::combat_backend(tagwise), "estimateGLMTagwiseDisp")

  mk <- function(apply_src) {
    f <- function(counts, batch, group = NULL, covar_mod = NULL, full_mod = TRUE,
                  shrink = FALSE, shrink.disp = FALSE, gene.subset.n = NULL) NULL
    body(f) <- str2lang(sprintf(
      "{ match_quantiles(a, b, c, d, e); glmFit(y); glmFit.default(y); estimateGLMTagwiseDisp(y); %s }",
      apply_src))
    e <- new.env(parent = globalenv())
    assign("match_quantiles",
           function(counts_sub, old_mu, old_phi, new_mu, new_phi) counts_sub, envir = e)
    environment(f) <- e
    f
  }
  tag <- "lapply(z, function(j) estimateGLMTagwiseDisp(j))"
  expect_silent(rnaparallel:::combat_backend(mk(paste(
    "n <- sapply(z, length); sapply(z, function(i) estimateGLMCommonDisp(i));", tag))))
  expect_silent(rnaparallel:::combat_backend(mk(paste(
    "g <- function(i) estimateGLMCommonDisp(i); sapply(X = z, FUN = g);", tag))))
  expect_error(rnaparallel:::combat_backend(mk(paste("sapply(z, length);", tag))),
               "sapply over estimateGLMCommonDisp")
  expect_error(rnaparallel:::combat_backend(mk(paste(
    "sapply(z, function(i) estimateGLMCommonDisp(i));",
    "lapply(z, function(j) { sample(3); estimateGLMTagwiseDisp(j) })"))),
    "lapply over estimateGLMTagwiseDisp")
})

# assignInNamespace replaces both the namespace binding and the registered S3 method, as a new edgeR would, and finds the generic in its caller's frame.
rp_swap_tagwise_default <- function(fn) {
  estimateGLMTagwiseDisp <- edgeR::estimateGLMTagwiseDisp
  utils::assignInNamespace("estimateGLMTagwiseDisp.default", fn, ns = "edgeR")
}

test_that("the tagwise split stands down when edgeR changes any default it transcribes", {
  withr::local_options(combat.progress = FALSE, combat.min.disp.cells = 0)
  ns <- asNamespace("edgeR")
  orig <- utils::getS3method("estimateGLMTagwiseDisp", "default", envir = ns)
  set.seed(5)
  y <- matrix(rnbinom(200 * 12, mu = 60, size = 4), 200, 12)
  design <- stats::model.matrix(~ stats::rnorm(12))
  calls <- 0L
  spy <- function(idx, f, workers) { calls <<- calls + 1L; lapply(idx, f) }
  run_split <- function() rnaparallel:::estimateGLMTagwiseDisp_rows_parallel(
    y, design = design, dispersion = 0.1, prior.df = 0, workers = 2L, chunks = 4L,
    parallel_backend = spy)

  rnaparallel:::rp_count_reset()
  got0 <- run_split()
  expect_identical(got0, edgeR::estimateGLMTagwiseDisp(y, design = design, dispersion = 0.1,
                                                       prior.df = 0))
  expect_length(rnaparallel:::.rp_dispatch$fallback, 0L)
  expect_identical(calls, 1L)

  txt <- paste(deparse(body(orig), width.cutoff = 500L), collapse = "\n")
  drifts <- c(
    "offset <- log(colSums(y))" = "offset <- log(colSums(y) * calcNormFactors(y))",
    "span <- (10/ntags)^0.23" = "span <- (10/ntags)^0.3",
    "AveLogCPM <- aveLogCPM(y, offset = offset, weights = weights)" =
      "AveLogCPM <- aveLogCPM(y, offset = offset, weights = weights, prior.count = 0.5)")
  withr::defer(rp_swap_tagwise_default(orig))
  for (from in names(drifts)) {
    stopifnot(grepl(from, txt, fixed = TRUE))
    edited <- orig
    body(edited) <- str2lang(sub(from, drifts[[from]], txt, fixed = TRUE))
    rp_swap_tagwise_default(edited)
    ref <- edgeR::estimateGLMTagwiseDisp(y, design = design, dispersion = 0.1, prior.df = 0)
    rnaparallel:::rp_count_reset()
    calls <- 0L
    expect_identical(run_split(), ref, info = from)
    expect_true("estimateGLMTagwiseDisp" %in% rnaparallel:::.rp_dispatch$fallback, info = from)
    expect_identical(calls, 0L, info = from)
  }
})


# ---- foreach state ---------------------------------------------------------------

rp_foreach_state <- function() {
  g <- utils::getFromNamespace(".foreachGlobals", "foreach")
  keys <- intersect(c("fun", "data", "info"), ls(g, all.names = TRUE))
  list(g = g, saved = mget(keys, envir = g))
}

rp_foreach_restore <- function(st) {
  rm(list = intersect(c("fun", "data", "info"), ls(st$g, all.names = TRUE)), envir = st$g)
  for (k in names(st$saved)) assign(k, st$saved[[k]], envir = st$g)
}

test_that("a foreach dispatch leaves an unregistered session unregistered", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("doParallel")
  withr::local_options(combat.progress = FALSE)
  st <- rp_foreach_state()
  withr::defer(rp_foreach_restore(st))
  rm(list = names(st$saved), envir = st$g)
  expect_false(foreach::getDoParRegistered())

  idx <- rnaparallel:::combat_row_chunks(12L, chunks = 3L)
  out <- rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), workers = 2L,
                                              parallel_backend = "foreach",
                                              cells = Inf, min_cells = 0)
  expect_identical(out, lapply(idx, sum))
  expect_false(foreach::getDoParRegistered())

  foreach::registerDoSEQ()
  rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), workers = 2L,
                                       parallel_backend = "foreach", cells = Inf, min_cells = 0)
  expect_identical(foreach::getDoParName(), "doSEQ")
})

test_that("the package's own foreach pool keeps no dispatch payload after it returns", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("doParallel")
  withr::local_options(combat.progress = FALSE)
  st <- rp_foreach_state()
  withr::defer(rp_foreach_restore(st))
  foreach::registerDoSEQ()

  idx <- rnaparallel:::combat_row_chunks(12L, chunks = 3L)
  out <- rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), workers = 2L,
                                              parallel_backend = "foreach",
                                              cells = Inf, min_cells = 0)
  expect_identical(out, lapply(idx, sum))
  cl <- rnaparallel:::.combat_clusters$FORK$cl
  held <- unlist(parallel::clusterCall(cl, function() {
    g <- get(".doSnowGlobals", envir = asNamespace("doParallel"))
    exists("exportenv", envir = g, inherits = FALSE)
  }))
  expect_false(any(held))
  rnaparallel::combat_cluster_stop()
})


# ---- stale progress rows ---------------------------------------------------------

test_that("the reporter reads only rows written after its dispatch began", {
  d <- withr::local_tempdir()
  writeLines("1\told stage\t1\tstart", file.path(d, "rnaparallel-1.tsv"))
  cache <- rnaparallel:::rp_progress_mark(d)
  cat("2\tnew\t1\tstart\n2\tnew\t1\tdone\n", file = file.path(d, "rnaparallel-1.tsv"), append = TRUE)
  cat("3\tnew\t2\tstart\n3\tnew\t2\tdone\n", file = file.path(d, "rnaparallel-2.tsv"))
  s <- rnaparallel:::rp_progress_summarise(rnaparallel:::rp_progress_read(d, cache = cache))
  expect_identical(c(s$done, s$started, s$stalled), c(2L, 2L, 0L))
  expect_identical(s$stage, "new")
})


# ---- input gates and memory limits -----------------------------------------------

test_that("the match_quantiles gate still refuses negative means and non-positive dispersions", {
  cs <- matrix(5L, 3, 2); om <- matrix(2, 3, 2)
  expect_null(rnaparallel:::combat_mq_dispatch(identity, cs, replace(om, 1L, -1), c(1, 1, 1)))
  expect_null(rnaparallel:::combat_mq_dispatch(identity, cs, om, c(1, 0, 1)))
  expect_null(rnaparallel:::combat_mq_dispatch(identity, cs, om, c(1, -2, 1)))
})

test_that("no chunk exceeds combat.mem.chunk.cells when the budget can hold min_rows", {
  for (case in list(c(10, 3, 10), c(100, 50, 500), c(1000, 7, 333), c(17, 4, 9), c(50, 1, 7))) {
    withr::local_options(combat.mem.chunk.cells = case[3])
    idx <- rnaparallel:::combat_row_chunks(case[1], workers = 2L, ncol = case[2])
    expect_lte(max(lengths(idx)) * case[2], case[3])
    expect_identical(sort(unlist(idx)), seq_len(case[1]))
  }
})

test_that("the R_MAX_VSIZE tier never rounds above the requested fraction", {
  testthat::local_mocked_bindings(rp_mem_total = function() 100 * 2^30, .package = "rnaparallel")
  p <- withr::local_tempfile()
  v <- suppressMessages(rnaparallel_set_mem_limit(fraction = 0.5, path = p))
  expect_identical(readLines(p), "R_MAX_VSIZE=32Gb")
  expect_identical(v, 32 * 2^30)

  testthat::local_mocked_bindings(rp_mem_total = function() 12 * 2^30, .package = "rnaparallel")
  p2 <- withr::local_tempfile()
  v2 <- suppressMessages(rnaparallel_set_mem_limit(fraction = 0.5, path = p2))
  expect_identical(readLines(p2), "R_MAX_VSIZE=6144Mb")
  expect_identical(v2, 6 * 2^30)
})
