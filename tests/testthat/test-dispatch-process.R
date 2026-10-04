## One dispatch as a process: interrupts, orphaned forks, parallel's RNG stream, the reporter and row-split stand-downs.

dp_rscript <- function(lines, timeout = 120) {
  script <- tempfile(fileext = ".R")
  writeLines(c(rp_load_line(), lines), script)
  libs <- paste0("R_LIBS=", shQuote(paste(.libPaths(), collapse = .Platform$path.sep)))
  suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)),
                           stdout = TRUE, stderr = TRUE, env = c("NOT_CRAN=true", libs),
                           timeout = timeout))
}

dp_alive <- function(pids) vapply(pids, function(p) isTRUE(tools::pskill(p, 0L)), logical(1))

dp_kill <- function(pids) for (p in pids) try(tools::pskill(p, tools::SIGKILL), silent = TRUE)

dp_conds <- function(expr) {
  w <- character()
  value <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}

dp_backends <- function() {
  b <- list(serial = "serial", custom = function(idx, f, workers) lapply(idx, f))
  if (.Platform$OS.type != "windows") b$mclapply <- "mclapply"
  b
}

# assignInNamespace also rewrites a registered S3 method, and finds the method's generic in its caller's frame.
dp_swap_edger <- function(name, fn) {
  glmFit <- edgeR::glmFit
  estimateGLMTagwiseDisp <- edgeR::estimateGLMTagwiseDisp
  utils::assignInNamespace(name, fn, ns = "edgeR")
}

# The chunks sleep 8 s; a SIGINT sent once both have started must reach the caller well before they finish.
dp_interrupt_script <- function(backend, work) c(
  "options(combat.progress = FALSE)",
  sprintf("work <- %s", deparse(work)),
  "f <- function(ii) { writeLines('', file.path(work, paste0('pid-', Sys.getpid()))); Sys.sleep(8); ii }",
  "environment(f) <- list2env(list(work = work), parent = baseenv())",
  "master <- Sys.getpid()",
  "killer <- parallel::mcparallel({",
  "  t <- Sys.time()",
  "  while (length(list.files(work, '^pid-')) < 2L && Sys.time() - t < 60) Sys.sleep(0.05)",
  "  Sys.sleep(0.5)",
  "  writeLines(format(as.numeric(Sys.time()), digits = 15), file.path(work, 'sent'))",
  "  tools::pskill(master, tools::SIGINT)",
  "}, detached = TRUE)",
  sprintf(paste0("r <- tryCatch(rnaparallel:::combat_parallel_lapply(list(1:2, 3:4), f, 2L, %s, ",
                 "cells = Inf, min_cells = 0), interrupt = function(i) 'INTERRUPTED')"),
          deparse(backend)),
  "back <- as.numeric(Sys.time())",
  "sent <- as.numeric(readLines(file.path(work, 'sent')))",
  "cat(sprintf('\\nRESULT %s %.2f\\n', if (is.character(r)) r else 'value', back - sent))",
  "g <- function(ii) ii * 10L",
  "environment(g) <- baseenv()",
  sprintf("r2 <- rnaparallel:::combat_parallel_lapply(list(1:2, 3:4, 5:6), g, 2L, %s, cells = Inf, min_cells = 0)",
          deparse(backend)),
  "cat('\\nSECOND', identical(r2, list(1:2 * 10L, 3:4 * 10L, 5:6 * 10L)), '\\n')",
  "rnaparallel::combat_cluster_stop()")

dp_interrupt_check <- function(backend) {
  work <- withr::local_tempdir(.local_envir = parent.frame())
  withr::defer(dp_kill(as.integer(sub("^pid-", "", list.files(work, "^pid-")))),
               envir = parent.frame())
  out <- dp_rscript(dp_interrupt_script(backend, work))
  res <- strsplit(grep("^RESULT ", out, value = TRUE), " ", fixed = TRUE)
  info <- paste(tail(out, 8), collapse = "\n")
  expect_length(list.files(work, "^pid-"), 2L)
  expect_length(res, 1L)
  if (length(res) == 1L) {
    expect_identical(res[[1L]][2L], "INTERRUPTED", info = info)
    expect_lt(as.numeric(res[[1L]][3L]), 4)
  }
  expect_true(any(grepl("^SECOND TRUE", out)), info = info)
}


# ---- interrupts ------------------------------------------------------------------------

test_that("an interrupt reaches the caller during a foreach dispatch on the package's own pool", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("doParallel")
  dp_interrupt_check("foreach")
})

test_that("an interrupt reaches the caller during a BiocParallel dispatch with no attached children", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("BiocParallel")
  dp_interrupt_check("BiocParallel")
})


# ---- orphaned forks --------------------------------------------------------------------

test_that("the master loads ps itself before forking, so no fork has to", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("ps")
  skip_if(exists("Sys.getppid", envir = baseenv(), mode = "function"), "base R reads the ppid")
  out <- dp_rscript(c(
    "options(combat.progress = FALSE)",
    "cat('\\nBEFORE', isNamespaceLoaded('ps'), '\\n')",
    "idx <- rnaparallel:::combat_row_chunks(12L, chunks = 2L)",
    "r <- rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), 2L, 'mclapply', cells = Inf, min_cells = 0)",
    "cat('\\nAFTER', isNamespaceLoaded('ps'), identical(r, lapply(idx, sum)), '\\n')"))
  info <- paste(tail(out, 5), collapse = "\n")
  expect_true(any(grepl("^BEFORE FALSE", out)), info = info)
  expect_true(any(grepl("^AFTER TRUE TRUE", out)), info = info)
})

test_that("a fork orphaned by its master exits without deleting the master's tempdir", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if(!requireNamespace("ps", quietly = TRUE) &&
            !exists("Sys.getppid", envir = baseenv(), mode = "function"),
          "nothing on this build reads the ppid")
  work <- withr::local_tempdir()
  script <- tempfile(fileext = ".R")
  writeLines(c(
    rp_load_line(),
    "options(combat.progress = FALSE)",
    sprintf("work <- %s", deparse(work)),
    "writeLines(tempdir(), file.path(work, 'tempdir'))",
    "writeLines(as.character(Sys.getpid()), file.path(work, 'master'))",
    "f <- function(ii) { writeLines('', file.path(work, paste0('chunk-', Sys.getpid(), '-', ii))); Sys.sleep(2); ii }",
    "environment(f) <- list2env(list(work = work), parent = baseenv())",
    "invisible(rnaparallel:::combat_parallel_lapply(as.list(1:4), f, 2L, 'mclapply', cells = Inf,",
    "                                               min_cells = 0, preschedule = TRUE))"), script)
  libs <- paste0("R_LIBS=", shQuote(paste(.libPaths(), collapse = .Platform$path.sep)))
  system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)),
          stdout = FALSE, stderr = FALSE, env = c("NOT_CRAN=true", libs), wait = FALSE)

  chunks <- function() list.files(work, "^chunk-")
  deadline <- Sys.time() + 60
  while (length(chunks()) < 2L && Sys.time() < deadline) Sys.sleep(0.05)
  expect_length(chunks(), 2L)
  master <- as.integer(readLines(file.path(work, "master")))
  master_tmp <- readLines(file.path(work, "tempdir"))
  withr::defer(unlink(master_tmp, recursive = TRUE))
  forks <- unique(as.integer(sub("^chunk-([0-9]+)-.*$", "\\1", chunks())))
  withr::defer(dp_kill(forks))
  expect_true(dir.exists(master_tmp))
  tools::pskill(master, tools::SIGKILL)

  deadline <- Sys.time() + 20
  while (any(dp_alive(forks)) && Sys.time() < deadline) Sys.sleep(0.1)
  expect_false(any(dp_alive(forks)))
  expect_length(chunks(), 2L)
  expect_true(dir.exists(master_tmp))
})


# ---- parallel's RNG stream ----------------------------------------------------------------

test_that("a forking dispatch leaves parallel's L'Ecuyer stream where the caller had it", {
  skip_on_os("windows")
  withr::local_options(combat.progress = FALSE)
  withr::local_seed(1, .rng_kind = "L'Ecuyer-CMRG")
  rngenv <- utils::getFromNamespace("RNGenv", "parallel")
  had <- exists("LEcuyer.seed", envir = rngenv, inherits = FALSE)
  saved <- get0("LEcuyer.seed", envir = rngenv, inherits = FALSE)
  withr::defer(if (had) assign("LEcuyer.seed", saved, envir = rngenv)
               else if (exists("LEcuyer.seed", envir = rngenv, inherits = FALSE)) rm("LEcuyer.seed", envir = rngenv))

  idx <- rnaparallel:::combat_row_chunks(12L, chunks = 3L)
  backends <- "mclapply"
  if (requireNamespace("BiocParallel", quietly = TRUE)) backends <- c(backends, "BiocParallel")
  for (b in backends) {
    parallel::mc.reset.stream()
    l0 <- get("LEcuyer.seed", envir = rngenv, inherits = FALSE)
    out <- rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), 2L, b,
                                                cells = Inf, min_cells = 0)
    expect_identical(out, lapply(idx, sum), info = b)
    expect_identical(get0("LEcuyer.seed", envir = rngenv, inherits = FALSE), l0, info = b)

    rm("LEcuyer.seed", envir = rngenv)
    rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), 2L, b, cells = Inf, min_cells = 0)
    expect_false(exists("LEcuyer.seed", envir = rngenv, inherits = FALSE), info = b)
  }
})


# ---- the progress reporter ---------------------------------------------------------------

test_that("an error while the reporter starts still runs the reporter's stop", {
  skip_on_os("windows")
  withr::local_options(combat.progress = TRUE, combat.progress.dir = withr::local_tempdir())
  rd <- rnaparallel:::.rp_dispatch
  withr::defer(rd$reporter <- FALSE)
  local_mocked_bindings(rp_reporter_start = function(dir) {
    rd$reporter <- TRUE
    stop("reporter start failed")
  })
  idx <- rnaparallel:::combat_row_chunks(12L, chunks = 2L)
  for (b in list(function(idx, f, workers) lapply(idx, f), "mclapply")) {
    rd$reporter <- FALSE
    expect_error(rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), 2L, b,
                                                      cells = Inf, min_cells = 0),
                 "reporter start failed")
    expect_false(isTRUE(rd$reporter), info = if (is.character(b)) b else "custom")
  }
})

test_that("a future plan that resolves in one process starts no reporter", {
  skip_on_os("windows")
  skip_if_not_installed("future.apply")
  withr::local_options(combat.progress = TRUE, combat.progress.dir = withr::local_tempdir())
  local_mocked_bindings(rp_reporter_visible = function() TRUE, .package = "rnaparallel")
  old <- future::plan(future::sequential)
  withr::defer(future::plan(old))
  idx <- rnaparallel:::combat_row_chunks(12L, chunks = 2L)
  got <- dp_conds(rnaparallel:::combat_parallel_lapply(
    idx, function(i) isTRUE(rnaparallel:::.rp_dispatch$reporter), 2L, "future",
    cells = Inf, min_cells = 0))
  expect_identical(got$value, list(FALSE, FALSE))
  expect_true(any(grepl("resolves in one process", got$warnings)))
})


# ---- row splits that stand down --------------------------------------------------------

dp_counts <- function() {
  set.seed(5)
  matrix(rnbinom(40 * 12, mu = 60, size = 4), 40, 12,
         dimnames = list(paste0("g", 1:40), paste0("s", 1:12)))
}

test_that("a failed gene's glmFit refit raises each warning once and says it stood down", {
  orig <- edgeR::glmFit.default
  fake <- function(y, ...) {
    warning("dp glmFit warning")
    fit <- orig(y, ...)
    hit <- rownames(y) == "g7"
    if (any(hit)) fit$failed[hit] <- TRUE
    fit
  }
  dp_swap_edger("glmFit.default", fake)
  withr::defer(dp_swap_edger("glmFit.default", orig))

  y <- dp_counts()
  design <- stats::model.matrix(~ stats::rnorm(12))
  off <- log(colSums(y))
  ref <- dp_conds(edgeR::glmFit.default(y, design = design, dispersion = 0.1, offset = off,
                                        lib.size = NULL, weights = NULL, prior.count = 0.125,
                                        start = NULL))
  expect_gt(length(ref$warnings), 0L)
  for (b in names(dp_backends())) {
    rnaparallel:::rp_count_reset()
    got <- dp_conds(rnaparallel:::glmFit_rows_parallel(
      y, design = design, dispersion = 0.1, offset = off, workers = 2L, chunks = 4L,
      parallel_backend = dp_backends()[[b]]))
    expect_identical(got$value, ref$value, info = b)
    expect_identical(got$warnings, ref$warnings, info = b)
    expect_true("glmFit" %in% rnaparallel:::.rp_dispatch$fallback, info = b)
  }
})

test_that("a non-finite tagwise split recomputes unsplit, raises each warning once and says so", {
  orig <- utils::getS3method("estimateGLMTagwiseDisp", "default", envir = asNamespace("edgeR"))
  fake <- function(y, ...) {
    warning("dp tagwise warning")
    d <- orig(y, ...)
    d[rownames(y) == "g7"] <- NaN
    d
  }
  dp_swap_edger("estimateGLMTagwiseDisp.default", fake)
  withr::defer(dp_swap_edger("estimateGLMTagwiseDisp.default", orig))

  y <- dp_counts()
  design <- stats::model.matrix(~ stats::rnorm(12))
  off <- log(colSums(y))
  span <- (10 / nrow(y))^0.23
  alc <- edgeR::aveLogCPM(y, offset = off)
  ref <- dp_conds(edgeR::estimateGLMTagwiseDisp(y, design = design, offset = off, dispersion = 0.1,
                                                prior.df = 0, span = span, AveLogCPM = alc))
  expect_gt(length(ref$warnings), 0L)
  expect_true(is.nan(ref$value[7L]))
  for (b in names(dp_backends())) {
    rnaparallel:::rp_count_reset()
    got <- dp_conds(rnaparallel:::estimateGLMTagwiseDisp_rows_parallel(
      y, design = design, dispersion = 0.1, offset = off, prior.df = 0, span = span,
      AveLogCPM = alc, workers = 2L, chunks = 4L, parallel_backend = dp_backends()[[b]]))
    expect_identical(got$value, ref$value, info = b)
    expect_identical(got$warnings, ref$warnings, info = b)
    expect_true("estimateGLMTagwiseDisp" %in% rnaparallel:::.rp_dispatch$fallback, info = b)
  }
})
