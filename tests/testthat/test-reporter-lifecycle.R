## The forked progress reporter: where a dispatch's rows go, what may end it, and what rp_reporter_stop() may signal.

rp_lc_running <- function(pid) {
  st <- tryCatch(ps::ps_status(ps::ps_handle(as.integer(pid))), error = function(e) NA_character_)
  !is.na(st) && st != "zombie"
}

rp_lc_gone <- function(pid, secs = 3) {
  deadline <- Sys.time() + secs
  while (rp_lc_running(pid) && Sys.time() < deadline) Sys.sleep(0.05)
  !rp_lc_running(pid)
}

rp_lc_dispatch <- function(exec) {
  idx <- rnaparallel:::combat_row_chunks(12L, chunks = 3L)
  rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), workers = 2L,
                                       parallel_backend = exec, cells = Inf, min_cells = 0)
}

rp_lc_reporter <- function(d) {
  h <- NULL
  invisible(utils::capture.output(h <- with_mocked_bindings(
    rnaparallel:::rp_reporter_start(d), rp_reporter_visible = function() TRUE,
    .package = "rnaparallel")))
  h
}


# ---- where the rows go ---------------------------------------------------------------

test_that("an unnamed progress directory is new for each dispatch and gone once it returns", {
  skip_on_os("windows")
  local_mocked_bindings(rp_reporter_visible = function() TRUE, .package = "rnaparallel")
  for (progress in c(FALSE, TRUE)) {
    withr::local_options(combat.progress.dir = NULL, combat.progress = progress)
    seen <- character()
    rows <- integer()
    exec <- function(idx, f, w) {
      pd <- get("progress_dir", envir = environment(f))
      out <- lapply(idx, f)
      seen <<- c(seen, pd)
      rows <<- c(rows, length(unlist(lapply(
        list.files(pd, "^rnaparallel-.*\\.tsv$", full.names = TRUE), readLines))))
      out
    }
    rp_lc_dispatch(exec)
    rp_lc_dispatch(exec)
    info <- paste("combat.progress =", progress)
    expect_identical(rows, c(6L, 6L), info = info)
    expect_length(unique(seen), 2L)
    expect_identical(dirname(seen), rep(file.path(tempdir(), "rnaparallel-progress"), 2L),
                     info = info)
    expect_false(any(dir.exists(seen)), info = info)
  }
})

test_that("a companion call that stays serial leaves no unnamed progress directory behind", {
  skip_if_not_installed("limma")
  withr::local_options(combat.progress.dir = NULL, combat.progress = FALSE)
  set.seed(1)
  y <- matrix(rnorm(200), 20, 10)
  out <- removeBatchEffect_parallel(y, batch = rep(1:2, 5), workers = 1L)
  expect_identical(out, limma::removeBatchEffect(y, batch = rep(1:2, 5)))
  expect_false(dir.exists(rnaparallel:::rp_progress_own_dir()))
})

test_that("a straggler writing into a removed dispatch directory stays silent", {
  withr::local_options(warn = 1)
  gone <- file.path(withr::local_tempdir(), "gone")
  expect_no_warning(rnaparallel:::rp_progress_file_write(gone, "stage", 1L, "start"))
  expect_false(dir.exists(gone))
})

test_that("a caller-named progress directory keeps every row after the dispatch", {
  d <- withr::local_tempdir()
  withr::local_options(combat.progress.dir = d, combat.progress = FALSE)
  rp_lc_dispatch(function(idx, f, w) lapply(idx, f))
  files <- list.files(d, "^rnaparallel-.*\\.tsv$", full.names = TRUE)
  expect_length(files, 1L)
  expect_length(readLines(files[1L]), 6L)
})


# ---- who can see the reporter -----------------------------------------------------------

rp_lc_seen <- function(visible) {
  local_mocked_bindings(rp_reporter_visible = function() visible, .package = "rnaparallel")
  seen <- NULL
  exec <- function(idx, f, w) {
    out <- lapply(idx, f)
    pd <- get("progress_dir", envir = environment(f))
    seen <<- list(reporter = isTRUE(rnaparallel:::.rp_dispatch$reporter), dir = pd,
                  rows = length(unlist(lapply(
                    list.files(pd, "^rnaparallel-.*\\.tsv$", full.names = TRUE), readLines))))
    out
  }
  out <- rp_lc_dispatch(exec)
  c(seen, list(out = out))
}

test_that("a dispatch nobody can see starts no reporter, with the same result and progress rows", {
  skip_on_os("windows")
  withr::local_envvar(RNAPARALLEL_IN_WORKER = NA)
  withr::local_options(combat.progress = TRUE, combat.fork = TRUE, combat.progress.dir = NULL)
  hidden <- rp_lc_seen(FALSE)
  shown <- rp_lc_seen(TRUE)
  expect_false(hidden$reporter)
  expect_true(shown$reporter)
  expect_identical(hidden$out, shown$out)
  expect_identical(hidden$out, lapply(rnaparallel:::combat_row_chunks(12L, chunks = 3L), sum))
  expect_identical(c(hidden$rows, shown$rows), c(6L, 6L))
  expect_false(dir.exists(hidden$dir))
  expect_false(dir.exists(shown$dir))
})

test_that("an Rscript run whose output goes to a log gets progress rows and no reporter frames", {
  skip_on_cran()
  skip_on_os("windows")
  d <- withr::local_tempdir()
  dir.create(file.path(d, "progress"))
  script <- file.path(d, "run.R")
  log <- file.path(d, "run.log")
  writeLines(c(
    rp_load_line(),
    sprintf("options(combat.progress = TRUE, combat.progress.dir = %s)",
            deparse(file.path(d, "progress"))),
    "exec <- function(idx, f, w) { out <- lapply(idx, f); Sys.sleep(1.5); out }",
    "r <- rnaparallel:::combat_parallel_lapply(list(1:3, 4:6), function(i) sum(i), 2L, exec,",
    "                                          cells = Inf, min_cells = 0)",
    "cat('\\nRESULT', identical(r, list(6L, 15L)), '\\n')"), script)
  libs <- paste0("R_LIBS=", shQuote(paste(.libPaths(), collapse = .Platform$path.sep)))
  suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)),
                           stdout = log, stderr = log, env = c("NOT_CRAN=true", libs), timeout = 60))
  out <- readChar(log, file.size(log), useBytes = TRUE)
  expect_match(out, "RESULT TRUE", fixed = TRUE)
  expect_false(grepl("\\|[=-]{50}\\||waiting for", out), info = out)
  rows <- unlist(lapply(list.files(file.path(d, "progress"), "^rnaparallel-.*\\.tsv$",
                                   full.names = TRUE), readLines))
  expect_length(rows, 4L)
})


# ---- the reporter's view ---------------------------------------------------------------

test_that("the reporter's watch starts from a cache that marks every file already present", {
  skip_on_os("windows")
  withr::local_envvar(RNAPARALLEL_IN_WORKER = NA)
  withr::local_options(combat.progress = TRUE, combat.fork = TRUE)
  d <- withr::local_tempdir()
  stale <- file.path(d, "rnaparallel-1.tsv")
  writeLines("1\tearlier stage\t1\tstart", stale)
  out <- withr::local_tempfile()
  local_mocked_bindings(rp_progress_watch = function(dir, interval, stall_after,
                                                     cache = NULL, alive = NULL) {
    sizes <- if (is.environment(cache)) unlist(eapply(cache, function(e) e$size)) else numeric()
    saveRDS(sizes, paste0(out, ".part"))
    file.rename(paste0(out, ".part"), out)
    invisible(NULL)
  }, .package = "rnaparallel")
  h <- rp_lc_reporter(d)
  withr::defer(rnaparallel:::rp_reporter_stop(h))
  expect_false(is.null(h))
  deadline <- Sys.time() + 5
  while (!file.exists(out) && Sys.time() < deadline) Sys.sleep(0.05)
  expect_true(file.exists(out))
  expect_identical(readRDS(out), stats::setNames(file.size(stale), stale))
})

test_that("a watch whose view holds more done rows than start rows draws and returns", {
  d <- withr::local_tempdir()
  writeLines(c("1\ta\t1\tstart", "2\ta\t1\tdone", "2\tprev\t7\tdone"), file.path(d, "rnaparallel-1.tsv"))
  res <- NULL
  invisible(utils::capture.output(err <- tryCatch({
    res <- rnaparallel_progress(d, watch = TRUE, interval = 0.05, stall_after = 10)
    NULL
  }, error = function(e) conditionMessage(e))))
  expect_null(err)
  expect_identical(c(res$done, res$started), c(2L, 1L))
})

test_that("the reporter mode of the watch leaves a finished bar on its line, without a newline", {
  d <- withr::local_tempdir()
  writeLines(c("1\ta\t1\tstart", "2\ta\t1\tdone"), file.path(d, "rnaparallel-1.tsv"))
  watch_bytes <- function(alive) {
    f <- withr::local_tempfile()
    con <- file(f, "w")
    sink(con)
    tryCatch(rnaparallel:::rp_progress_watch(d, 0.05, 30, alive = alive),
             finally = { sink(); close(con) })
    readChar(f, file.size(f), useBytes = TRUE)
  }
  in_reporter <- watch_bytes(function() TRUE)
  expect_match(in_reporter, "100%", fixed = TRUE)
  expect_false(grepl("\n", in_reporter, fixed = TRUE))
  standalone <- watch_bytes(NULL)
  expect_match(standalone, "100%", fixed = TRUE)
  expect_true(endsWith(standalone, "\n"))
})


# ---- what may end the reporter -----------------------------------------------------------

test_that("an error inside the reporter's watch does not end the reporter", {
  skip_on_os("windows")
  skip_if_not_installed("ps")
  withr::local_envvar(RNAPARALLEL_IN_WORKER = NA)
  withr::local_options(combat.progress = TRUE, combat.fork = TRUE)
  local_mocked_bindings(rp_progress_watch = function(...) stop("watch failed"),
                        .package = "rnaparallel")
  h <- rp_lc_reporter(withr::local_tempdir())
  withr::defer(rnaparallel:::rp_reporter_stop(h))
  expect_false(is.null(h))
  Sys.sleep(1)
  expect_true(rp_lc_running(h$pid))
})

test_that("a done-only row from an earlier dispatch does not end the reporter", {
  skip_on_os("windows")
  skip_if_not_installed("ps")
  withr::local_envvar(RNAPARALLEL_IN_WORKER = NA)
  withr::local_options(combat.progress = TRUE, combat.fork = TRUE)
  d <- withr::local_tempdir()
  writeLines(c("1\tprev\t1\tstart", "1\tprev\t2\tstart"), file.path(d, "rnaparallel-1.tsv"))
  h <- rp_lc_reporter(d)
  withr::defer(rnaparallel:::rp_reporter_stop(h))
  expect_false(is.null(h))
  cat("2\tprev\t1\tdone\n2\tprev\t2\tdone\n", file = file.path(d, "rnaparallel-1.tsv"), append = TRUE)
  cat("2\tcur\t1\tstart\n", file = file.path(d, "rnaparallel-2.tsv"))
  Sys.sleep(1.5)
  expect_true(rp_lc_running(h$pid))
})

test_that("an interrupt does not end the reporter, in its watch or after it", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("ps")
  withr::local_envvar(RNAPARALLEL_IN_WORKER = NA)
  withr::local_options(combat.progress = TRUE, combat.fork = TRUE)
  local_mocked_bindings(rp_progress_watch = function(...) Sys.sleep(1),
                        .package = "rnaparallel")
  h <- rp_lc_reporter(withr::local_tempdir())
  withr::defer(rnaparallel:::rp_reporter_stop(h))
  expect_false(is.null(h))
  Sys.sleep(0.6)
  tools::pskill(h$pid, tools::SIGINT)
  Sys.sleep(1.2)
  expect_true(rp_lc_running(h$pid))
  tools::pskill(h$pid, tools::SIGINT)
  Sys.sleep(0.5)
  expect_true(rp_lc_running(h$pid))
})


# ---- what rp_reporter_stop may signal -----------------------------------------------------

test_that("rp_reporter_stop leaves alone a pid that no longer names the reporter it started", {
  skip_on_os("windows")
  skip_if_not_installed("ps")
  withr::local_envvar(RNAPARALLEL_IN_WORKER = NA)
  withr::local_options(combat.progress = TRUE, combat.fork = TRUE)
  local_mocked_bindings(rp_progress_watch = function(...) invisible(NULL), .package = "rnaparallel")
  h <- rp_lc_reporter(withr::local_tempdir())
  withr::defer(rnaparallel:::rp_reporter_stop(h))
  expect_false(is.null(h))
  victim <- parallel::mcparallel(Sys.sleep(30))
  withr::defer({
    tools::pskill(victim$pid, tools::SIGKILL)
    suppressWarnings(parallel::mccollect(victim, wait = TRUE))
  })
  fake <- h
  fake$pid <- victim$pid
  rnaparallel:::rp_reporter_stop(fake)
  Sys.sleep(0.3)
  expect_true(rp_lc_running(victim$pid))
  rnaparallel:::rp_reporter_stop(h)
  expect_true(rp_lc_gone(h$pid))
})

test_that("rp_reporter_stop(NULL) stops the reporter that rp_reporter_start recorded", {
  skip_on_os("windows")
  skip_if_not_installed("ps")
  withr::local_envvar(RNAPARALLEL_IN_WORKER = NA)
  withr::local_options(combat.progress = TRUE, combat.fork = TRUE)
  local_mocked_bindings(rp_progress_watch = function(...) invisible(NULL), .package = "rnaparallel")
  h <- rp_lc_reporter(withr::local_tempdir())
  withr::defer(if (rp_lc_running(h$pid)) tools::pskill(h$pid, tools::SIGKILL))
  expect_false(is.null(h))
  expect_true(rp_lc_running(h$pid))
  rnaparallel:::rp_reporter_stop(NULL)
  expect_true(rp_lc_gone(h$pid))
  expect_false(isTRUE(rnaparallel:::.rp_dispatch$reporter))
})
