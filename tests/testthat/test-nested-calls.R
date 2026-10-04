## Each nested companion call is capped once, hides no option change, and reports the worker count it used.

skip_r1 <- function() {
  testthat::skip_if_not_installed("limma")
  testthat::skip_if_not_installed("edgeR")
  testthat::skip_if_not_installed("statmod")
}

collect <- function(expr) {
  warns <- character()
  msgs <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  list(value = value, warnings = warns, messages = msgs)
}

test_that("a nested companion call leaves the caller's combat.mem.guard visible to its own code", {
  skip_r1()
  set.seed(7)
  y <- matrix(stats::rnorm(400 * 24, 8, 2), 400, 24)
  b <- factor(rep(1:4, length.out = 24))
  set.seed(12)
  G <- 240L
  S <- 12L
  grp <- factor(rep(c("A", "B"), each = S / 2))
  des <- stats::model.matrix(~ grp)
  blk <- factor(rep(seq_len(S / 2), each = 2))
  M <- matrix(stats::rnorm(G * S, 8, 2), G, S)
  ref_rbe <- limma::removeBatchEffect(y, batch = b)
  ref_fit <- limma::lmFit(M, des, block = blk, correlation = NULL)
  local_mocked_bindings(rp_mem_available = function() 1e15, rp_mem_rss = function() 1)
  for (caller in list(NULL, TRUE)) {
    withr::local_options(combat.mem.guard = caller, combat.progress = FALSE)
    seen <- list()
    spy <- function(idx, f, workers) {
      seen[length(seen) + 1L] <<- list(getOption("combat.mem.guard"))
      lapply(idx, f)
    }
    tag <- paste("caller value", deparse(caller))
    expect_identical(removeBatchEffect_parallel(y, batch = b, workers = 2L, chunks = 3L,
                                                parallel_backend = spy), ref_rbe, info = tag)
    expect_true(length(seen) > 0L, info = tag)
    expect_true(all(vapply(seen, identical, logical(1), caller)), info = tag)
    seen <- list()
    expect_identical(lmFit_parallel(M, des, block = blk, correlation = NULL, workers = 2L,
                                    chunks = 3L, parallel_backend = spy), ref_fit, info = tag)
    expect_identical(length(seen), 2L, info = tag)
    expect_true(all(vapply(seen, identical, logical(1), caller)), info = tag)
    expect_identical(getOption("combat.mem.guard"), caller, info = tag)
  }
})

test_that("rp_mem_cap passes workers through inside rp_uncapped, nested or not, and caps after", {
  local_mocked_bindings(rp_mem_available = function() 250 / 0.8, rp_mem_rss = function() 100)
  withr::local_options(combat.mem.guard = TRUE, combat.mem.divergence = 1)
  inner <- rp_uncapped({ rp_uncapped(NULL); rp_mem_cap(3L) })
  expect_identical(inner, 3L)
  expect_warning(after <- rp_mem_cap(3L), "workers need")
  expect_identical(after, 2L)
  withr::local_options(combat.mem.guard = "yes")
  expect_error(rp_uncapped(rp_mem_cap(3L)), "must be TRUE or FALSE")
})

test_that("calcNormFactors_parallel on a DGEList caps workers once and reports the count it used", {
  skip_r1()
  set.seed(5)
  counts <- matrix(stats::rnbinom(400 * 24, mu = 50, size = 5), 400, 24)
  dge <- edgeR::DGEList(counts)
  edger_norm <- if (exists("normLibSizes.default", envir = asNamespace("edgeR"),
                           inherits = FALSE)) edgeR::normLibSizes else edgeR::calcNormFactors
  ref <- edger_norm(dge)
  rss_calls <- 0L
  local_mocked_bindings(
    rp_mem_available = function() 250 / 0.8,
    rp_mem_rss = function() {
      rss_calls <<- rss_calls + 1L
      100 + 30 * (rss_calls - 1L)
    })
  seen <- integer()
  spy <- function(idx, f, workers) { seen <<- c(seen, workers); lapply(idx, f) }
  withr::local_options(combat.timing = TRUE, combat.timing.min = 0, combat.mem.guard = TRUE,
                       combat.mem.divergence = 1, combat.min.norm.cells = 0,
                       combat.progress = FALSE)
  got <- collect(calcNormFactors_parallel(dge, workers = 3L, chunks = 3L, parallel_backend = spy))
  expect_identical(got$value, ref)
  expect_length(grep("workers need", got$warnings), 1L)
  expect_identical(rss_calls, 1L)
  expect_true(length(seen) > 0L && all(seen == 2L))
  expect_true(any(grepl("custom x2", got$messages, fixed = TRUE)))
  withr::local_options(combat.mem.guard = "yes")
  expect_error(calcNormFactors_parallel(dge, workers = 2L), "must be TRUE or FALSE")
})

test_that("the timing line after a capped nested call names the worker count it used", {
  withr::local_options(combat.timing = TRUE, combat.timing.min = 0, combat.progress = FALSE,
                       combat.quiet = FALSE)
  expect_identical(rp_or0(.rp_dispatch$depth), 0L)
  be <- function(idx, f, workers) NULL
  got <- collect({
    outer <- rp_step_begin("r1 outer", "unit", matrix(1, 2, 2), be, 3L)
    inner <- rp_step_begin(NULL, "unit", matrix(1, 2, 2), be, 2L)
    rp_count(TRUE)
    rp_step_end(inner)
    rp_step_end(outer)
  })
  line <- grep("r1 outer", got$messages, value = TRUE)
  expect_length(line, 1L)
  expect_match(line, "custom x2", fixed = TRUE)
  got <- collect({
    outer <- rp_step_begin("r1 alone", "unit", matrix(1, 2, 2), be, 4L)
    rp_count(TRUE)
    rp_step_end(outer)
  })
  expect_match(grep("r1 alone", got$messages, value = TRUE), "custom x4", fixed = TRUE)
  expect_identical(rp_or0(.rp_dispatch$depth), 0L)
})

test_that("an unwritable combat.progress.dir is refused without a base R warning", {
  skip_on_os("windows")
  base <- withr::local_tempdir()
  target <- file.path(base, "locked")
  dir.create(target)
  Sys.chmod(target, mode = "0500")
  withr::defer(Sys.chmod(target, mode = "0700"))
  skip_if(file.access(target, 2L) == 0L, "this user can write to a mode 0500 directory")
  withr::local_options(combat.progress.dir = target)
  got <- collect(tryCatch(rp_step_begin(NULL, "unit", matrix(1, 2, 2), "serial", 1L),
                          error = function(e) conditionMessage(e)))
  expect_match(got$value, "exists but is not writable", fixed = TRUE)
  expect_identical(got$warnings, character())
})

test_that("combat.progress = FALSE starts no forked progress reporter", {
  skip_on_os("windows")
  withr::local_envvar(RNAPARALLEL_IN_WORKER = NA)
  withr::local_options(combat.progress = FALSE, combat.fork = TRUE)
  local_mocked_bindings(rp_reporter_visible = function() TRUE, .package = "rnaparallel")
  dir <- withr::local_tempdir()
  h <- rnaparallel:::rp_reporter_start(dir)
  withr::defer(if (!is.null(h)) rnaparallel:::rp_reporter_stop(h))
  expect_null(h)
})
