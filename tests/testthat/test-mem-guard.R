# The memory readers use /proc on Linux and ps elsewhere, so they are NA only with neither, and the cap arithmetic below mocks them.

test_that("rp_mem_available and rp_mem_rss return NA only with neither /proc nor ps", {
  skip_if(file.exists("/proc/meminfo"), "this machine has /proc; not exercising the NA path")
  skip_if(requireNamespace("ps", quietly = TRUE),
         "ps is installed here, so real values are expected. See the next test")
  expect_true(is.na(rnaparallel:::rp_mem_available()))
  expect_true(is.na(rnaparallel:::rp_mem_rss()))
})

test_that("rp_mem_available and rp_mem_rss return real numbers via ps off Linux", {
  skip_if(file.exists("/proc/meminfo"), "this machine has /proc; ps fallback not exercised here")
  skip_if_not_installed("ps")
  avail <- rnaparallel:::rp_mem_available()
  rss <- rnaparallel:::rp_mem_rss()
  expect_true(is.numeric(avail) && !is.na(avail) && avail > 0)
  expect_true(is.numeric(rss) && !is.na(rss) && rss > 0)
})

test_that("rp_mem_peak returns NA only with neither /proc nor ps", {
  skip_if(file.exists("/proc/self/status"), "this machine has /proc; not exercising the NA path")
  skip_if(requireNamespace("ps", quietly = TRUE),
         "ps is installed here, so a real value is expected. See the next test")
  expect_true(is.na(rnaparallel:::rp_mem_peak()))
})

test_that("rp_mem_peak returns a real number via ps off Linux", {
  skip_if(file.exists("/proc/self/status"), "this machine has /proc; ps fallback not exercised here")
  skip_if_not_installed("ps")
  peak <- rnaparallel:::rp_mem_peak()
  expect_true(is.numeric(peak) && !is.na(peak) && peak > 0)
})

test_that("rp_mem_cap is a no-op when it cannot read memory (NA means proceed)", {
  withr::local_options(combat.mem.guard = TRUE)
  testthat::local_mocked_bindings(
    rp_mem_available = function() NA_real_,
    rp_mem_rss = function() NA_real_,
    .package = "rnaparallel"
  )
  expect_identical(rnaparallel:::rp_mem_cap(16L), 16L)
})

test_that("rp_mem_cap is a no-op for workers = 1 regardless of memory", {
  withr::local_options(combat.mem.guard = TRUE)
  testthat::local_mocked_bindings(
    rp_mem_available = function() 1 * 2^30,   # 1 GB available
    rp_mem_rss = function() 100 * 2^30,       # a 100 GB parent
    .package = "rnaparallel"
  )
  # one worker cannot fork-multiply anything; nothing to cap
  expect_identical(rnaparallel:::rp_mem_cap(1L), 1L)
})

test_that("rp_mem_cap passes workers through when the fit clears headroom", {
  withr::local_options(combat.mem.guard = TRUE)
  testthat::local_mocked_bindings(
    rp_mem_available = function() 100 * 2^30,  # 100 GB available
    rp_mem_rss = function() 1 * 2^30,          # 1 GB parent
    .package = "rnaparallel"
  )
  withr::local_options(combat.mem.divergence = 0.25)
  # 8 workers * 1 GB * 0.25 = 2 GB needed, well under 80 GB headroom (0.8 * 100 GB)
  expect_identical(rnaparallel:::rp_mem_cap(8L), 8L)
})

test_that("default combat.mem.divergence is 1, not the guard's earlier 0.25", {
  withr::local_options(combat.mem.guard = TRUE)
# A measured run of 4 workers off a 23 GB parent needed 111 GB, a per-worker divergence of about 0.96.
  testthat::local_mocked_bindings(
    rp_mem_available = function() 30 * 2^30,   # 30 GB available (headroom 24 GB)
    rp_mem_rss = function() 23 * 2^30,         # the PR's own measured 23 GB parent
    .package = "rnaparallel"
  )
  withr::local_options(combat.mem.divergence = NULL)   # no explicit option: use the default
  # at divergence = 1: 4 * 23 GB = 92 GB needed, far past 24 GB headroom -> must degrade
  expect_warning(fit <- rnaparallel:::rp_mem_cap(4L), "need")
  expect_true(fit < 4L)
# At divergence 0.25 the same case needs 4 * 23 * 0.25 = 23 GB, under the 24 GB headroom, so it must not degrade.
  withr::local_options(combat.mem.divergence = 0.25)
  expect_identical(rnaparallel:::rp_mem_cap(4L), 4L)
})

test_that("rp_mem_cap degrades workers instead of leaving a fit that would be killed", {
  withr::local_options(combat.mem.guard = TRUE)
  testthat::local_mocked_bindings(
    rp_mem_available = function() 10 * 2^30,   # 10 GB available
    rp_mem_rss = function() 5 * 2^30,          # 5 GB parent
    .package = "rnaparallel"
  )
  withr::local_options(combat.mem.divergence = 0.5)
  # 16 workers * 5 GB * 0.5 = 40 GB needed against 8 GB headroom (0.8 * 10 GB): must degrade
  expect_warning(fit <- rnaparallel:::rp_mem_cap(16L), "need")
  expect_true(fit < 16L)
  expect_true(fit >= 1L)
  # the fit itself must actually clear headroom, not just be smaller
  expect_true(fit * 5 * 2^30 * 0.5 <= 10 * 2^30 * 0.8)
})

test_that("rp_mem_cap never degrades below 1 worker", {
  withr::local_options(combat.mem.guard = TRUE)
  testthat::local_mocked_bindings(
    rp_mem_available = function() 1 * 2^20,    # 1 MB available: nothing fits
    rp_mem_rss = function() 50 * 2^30,         # 50 GB parent
    .package = "rnaparallel"
  )
  withr::local_options(combat.mem.divergence = 1)
  expect_warning(fit <- rnaparallel:::rp_mem_cap(16L), "need")
  expect_identical(fit, 1L)
})

test_that("combat.mem.guard = FALSE disables the cap even on a fit that would be killed", {
  testthat::local_mocked_bindings(
    rp_mem_available = function() 1 * 2^30,
    rp_mem_rss = function() 50 * 2^30,
    .package = "rnaparallel"
  )
  withr::local_options(combat.mem.guard = FALSE, combat.mem.divergence = 1)
  expect_identical(rnaparallel:::rp_mem_cap(16L), 16L)
})

test_that("combat.mem.divergence = 0 disables the cap", {
  testthat::local_mocked_bindings(
    rp_mem_available = function() 1 * 2^30,
    rp_mem_rss = function() 50 * 2^30,
    .package = "rnaparallel"
  )
  withr::local_options(combat.mem.guard = TRUE, combat.mem.divergence = 0)
  expect_identical(rnaparallel:::rp_mem_cap(16L), 16L)
})

test_that("a garbage combat.mem.divergence is refused, not silently ignored", {
  withr::local_options(combat.mem.guard = TRUE)
  testthat::local_mocked_bindings(
    rp_mem_available = function() 100 * 2^30,
    rp_mem_rss = function() 1 * 2^30,
    .package = "rnaparallel"
  )
  withr::local_options(combat.mem.divergence = "lots")
  expect_error(rnaparallel:::rp_mem_cap(8L), "must be a single non-negative number")
})

test_that("a garbage combat.mem.guard is refused, not silently ignored", {
  withr::local_options(combat.mem.guard = "yes")
  expect_error(rnaparallel:::rp_mem_cap(8L), "must be TRUE or FALSE")
})

test_that("rp_mem_cap is reached from rp_prologue, not just directly", {
  withr::local_options(combat.mem.guard = TRUE)
  testthat::local_mocked_bindings(
    rp_mem_available = function() 10 * 2^30,
    rp_mem_rss = function() 5 * 2^30,
    .package = "rnaparallel"
  )
  withr::local_options(combat.mem.divergence = 0.5)
  expect_warning(w <- rnaparallel:::rp_prologue(16L), "need")
  expect_true(w < 16L)
})

test_that("the mem-guard warning names all three numbers, not just the outcome", {
  withr::local_options(combat.mem.guard = TRUE)
  testthat::local_mocked_bindings(
    rp_mem_available = function() 10 * 2^30,
    rp_mem_rss = function() 5 * 2^30,
    .package = "rnaparallel"
  )
  withr::local_options(combat.mem.divergence = 0.5)
  expect_warning(rnaparallel:::rp_mem_cap(16L),
                regexp = "16 workers.*GB.*GB.*GB", perl = TRUE)
})

test_that("worker survives dispatch on a healthy run even when its own ppid is 1", {
  skip_if_not_installed("sva")
# Every Unix PSOCK worker has ppid 1 from birth, so this runs each real backend rather than a mock.
  d <- make_counts(22, G = 60, n_per_batch = c(4, 4))
  ref <- quietly(sva::ComBat_seq(d$counts, d$batch, group = NULL))
  needs <- c(mclapply = "parallel", future = "future.apply",
             BiocParallel = "BiocParallel", foreach = "doParallel")
  dropped <- character()
  for (be in setdiff(combat_backends(), "serial")) {
    pkg <- needs[[be]]
    if (!requireNamespace(pkg, quietly = TRUE)) {
      dropped <- c(dropped, be)
      next
    }
    got <- local({
      if (be == "future") local_socket_plan(2L)
      quietly(ComBat_seq_parallel(d$counts, d$batch, group = NULL,
                                  workers = 2L, parallel_backend = be))
    })
    expect_identical(got, ref, info = be)
  }
  if (length(dropped)) skip(paste("not installed, so not exercised:", paste(dropped, collapse = ", ")))
})

# The orphan-fork exit check runs only in parallel's own forks (isChild()) and quits when rp_getppid() differs from master_pid, so an NA ppid skips it.

test_that("rp_getppid never errors, even when Sys.getppid does not exist on this build", {
  expect_no_error(v <- rnaparallel:::rp_getppid())
  expect_true(is.na(v) || (is.numeric(v) && v > 0))
})

test_that("rp_getppid returns a real value when Sys.getppid is present", {
  skip_if_not(exists("Sys.getppid", where = baseenv(), mode = "function"),
             "this R build has no Sys.getppid; NA path covered by the test above")
  v <- rnaparallel:::rp_getppid()
  expect_true(is.numeric(v) && v > 0)
})

test_that("rp_getppid returns NA, not an error, when both Sys.getppid and ps are unavailable", {
  skip_if(exists("Sys.getppid", where = baseenv(), mode = "function"),
         "cannot hide a real base function that already exists on this build")
  testthat::local_mocked_bindings(
    requireNamespace = function(...) FALSE,
    .package = "base"
  )
  expect_true(is.na(rnaparallel:::rp_getppid()))
})

test_that("a fork whose ppid reads NA skips the orphan check and finishes its chunk", {
  skip_on_os("windows")
  withr::local_options(combat.progress = FALSE)
  d <- withr::local_tempdir()
  testthat::local_mocked_bindings(
    rp_getppid = function() { file.create(file.path(d, Sys.getpid())); NA_integer_ },
    .package = "rnaparallel"
  )
  exec <- function(idx, f, w) {
    jobs <- lapply(idx, function(i) parallel::mcparallel(f(i)))
    unname(parallel::mccollect(jobs))
  }
  idx <- rnaparallel:::combat_row_chunks(12L, chunks = 2L)
  out <- rnaparallel:::combat_parallel_lapply(idx, function(i) sum(i), workers = 2L,
                                              parallel_backend = exec, cells = Inf, min_cells = 0)
  expect_identical(out, lapply(idx, sum))
  expect_length(list.files(d), length(idx))
})
