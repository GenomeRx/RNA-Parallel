## The memory guard: rp_mem_available(), rp_mem_rss(), rp_mem_cap(), rp_mem_peak().
##
## /proc is the Linux source; the `ps` package is the cross-platform fallback (same pattern
## as rp_getppid() below), so these are NA only when NEITHER is available, not simply "off
## Linux" -- this suite runs on Windows too, and with `ps` installed (it is a Suggests
## dependency) these now return real numbers there, which is the whole point of having the
## fallback. Real arithmetic elsewhere in this file is exercised by mocking
## rp_mem_available()/rp_mem_rss() rather than depending on actual system memory pressure,
## which is not reproducible in CI.

test_that("rp_mem_available and rp_mem_rss return NA only with neither /proc nor ps", {
  skip_if(file.exists("/proc/meminfo"), "this machine has /proc; not exercising the NA path")
  skip_if(requireNamespace("ps", quietly = TRUE),
         "ps is installed here, so real values are expected -- see the next test")
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
         "ps is installed here, so a real value is expected -- see the next test")
  expect_true(is.na(rnaparallel:::rp_mem_peak()))
})

test_that("rp_mem_peak returns a real number via ps off Linux", {
  skip_if(file.exists("/proc/self/status"), "this machine has /proc; ps fallback not exercised here")
  skip_if_not_installed("ps")
  peak <- rnaparallel:::rp_mem_peak()
  expect_true(is.numeric(peak) && !is.na(peak) && peak > 0)
})

test_that("rp_mem_cap is a no-op when it cannot read memory (NA means proceed)", {
  # Works on every platform: on Linux the mock stands in for a bad reading; off Linux
  # rp_mem_available()/rp_mem_rss() already return NA on their own.
  testthat::local_mocked_bindings(
    rp_mem_available = function() NA_real_,
    rp_mem_rss = function() NA_real_,
    .package = "rnaparallel"
  )
  expect_identical(rnaparallel:::rp_mem_cap(16L), 16L)
})

test_that("rp_mem_cap is a no-op for workers = 1 regardless of memory", {
  testthat::local_mocked_bindings(
    rp_mem_available = function() 1 * 2^30,   # 1 GB available
    rp_mem_rss = function() 100 * 2^30,       # a 100 GB parent
    .package = "rnaparallel"
  )
  # one worker cannot fork-multiply anything; nothing to cap
  expect_identical(rnaparallel:::rp_mem_cap(1L), 1L)
})

test_that("rp_mem_cap passes workers through when the fit clears headroom", {
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
  # The PR that added this guard measured 4 workers off a 23 GB parent needing 111 GB, a
  # real per-worker divergence of about 0.96: at the old 0.25 default this exact case
  # computed only 23 GB needed and would have proceeded unwarned into the same kill.
  testthat::local_mocked_bindings(
    rp_mem_available = function() 30 * 2^30,   # 30 GB available (headroom 24 GB)
    rp_mem_rss = function() 23 * 2^30,         # the PR's own measured 23 GB parent
    .package = "rnaparallel"
  )
  withr::local_options(combat.mem.divergence = NULL)   # no explicit option: use the default
  # at divergence = 1: 4 * 23 GB = 92 GB needed, far past 24 GB headroom -> must degrade
  expect_warning(fit <- rnaparallel:::rp_mem_cap(4L), "need")
  expect_true(fit < 4L)
  # at the OLD 0.25 default this same case computed 4 * 23 * 0.25 = 23 GB needed, which is
  # under 24 GB headroom and would NOT have degraded -- confirming the default actually
  # changed behavior on the case it exists to catch, not just the option's stated value
  withr::local_options(combat.mem.divergence = 0.25)
  expect_identical(rnaparallel:::rp_mem_cap(4L), 4L)
})

test_that("rp_mem_cap degrades workers instead of leaving a fit that would be killed", {
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
  testthat::local_mocked_bindings(
    rp_mem_available = function() 10 * 2^30,
    rp_mem_rss = function() 5 * 2^30,
    .package = "rnaparallel"
  )
  withr::local_options(combat.mem.divergence = 0.5)
  expect_warning(rnaparallel:::rp_mem_cap(16L),
                regexp = "16 workers.*GB.*GB.*GB", perl = TRUE)
})

## rp_getppid(): the orphan-fork exit check reads this. It does not exist as a base R
## function on every build (confirmed FALSE via exists() on the R 4.6.1 UCRT Windows build
## this package is tested on), which crashed every single dispatch under a real
## future::multisession run before this fallback existed -- caught fixing the bar, not a
## hypothetical.

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
