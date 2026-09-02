## The memory guard: rp_mem_available(), rp_mem_rss(), rp_mem_cap(), rp_mem_peak().
##
## /proc only exists on Linux, so most of this asserts the NA/pass-through path on any other
## platform (this suite runs on Windows too) and exercises the real arithmetic by mocking
## rp_mem_available()/rp_mem_rss() rather than depending on actual system memory pressure,
## which is not reproducible in CI.

test_that("rp_mem_available and rp_mem_rss return NA off Linux", {
  skip_if(file.exists("/proc/meminfo"), "this machine has /proc; NA path not exercised here")
  expect_true(is.na(rnaparallel:::rp_mem_available()))
  expect_true(is.na(rnaparallel:::rp_mem_rss()))
})

test_that("rp_mem_peak returns NA off Linux", {
  skip_if(file.exists("/proc/self/status"), "this machine has /proc; NA path not exercised here")
  expect_true(is.na(rnaparallel:::rp_mem_peak()))
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
