## rnaparallel_set_mem_limit(): rp_mem_total() and the .Renviron writer.
##
## All file-touching tests pass path= explicitly to a tempfile. path.expand("~") on Windows
## resolves via USERPROFILE, not HOME, so a test that tried to redirect this by setting HOME
## still wrote to the real ~/.Renviron -- caught building this, which is why path= exists and
## every test below uses it instead of relying on the default.

test_that("rp_mem_total returns NA when neither /proc nor a Windows reading is available", {
  skip_if(file.exists("/proc/meminfo"), "this machine has /proc; NA path not exercised here")
  skip_if(!identical(.Platform$OS.type, "windows"),
          "off Windows with no /proc there is no second source to fail")
  # cannot force the Windows PowerShell call itself to fail portably; this only asserts the
  # function returns a number type consistent with the contract when it CAN read something
  v <- rnaparallel:::rp_mem_total()
  expect_true(is.na(v) || (is.numeric(v) && v > 0))
})

test_that("a garbage fraction is refused, not silently ignored", {
  expect_error(rnaparallel_set_mem_limit(fraction = 0), "must be a single number in")
  expect_error(rnaparallel_set_mem_limit(fraction = 1.5), "must be a single number in")
  expect_error(rnaparallel_set_mem_limit(fraction = "half"), "must be a single number in")
  expect_error(rnaparallel_set_mem_limit(fraction = NA_real_), "must be a single number in")
})

test_that("a garbage path is refused, not silently ignored", {
  expect_error(rnaparallel_set_mem_limit(path = ""), "must be a single non-empty string")
  expect_error(rnaparallel_set_mem_limit(path = NA_character_), "must be a single non-empty string")
  expect_error(rnaparallel_set_mem_limit(path = character(0)), "must be a single non-empty string")
})

test_that("dry_run writes nothing", {
  testthat::local_mocked_bindings(rp_mem_total = function() 64 * 2^30, .package = "rnaparallel")
  p <- withr::local_tempfile()
  expect_false(file.exists(p))
  suppressMessages(rnaparallel_set_mem_limit(dry_run = TRUE, path = p))
  expect_false(file.exists(p))
})

test_that("NA total RAM writes nothing and returns NA invisibly", {
  testthat::local_mocked_bindings(rp_mem_total = function() NA_real_, .package = "rnaparallel")
  p <- withr::local_tempfile()
  expect_message(v <- rnaparallel_set_mem_limit(path = p), "could not read total RAM")
  expect_true(is.na(v))
  expect_false(file.exists(p))
})

test_that("halves total RAM and rounds to the nearest tier, on a fresh file", {
  testthat::local_mocked_bindings(rp_mem_total = function() 100 * 2^30, .package = "rnaparallel")
  p <- withr::local_tempfile()
  suppressMessages(rnaparallel_set_mem_limit(fraction = 0.5, path = p))
  # 100 GB * 0.5 = 50 GB, nearer to 64 than 32 in log space
  expect_identical(readLines(p), "R_MAX_VSIZE=64Gb")
})

test_that("preserves other lines and replaces only a pre-existing R_MAX_VSIZE line", {
  testthat::local_mocked_bindings(rp_mem_total = function() 32 * 2^30, .package = "rnaparallel")
  p <- withr::local_tempfile()
  writeLines(c("SOME_VAR=1", "R_MAX_VSIZE=999Gb", "# mentions R_MAX_VSIZE in a comment"), p)
  suppressMessages(rnaparallel_set_mem_limit(fraction = 0.5, path = p))
  out <- readLines(p)
  expect_identical(out[1], "SOME_VAR=1")
  expect_identical(out[2], "R_MAX_VSIZE=16Gb")           # replaced, not duplicated
  expect_identical(out[3], "# mentions R_MAX_VSIZE in a comment")  # comment untouched
  expect_length(out, 3L)
})

test_that("appends R_MAX_VSIZE when the file has no existing line for it", {
  testthat::local_mocked_bindings(rp_mem_total = function() 16 * 2^30, .package = "rnaparallel")
  p <- withr::local_tempfile()
  writeLines("OTHER=1", p)
  suppressMessages(rnaparallel_set_mem_limit(fraction = 0.5, path = p))
  out <- readLines(p)
  expect_identical(out[1], "OTHER=1")
  expect_identical(out[2], "R_MAX_VSIZE=8Gb")
})

test_that("creates the file when path does not exist yet", {
  testthat::local_mocked_bindings(rp_mem_total = function() 128 * 2^30, .package = "rnaparallel")
  p <- withr::local_tempfile()
  expect_false(file.exists(p))
  suppressMessages(rnaparallel_set_mem_limit(fraction = 0.5, path = p))
  expect_true(file.exists(p))
  expect_identical(readLines(p), "R_MAX_VSIZE=64Gb")
})

test_that("returns the chosen value invisibly in bytes", {
  testthat::local_mocked_bindings(rp_mem_total = function() 64 * 2^30, .package = "rnaparallel")
  p <- withr::local_tempfile()
  v <- suppressMessages(rnaparallel_set_mem_limit(fraction = 0.5, path = p))
  expect_identical(v, 32 * 2^30)
})

test_that("rp_mem_peak returns NA off Linux", {
  skip_if(file.exists("/proc/self/status"), "this machine has /proc; NA path not exercised here")
  expect_true(is.na(rnaparallel:::rp_mem_peak()))
})
