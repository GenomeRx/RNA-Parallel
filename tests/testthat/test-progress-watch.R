## rp_progress_watch()'s activity tracking and return contract, plus rnaparallel_progress()'s
## argument validation on the watch = TRUE path.

test_that("rnaparallel_progress refuses a garbage interval instead of busy-polling or erroring later", {
  d <- withr::local_tempdir()
  expect_error(rnaparallel_progress(d, watch = TRUE, interval = 0), "positive number of seconds")
  expect_error(rnaparallel_progress(d, watch = TRUE, interval = -1), "positive number of seconds")
  expect_error(rnaparallel_progress(d, watch = TRUE, interval = "fast"), "positive number of seconds")
  expect_error(rnaparallel_progress(d, watch = TRUE, interval = NA_real_), "positive number of seconds")
})

test_that("rnaparallel_progress refuses a garbage stall_after", {
  d <- withr::local_tempdir()
  expect_error(rnaparallel_progress(d, watch = TRUE, stall_after = 0), "positive number of seconds")
  expect_error(rnaparallel_progress(d, watch = TRUE, stall_after = -5), "positive number of seconds")
})

test_that("watch = FALSE (the default) is unaffected by the new validation", {
  d <- withr::local_tempdir()
  # no files yet: one-shot path prints "no files" and returns quietly, same as before
  expect_message(v <- rnaparallel_progress(d), "no rnaparallel")
  expect_identical(v$done, 0L)
})

test_that("watch mode treats a done-only change as activity, not just a started-only change", {
  skip_on_cran()
  skip_on_os("windows")
  d <- withr::local_tempdir()
  path <- file.path(d, "rnaparallel-1.tsv")
  writeLines(sprintf("1\ta\t%d\tstart", 1:8), path)
# Started rows never change while done rows land every 0.6 s over 4.8 s against a 2.5 s stall limit, so a watch that tracks only started rows stops before the last of them.
  writer <- parallel::mcparallel({
    for (k in 1:8) {
      Sys.sleep(0.6)
      cat(sprintf("2\ta\t%d\tdone\n", k), file = path, append = TRUE)
    }
  })
  withr::defer(parallel::mccollect(writer, wait = TRUE))
  invisible(utils::capture.output(
    result <- rnaparallel_progress(d, watch = TRUE, interval = 0.05, stall_after = 2.5)))
  expect_identical(result$started, 8L)
  expect_identical(result$done, 8L)
})

test_that("watch mode returns the last real summary on a stall, not NULL", {
  # Before this fix, a stall returned invisible(NULL), throwing away data the summarize
  # step had already computed. One start row, no matching done: guaranteed to stall within
  # a short stall_after, and the return should still carry started = 1, stalled = 1.
  d <- withr::local_tempdir()
  writeLines("1\ta\t1\tstart", file.path(d, "rnaparallel-1.tsv"))
  result <- rnaparallel_progress(d, watch = TRUE, interval = 0.05, stall_after = 0.2)
  expect_false(is.null(result))
  expect_identical(result$started, 1L)
  expect_identical(result$stalled, 1L)
  expect_identical(result$done, 0L)
})

test_that("watch mode on an empty directory (never written to) still returns a list, not NULL", {
  d <- withr::local_tempdir()   # no worker files at all
  result <- rnaparallel_progress(d, watch = TRUE, interval = 0.05, stall_after = 0.2)
  expect_false(is.null(result))
  expect_identical(result$done, 0L)
  expect_identical(result$started, 0L)
})
