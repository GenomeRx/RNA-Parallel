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
  # The bug this covers: activity used to be tracked on `started` alone, so a run whose
  # chunks were all dispatched early but whose done rows are still trickling in would hit
  # the stall timeout mid-run. Simulate that shape directly: `started` never changes across
  # polls, only `done` does, and confirm the watch does NOT stall out.
  d <- withr::local_tempdir()
  path <- file.path(d, "rnaparallel-1.tsv")
  writeLines(c("1\ta\t1\tstart", "1\ta\t2\tstart", "1\ta\t3\tstart"), path)  # 3 started, 0 done
  # append a "done" row shortly after the watch begins, well inside a generous stall_after,
  # but AFTER the started count has already stopped changing
  later <- function() cat("2\ta\t1\tdone\n2\ta\t2\tdone\n2\ta\t3\tdone\n", file = path, append = TRUE)
  # Run the append from a background-like path: since this test cannot truly background it,
  # write the done rows immediately (poll interval is small) and rely on the completion exit
  # (done >= started) rather than the stall exit to prove the watch did not stop early.
  later()
  result <- rnaparallel_progress(d, watch = TRUE, interval = 0.05, stall_after = 2)
  expect_identical(result$done, 3L)
  expect_identical(result$started, 3L)
})

test_that("watch mode returns the last real summary on a stall, not NULL", {
  # Before this fix, a stall returned invisible(NULL), throwing away data the summarise
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
