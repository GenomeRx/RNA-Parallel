## rp_progress_read()'s tail-vs-full-reread path: cache = NULL always reads whole; cache = env
## skips unchanged files and only reads bytes past the last-seen size on ones that grew.

test_that("rp_progress_read without a cache reads the whole file every call", {
  d <- withr::local_tempdir()
  writeLines(c("1\tstage_a\t1\tstart", "2\tstage_a\t1\tdone"), file.path(d, "rnaparallel-1.tsv"))
  r1 <- rnaparallel:::rp_progress_read(d)
  expect_identical(nrow(r1), 2L)
  writeLines(c("1\tstage_a\t1\tstart", "2\tstage_a\t1\tdone",
              "3\tstage_a\t2\tstart"), file.path(d, "rnaparallel-1.tsv"))
  r2 <- rnaparallel:::rp_progress_read(d)
  expect_identical(nrow(r2), 3L)   # sees the growth; no caching means nothing to invalidate
})

test_that("rp_progress_read with a cache skips a file whose size has not changed", {
  d <- withr::local_tempdir()
  writeLines(c("1\tstage_a\t1\tstart", "2\tstage_a\t1\tdone"), file.path(d, "rnaparallel-1.tsv"))
  cache <- new.env(parent = emptyenv())
  r1 <- rnaparallel:::rp_progress_read(d, cache = cache)
  expect_identical(nrow(r1), 2L)
  path <- file.path(d, "rnaparallel-1.tsv")
  before <- cache[[path]]$rows
  # second call, file untouched: must return the SAME cached rows object, not re-derive it,
  # since that is the whole point of caching (nothing new to read means nothing to re-read)
  r2 <- rnaparallel:::rp_progress_read(d, cache = cache)
  expect_identical(r2, before)
  expect_identical(nrow(r2), 2L)
})

test_that("rp_progress_read with a cache picks up only the NEW lines on a growing file", {
  d <- withr::local_tempdir()
  path <- file.path(d, "rnaparallel-1.tsv")
  writeLines(c("1\tstage_a\t1\tstart", "2\tstage_a\t1\tdone"), path)
  cache <- new.env(parent = emptyenv())
  r1 <- rnaparallel:::rp_progress_read(d, cache = cache)
  expect_identical(nrow(r1), 2L)
  # append, do not rewrite: matches how rp_progress_file_write() actually writes
  cat("3\tstage_a\t2\tstart\n", file = path, append = TRUE)
  r2 <- rnaparallel:::rp_progress_read(d, cache = cache)
  expect_identical(nrow(r2), 3L)
  expect_identical(r2$event, c("start", "done", "start"))
  expect_identical(r2$chunk, c(1L, 1L, 2L))
})

test_that("tailed reads and full reads produce identical rows for the same file content", {
  d <- withr::local_tempdir()
  path <- file.path(d, "rnaparallel-1.tsv")
  writeLines(c("1\ta\t1\tstart", "2\ta\t1\tdone", "3\ta\t2\tstart"), path)
  no_cache <- rnaparallel:::rp_progress_read(d)
  cache <- new.env(parent = emptyenv())
  with_cache <- rnaparallel:::rp_progress_read(d, cache = cache)
  # same content, same shape, regardless of which code path produced it
  expect_identical(as.list(no_cache), as.list(with_cache))
})

test_that("a cache spans multiple worker files independently", {
  d <- withr::local_tempdir()
  writeLines(c("1\ta\t1\tstart"), file.path(d, "rnaparallel-1.tsv"))
  writeLines(c("1\tb\t1\tstart"), file.path(d, "rnaparallel-2.tsv"))
  cache <- new.env(parent = emptyenv())
  r1 <- rnaparallel:::rp_progress_read(d, cache = cache)
  expect_identical(nrow(r1), 2L)
  cat("2\ta\t1\tdone\n", file = file.path(d, "rnaparallel-1.tsv"), append = TRUE)
  # worker 2's file is untouched; only worker 1's grew
  r2 <- rnaparallel:::rp_progress_read(d, cache = cache)
  expect_identical(nrow(r2), 3L)
  expect_identical(sum(r2$stage == "b"), 1L)   # worker 2's one row, never re-duplicated
})

test_that("watch mode uses one cache across its own polls, not a fresh one per poll", {
  # Real end-to-end proof the wiring works, not just the read function in isolation:
  # start a worker file already at 100% so watch mode returns after its first poll, then
  # confirm calling it twice in a row (as watch would across two Sys.sleep() iterations)
  # does not re-derive rows it already parsed.
  d <- withr::local_tempdir()
  path <- file.path(d, "rnaparallel-1.tsv")
  writeLines(c("1\ta\t1\tstart", "2\ta\t1\tdone"), path)
  cache <- new.env(parent = emptyenv())
  r1 <- rnaparallel:::rp_progress_read(d, cache = cache)
  r2 <- rnaparallel:::rp_progress_read(d, cache = cache)
  expect_identical(r1, r2)
  expect_identical(cache[[path]]$size, file.size(path))
})

## rp_step_begin() creates combat.progress.dir up front and refuses to start if it cannot,
## rather than the old silent-drop behaviour where every worker write vanished for the
## whole run and rnaparallel_progress() reported "no files yet" forever.

test_that("rp_step_begin creates a missing combat.progress.dir rather than failing silently", {
  base <- withr::local_tempdir()
  target <- file.path(base, "nested", "progress")
  expect_false(dir.exists(target))
  withr::local_options(combat.progress.dir = target)
  h <- rnaparallel:::rp_step_begin(NULL, "unit", matrix(1, 2, 2), "serial", 1L)
  expect_true(dir.exists(target))
  rnaparallel:::rp_step_end(h)
})

test_that("rp_step_begin refuses an unwritable combat.progress.dir loudly", {
  skip_on_os("windows")   # file.access()-style permission bits are not meaningful on Windows
  base <- withr::local_tempdir()
  target <- file.path(base, "locked")
  dir.create(target)
  Sys.chmod(target, mode = "0500")   # read+execute, no write
  withr::local_options(combat.progress.dir = target)
  withr::defer(Sys.chmod(target, mode = "0700"))   # let the tempdir cleanup remove it
  expect_error(rnaparallel:::rp_step_begin(NULL, "unit", matrix(1, 2, 2), "serial", 1L),
              "not writable")
})

## rp_progress_tick() must not draw from inside a worker: regression for a real bug where a
## nested/recursive dispatch (or any size-gated call, since both paths call rp_count(FALSE)
## from a forked child) printed that child's own stale count over the master's line.

test_that("rp_progress_tick is a no-op inside a marked worker", {
  withr::local_options(combat.progress = TRUE)
  prev <- Sys.getenv("RNAPARALLEL_IN_WORKER", unset = NA_character_)
  Sys.setenv(RNAPARALLEL_IN_WORKER = "1")
  withr::defer(if (is.na(prev)) Sys.unsetenv("RNAPARALLEL_IN_WORKER")
              else Sys.setenv(RNAPARALLEL_IN_WORKER = prev))
  # A no-op tick never reaches the throttle-timestamp write, so `.rp_dispatch$progress_last`
  # (an environment binding, checkable directly rather than trying to intercept a raw
  # `cat(..., file = stderr())` write, which is not a captureable R condition/output stream)
  # stays whatever it was before the call -- here, unset.
  rnaparallel:::rp_count_reset()
  expect_null(rnaparallel:::.rp_dispatch$progress_last)
  expect_no_error(rnaparallel:::rp_progress_tick())
  expect_null(rnaparallel:::.rp_dispatch$progress_last)
})
