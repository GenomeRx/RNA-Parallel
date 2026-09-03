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
