## combat.mem.chunk.cells: an explicit per-chunk cell-count ceiling, raising the chunk count
## past what workers/chunks requested so no single chunk exceeds the budget. Opt-in (NULL by
## default = no effect), because converting "cells" to real bytes is workload-dependent and
## this package has no measured per-companion multiplier to derive one from live RAM safely.

test_that("unset (the default) changes nothing: identical chunk layout to before", {
  withr::local_options(combat.mem.chunk.cells = NULL)
  a <- rnaparallel:::combat_row_chunks(100L, workers = 4L, chunks = 4L, ncol = 50L)
  b <- rnaparallel:::combat_row_chunks(100L, workers = 4L, chunks = 4L)   # no ncol at all
  expect_identical(a, b)
  expect_length(a, 4L)
})

test_that("raises the chunk count when a chunk would exceed the cell budget", {
  # 100 rows x 50 cols = 5,000 cells/row-block at chunks=4 (25 rows/chunk -> 1,250 cells/chunk).
  # A 500-cell budget needs at least ceiling(100*50/500) = 10 chunks.
  withr::local_options(combat.mem.chunk.cells = 500)
  idx <- rnaparallel:::combat_row_chunks(100L, workers = 4L, chunks = 4L, ncol = 50L)
  expect_gte(length(idx), 10L)
  # every row still covered exactly once, regardless of how many chunks that took
  expect_identical(sort(unlist(idx)), 1:100)
})

test_that("never LOWERS the chunk count below what was requested", {
  # A generous budget must not merge chunks back down; chunks is a floor operation only.
  withr::local_options(combat.mem.chunk.cells = 1e12)
  idx <- rnaparallel:::combat_row_chunks(100L, workers = 4L, chunks = 8L, ncol = 50L)
  expect_length(idx, 8L)
})

test_that("has no effect when ncol is not supplied, even with the option set", {
  withr::local_options(combat.mem.chunk.cells = 1)   # would force max chunking if it applied
  idx <- rnaparallel:::combat_row_chunks(100L, workers = 4L, chunks = 4L)   # ncol omitted
  expect_length(idx, 4L)
})

test_that("respects min_rows even when the cell budget would want more chunks", {
  # budget wants far more than 50 chunks (100 rows / 2 min_rows), but min_rows still wins
  withr::local_options(combat.mem.chunk.cells = 1)
  idx <- rnaparallel:::combat_row_chunks(100L, workers = 4L, chunks = 4L, ncol = 50L, min_rows = 2L)
  expect_lte(length(idx), 50L)
  expect_identical(sort(unlist(idx)), 1:100)
})

test_that("a garbage combat.mem.chunk.cells is ignored, not an error and not a crash", {
  # This option is read with suppressWarnings(as.numeric()), the same tolerant pattern the
  # other opt-in memory options use, since a raw dispatch-path option here should not be
  # able to take down a whole run over a typo; NA/non-numeric falls back to "no effect".
  withr::local_options(combat.mem.chunk.cells = "lots")
  expect_no_error(idx <- rnaparallel:::combat_row_chunks(100L, workers = 4L, chunks = 4L, ncol = 50L))
  expect_length(idx, 4L)
})

test_that("a real companion call is bit-identical with the budget forcing extra chunks", {
  # The actual promise: this feature must never change a single output value, only how the
  # SAME work is split. Force it to a tiny budget (guaranteed to raise chunk count well past
  # what was asked) and confirm the answer is untouched.
  skip_if_not_installed("limma")
  skip_if_not_installed("edgeR")
  set.seed(9)
  G <- 200L; S <- 12L
  lib <- round(exp(stats::rnorm(S, log(2e6), 0.3)))
  p   <- 2^pmax(stats::rnorm(G, 3, 2), -1); p <- p / sum(p)
  cts <- matrix(stats::rnbinom(G * S, mu = outer(p, lib), size = 10), G, S)
  rownames(cts) <- sprintf("g%04d", seq_len(G))
  colnames(cts) <- sprintf("s%02d", seq_len(S))

  ns <- asNamespace("edgeR")
  nm <- if (exists("normLibSizes.default", envir = ns, inherits = FALSE)) "normLibSizes"
        else "calcNormFactors"
  ref <- get(nm, envir = ns, inherits = FALSE)(cts)

  withr::local_options(combat.mem.chunk.cells = 1)   # forces maximum chunking
  got <- calcNormFactors_parallel(cts, workers = 2L, chunks = 1L)
  expect_identical(got, ref)
})
