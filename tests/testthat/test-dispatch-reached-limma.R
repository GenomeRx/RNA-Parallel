## identical() cannot tell a working parallel layer from a wrapper that forwards to limma and
## does nothing: the forwarder passes every assertion in test-equivalence-limma.R. These
## tests exist to prove the split is load-bearing, and the first one cannot be faked.

skip_if_no_limma <- function() {
  testthat::skip_if_not_installed("limma")
  testthat::skip_if_not_installed("edgeR")
  testthat::skip_if_not_installed("statmod")
}

fixture <- function(G = 200L, S = 12L, seed = 5L) {
  set.seed(seed)
  grp <- factor(rep(c("A", "B"), each = S / 2))
  lib <- round(exp(stats::rnorm(S, log(2e6), 0.3)))
  p   <- 2^pmax(stats::rnorm(G, 3, 2), -1); p <- p / sum(p)
  cts <- matrix(stats::rnbinom(G * S, mu = outer(p, lib), size = 10), G, S)
  rownames(cts) <- sprintf("g%04d", seq_len(G))
  colnames(cts) <- sprintf("s%02d", seq_len(S))
  des <- stats::model.matrix(~ grp)
  list(counts = cts, design = des, block = factor(rep(seq_len(S / 2), each = 2)),
       v = limma::voom(edgeR_norm(edgeR::DGEList(cts)), des))
}

## A backend that runs the work and then corrupts one chunk. If the entry point actually
## dispatches through it, the corruption reaches the result. A silent forwarder returns the
## vendor's answer untouched and this test fails, which is the whole point.
poison_chunk <- function(field) {
  function(idx, f, workers) {
    r <- lapply(idx, f)
    v <- r[[1L]]$value
    if (is.list(v) && !is.null(v[[field]])) v[[field]][1L] <- v[[field]][1L] + 1
    r[[1L]]$value <- v
    r
  }
}

counting <- function() {
  n <- 0L
  list(backend = function(idx, f, workers) { n <<- n + 1L; lapply(idx, f) },
       count = function() n)
}

test_that("lmFit_parallel actually dispatches, proved by poisoning a chunk", {
  skip_if_no_limma()
  d <- fixture()
  ref <- limma::lmFit(d$v, d$design)
  got <- lmFit_parallel(d$v, d$design, workers = 2L, chunks = 4L,
                        parallel_backend = poison_chunk("sigma"))
  expect_false(identical(got, ref))
  expect_false(identical(got$sigma, ref$sigma))
  # and only the poisoned chunk moved: everything else still matches
  expect_equal(sum(got$sigma != ref$sigma), 1L)
})

test_that("duplicateCorrelation_parallel actually dispatches", {
  skip_if_no_limma()
  d <- fixture()
  ref <- limma::duplicateCorrelation(d$v, d$design, block = d$block)
  spy <- counting()
  got <- duplicateCorrelation_parallel(d$v, d$design, block = d$block, workers = 2L,
                                       chunks = 4L, parallel_backend = spy$backend)
  expect_identical(got, ref)
  expect_identical(spy$count(), 1L)
})

test_that("calcNormFactors_parallel actually dispatches", {
  skip_if_no_limma()
  d <- fixture()
  spy <- counting()
  got <- calcNormFactors_parallel(d$counts, workers = 2L, chunks = 4L,
                                  parallel_backend = spy$backend)
  expect_identical(got, edgeR_norm(d$counts))
  # Two dispatches, not one: the reference-column selection ranks across genes and the
  # per-sample factor loop are separate stages, and setup-parallel.R zeroes both gates.
  expect_identical(spy$count(), 2L)

  # The DGEList path rebinds a DIFFERENT symbol -- the generic edgeR's DGEList method calls --
  # so it needs its own dispatch count. Point that rebind at the wrong name and the vendor
  # runs fully serially, returns identical() output, and every equivalence assertion in this
  # suite stays green: measured 2 dispatches on a matrix and 0 on a DGEList.
  spy2 <- counting()
  got2 <- calcNormFactors_parallel(edgeR::DGEList(d$counts), workers = 2L, chunks = 4L,
                                   parallel_backend = spy2$backend)
  expect_identical(got2, edgeR_norm(edgeR::DGEList(d$counts)))
  expect_identical(spy2$count(), 2L)
})

test_that("the backend gate refuses a limma whose rebind target moved", {
  skip_if_no_limma()
  d <- fixture()
  # A backend that no longer calls lm.series as a bare symbol. Rebinding cannot reach it, so
  # the companion would run serially while every equivalence test still passed.
  fake <- function(object, design = NULL, ndups = NULL, spacing = NULL, block = NULL,
                   correlation, weights = NULL, method = "ls", ...) {
    limma::lmFit(object, design)
  }
  expect_error(lmFit_parallel(d$v, d$design, workers = 2L, backend = fake),
               "no longer calls|bare symbol")
})

test_that("the backend gate refuses a namespace-qualified call", {
  skip_if_no_limma()
  d <- fixture()
  # all.names() would flatten limma::lm.series to include lm.series and wave this through.
  # The call-head walk is what refuses it.
  fake <- function(object, design = NULL, ndups = NULL, spacing = NULL, block = NULL,
                   correlation, weights = NULL, method = "ls", ...) {
    limma::lm.series(object, design = design)
  }
  expect_error(lmFit_parallel(d$v, d$design, workers = 2L, backend = fake),
               "no longer calls|bare symbol")
})

test_that("a backend returning the wrong number of chunks is refused, not reassembled", {
  skip_if_no_limma()
  d <- fixture()
  short <- function(idx, f, workers) lapply(idx, f)[-1L]
  expect_error(lmFit_parallel(d$v, d$design, workers = 2L, chunks = 4L,
                              parallel_backend = short),
               "length|chunk|result")
})

test_that("blocks that took different branches are never assembled into one result", {
  skip_if_no_limma()
  d <- fixture()
  # rp_invariant's component-set backstop. Two blocks on opposite sides of NoProbeWts return
  # different names(), and assembling them mixes colMeans sigma with per-gene mean sigma.
  slow <- list(coefficients = matrix(0, 2, 2), stdev.unscaled = matrix(0, 2, 2),
               sigma = c(1, 1), df.residual = c(1, 1), cov.coefficients = matrix(0, 2, 2),
               pivot = 1:2, rank = 2L)
  fast <- c(slow, list(qr = list(), assign = 1:2))
  expect_error(rp_invariant(list(slow, fast), c("rank", "pivot"), "test"),
               "different component sets|different branches")
})

test_that("an unreachable option value cannot silently disable the split", {
  skip_if_no_limma()
  d <- fixture()
  ref <- limma::lmFit(d$v, d$design)
  # Above the gate the split runs; below it the vendor is called whole. Both are identical(),
  # which is exactly why the gate needs its own assertion rather than being inferred.
  # voom weights carry no arrayweights attribute, so this input takes the slow branch and
  # is gated by combat.min.cells rather than combat.min.ls.cells. Setting the wrong one of
  # the two looks like a passing gate while nothing is gated, so both are named here.
  withr::local_options(combat.min.cells = 0, combat.min.ls.cells = 0)
  expect_false(identical(
    lmFit_parallel(d$v, d$design, workers = 2L, chunks = 4L,
                   parallel_backend = poison_chunk("sigma")), ref))
  withr::local_options(combat.min.cells = Inf, combat.min.ls.cells = Inf)
  expect_identical(
    lmFit_parallel(d$v, d$design, workers = 2L, chunks = 4L,
                   parallel_backend = poison_chunk("sigma")), ref)
})

test_that("duplicateCorrelation_parallel results reach the returned object", {
  skip_if_no_limma()
  # A counting spy proves only that the backend was called. This proves its output is used.
  d <- fixture()
  ref <- limma::duplicateCorrelation(d$v, d$design, block = d$block)
  poison <- function(idx, f, workers) {
    r <- lapply(idx, f)
    v <- r[[1L]]$value
    v$value$atanh.correlations[1L] <- v$value$atanh.correlations[1L] + 1
    r[[1L]]$value <- v
    r
  }
  got <- duplicateCorrelation_parallel(d$v, d$design, block = d$block, workers = 2L,
                                       chunks = 4L, parallel_backend = poison)
  expect_false(identical(got$atanh.correlations, ref$atanh.correlations))
})
