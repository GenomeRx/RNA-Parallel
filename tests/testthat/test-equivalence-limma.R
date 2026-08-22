## The limma and edgeR companions return identical() output, or they are worthless.
## Assertions are on the WHOLE returned object, never on selected fields: the component set
## itself varies between limma's two lm.series branches ($qr and $assign are present on one
## and absent on the other), and dimnames carry shapes that values alone cannot show.

skip_if_no_limma <- function() {
  testthat::skip_if_not_installed("limma")
  testthat::skip_if_not_installed("edgeR")
  testthat::skip_if_not_installed("statmod")
}

sim <- function(G = 200L, S = 12L, seed = 1L) {
  set.seed(seed)
  grp <- factor(rep(c("A", "B"), each = S / 2))
  lib <- round(exp(stats::rnorm(S, log(2e6), 0.3)))
  p   <- 2^pmax(stats::rnorm(G, 3, 2), -1); p <- p / sum(p)
  cts <- matrix(stats::rnbinom(G * S, mu = outer(p, lib), size = 10), G, S)
  rownames(cts) <- sprintf("g%04d", seq_len(G))
  colnames(cts) <- sprintf("s%02d", seq_len(S))
  list(counts = cts, group = grp, design = stats::model.matrix(~ grp),
       block = factor(rep(seq_len(S / 2), each = 2)))
}

test_that("calcNormFactors_parallel is identical on every method", {
  skip_if_no_limma()
  d <- sim()
  for (m in c("TMM", "TMMwsp", "RLE", "upperquartile", "none")) {
    expect_identical(calcNormFactors_parallel(d$counts, method = m, workers = 2L),
                     edgeR::calcNormFactors(d$counts, method = m),
                     info = m)
  }
})

test_that("calcNormFactors_parallel is identical for a DGEList and across chunk layouts", {
  skip_if_no_limma()
  d <- sim()
  expect_identical(calcNormFactors_parallel(edgeR::DGEList(d$counts), workers = 2L),
                   edgeR::calcNormFactors(edgeR::DGEList(d$counts)))
  ref <- edgeR::calcNormFactors(d$counts)
  for (k in c(1L, 2L, 3L, 4L, 8L)) {
    expect_identical(calcNormFactors_parallel(d$counts, workers = 2L, chunks = k), ref,
                     info = paste("chunks", k))
  }
})

test_that("calcNormFactors_parallel survives all-zero gene rows", {
  skip_if_no_limma()
  d <- sim(); d$counts[c(3L, 17L, 88L), ] <- 0L
  expect_identical(calcNormFactors_parallel(d$counts, workers = 2L),
                   edgeR::calcNormFactors(d$counts))
})

test_that("lmFit_parallel is identical on both lm.series branches", {
  skip_if_no_limma()
  d <- sim()
  v <- limma::voom(edgeR::calcNormFactors(edgeR::DGEList(d$counts)), d$design)

  # slow branch: voom's probe weights carry no arrayweights attribute
  expect_identical(lmFit_parallel(v, d$design, workers = 2L), limma::lmFit(v, d$design))
  # fast branch: no weights at all
  expect_identical(lmFit_parallel(v$E, d$design, workers = 2L), limma::lmFit(v$E, d$design))
  # fast branch with array weights
  aw <- limma::arrayWeights(v$E, d$design)
  expect_identical(lmFit_parallel(v$E, d$design, weights = aw, workers = 2L),
                   limma::lmFit(v$E, d$design, weights = aw))
})

test_that("lmFit_parallel is identical across chunk layouts", {
  skip_if_no_limma()
  d <- sim()
  v   <- limma::voom(edgeR::calcNormFactors(edgeR::DGEList(d$counts)), d$design)
  ref <- limma::lmFit(v, d$design)
  for (k in c(1L, 2L, 3L, 4L, 7L, 64L)) {
    expect_identical(lmFit_parallel(v, d$design, workers = 2L, chunks = k), ref,
                     info = paste("chunks", k))
  }
})

test_that("lmFit_parallel keeps the list(NULL, NULL) dimnames shape rbind would drop", {
  skip_if_no_limma()
  # An unnamed matrix with an unnamed design. gls.series has no coef.names fallback, so the
  # vendor returns dimnames list(NULL, NULL) while rbind alone collapses that to NULL.
  set.seed(1)
  M <- matrix(stats::rnorm(80), 20, 4)
  des <- cbind(rep(1, 4), c(0, 0, 1, 1)); colnames(des) <- NULL
  w <- matrix(stats::runif(80, 0.5, 2), 20, 4)
  got <- lmFit_parallel(M, des, block = c(1, 1, 2, 2), correlation = 0.3, weights = w,
                        workers = 2L, chunks = 2L)
  want <- limma::lmFit(M, des, block = c(1, 1, 2, 2), correlation = 0.3, weights = w)
  expect_identical(got, want)
  expect_identical(dimnames(got$coefficients), list(NULL, NULL))
})

test_that("lmFit_parallel pins the NoProbeWts branch when a weight sits under the gls punch", {
  skip_if_no_limma()
  # A weight strictly inside (0, 1e-15) is punched to NA by gls.series and not by lm.series.
  # Reading the wrong punch rule put 10 of 40 sigma on the wrong branch, silently.
  set.seed(11); n <- 40L; q <- 8L
  M <- matrix(stats::rnorm(n * q), n, q,
              dimnames = list(paste0("g", seq_len(n)), paste0("s", seq_len(q))))
  des <- cbind(Int = 1, Grp = rep(0:1, each = 4))
  # row 5 is not a block-leading row at chunks = 4, so rp_arrayweights_uniform passes it and
  # the punch rule is what decides. On row 1 the uniformity check refused first and this test
  # passed through the serial fallback without ever exercising what it names.
  W <- matrix(1, n, q, dimnames = dimnames(M)); W[5L, 1L] <- 1e-16
  attr(W, "arrayweights") <- TRUE
  expect_identical(
    lmFit_parallel(M, des, block = rep(1:4, each = 2), correlation = 0.3, weights = W,
                   workers = 2L, chunks = 4L),
    limma::lmFit(M, des, block = rep(1:4, each = 2), correlation = 0.3, weights = W))
})

test_that("lmFit_parallel refuses the paths it cannot split exactly", {
  skip_if_no_limma()
  d <- sim()
  expect_error(lmFit_parallel(d$counts, d$design, method = "robust", workers = 2L),
               "robust|mrlm")
  expect_error(lmFit_parallel(d$counts, d$design, ndups = 2L, workers = 2L),
               "ndups|unwrapdups")
})

test_that("duplicateCorrelation_parallel is identical, including the pooled consensus", {
  skip_if_no_limma()
  d <- sim()
  v <- limma::voom(edgeR::calcNormFactors(edgeR::DGEList(d$counts)), d$design)
  ref <- limma::duplicateCorrelation(v, d$design, block = d$block)
  for (k in c(1L, 2L, 3L, 7L)) {
    expect_identical(
      duplicateCorrelation_parallel(v, d$design, block = d$block, workers = 2L, chunks = k),
      ref, info = paste("chunks", k))
  }
})

test_that("duplicateCorrelation_parallel refuses the ndups path", {
  skip_if_no_limma()
  d <- sim()
  expect_error(duplicateCorrelation_parallel(d$counts, d$design, block = NULL, workers = 2L),
               "block")
})

test_that("every backend gives the same answer", {
  skip_if_no_limma()
  d <- sim()
  v   <- limma::voom(edgeR::calcNormFactors(edgeR::DGEList(d$counts)), d$design)
  ref <- limma::lmFit(v, d$design)
  for (b in c("mclapply", "serial")) {
    expect_identical(lmFit_parallel(v, d$design, workers = 2L, parallel_backend = b), ref,
                     info = b)
  }
})

test_that("the final gene list is identical end to end", {
  skip_if_no_limma()
  # A one-ulp sigma difference is invisible upstream and becomes a different gene list after
  # BH, so the assertion that matters is on topTable, not on the fit.
  d <- sim(G = 500L, S = 16L, seed = 3L)
  cm <- suppressWarnings(limma::makeContrasts(grpB, levels = d$design))

  ref <- limma::topTable(limma::eBayes(limma::contrasts.fit(
    limma::lmFit(limma::voom(edgeR::calcNormFactors(edgeR::DGEList(d$counts)), d$design),
                 d$design), cm)), number = Inf, sort.by = "none")

  par <- limma::topTable(limma::eBayes(limma::contrasts.fit(
    lmFit_parallel(limma::voom(calcNormFactors_parallel(edgeR::DGEList(d$counts), workers = 2L),
                               d$design), d$design, workers = 2L), cm)),
    number = Inf, sort.by = "none")

  expect_identical(ref, par)
})

test_that("the duplicateCorrelation memo engages and stays identical, weighted or not", {
  skip_if_no_limma()
  # The old body-pinned lift refused to install whenever weights were present or any cell was
  # non-finite, which is every voom fixture. The memo keys on the arguments instead, so it
  # installs on all of them and a varying finite-value mask simply misses the cache.
  set.seed(31)
  M <- matrix(rnorm(400 * 24), 400, 24,
              dimnames = list(paste0("g", 1:400), paste0("s", 1:24)))
  des <- cbind(Int = 1, Grp = rep(0:1, each = 12))
  blk <- factor(rep(1:12, each = 2))
  memo <- rnaparallel:::rp_dupcor_memo(limma::duplicateCorrelation,
                                       environment(limma::duplicateCorrelation))
  expect_true(is.function(memo))

  expect_identical(duplicateCorrelation_parallel(M, des, block = blk, workers = 2L),
                   limma::duplicateCorrelation(M, des, block = blk))

  W <- matrix(runif(400 * 24, 0.5, 2), 400, 24)
  expect_identical(duplicateCorrelation_parallel(M, des, block = blk, weights = W,
                                                 workers = 2L),
                   limma::duplicateCorrelation(M, des, block = blk, weights = W))

  Mna <- M; Mna[5L, 3L] <- NA
  expect_identical(duplicateCorrelation_parallel(Mna, des, block = blk, workers = 2L),
                   limma::duplicateCorrelation(Mna, des, block = blk))
})

test_that("the memo rebinds the vendor's own mixedModel2Fit, not a copy of it", {
  skip_if_no_limma()
  # This is the invariant that replaced roughly a hundred lines of body-text pinning: the
  # package holds no statmod source, so the object in the rebind environment must BE the
  # vendor's, differing only in where its La.svd resolves.
  memo <- rnaparallel:::rp_dupcor_memo(limma::duplicateCorrelation,
                                       environment(limma::duplicateCorrelation))
  fast <- environment(memo)$mixedModel2Fit
  expect_identical(body(fast), body(statmod::mixedModel2Fit))
  expect_identical(formals(fast), formals(statmod::mixedModel2Fit))
  # the only difference is the shadowed primitive
  expect_false(identical(environment(fast), environment(statmod::mixedModel2Fit)))
  expect_true(is.function(environment(fast)$La.svd))
})
