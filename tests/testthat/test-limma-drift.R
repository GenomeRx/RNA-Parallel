## Each test passes the original or a drifted copy of it and asserts identical() on the whole result.

skip_g4 <- function() {
  testthat::skip_if_not_installed("limma")
  testthat::skip_if_not_installed("edgeR")
  testthat::skip_if_not_installed("statmod")
}

outcome <- function(f) {
  w <- character()
  v <- withCallingHandlers(
    tryCatch(f(), error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  list(value = v, warnings = w)
}

test_that("the model.matrix memo keys on every variable the formula reads", {
  skip_g4()
  f <- limma::duplicateCorrelation
  drifted <- f
  body(drifted) <- do.call(substitute, list(body(f), list(A = quote(B))))
  environment(drifted) <- environment(f)
  skip_if(!grepl("~0 + B", paste(deparse(body(drifted)), collapse = ""), fixed = TRUE),
          "limma's per-gene loop no longer builds ~0 + A")
  set.seed(31)
  M <- matrix(stats::rnorm(60 * 24), 60, 24)
  des <- cbind(1, rep(0:1, each = 12))
  blk <- factor(rep(1:12, each = 2))
  M[1, 3] <- NA
  M[2, 9] <- NA
  expect_identical(
    duplicateCorrelation_parallel(M, des, block = blk, workers = 1L, chunks = 1L,
                                  parallel_backend = "serial", backend = drifted),
    drifted(M, des, block = blk))
  expect_identical(
    duplicateCorrelation_parallel(M, des, block = blk, workers = 2L, chunks = 3L,
                                  parallel_backend = "serial", backend = drifted),
    drifted(M, des, block = blk))
  expect_identical(
    duplicateCorrelation_parallel(M, des, block = blk, workers = 2L, chunks = 3L,
                                  parallel_backend = "serial"),
    limma::duplicateCorrelation(M, des, block = blk))
})

test_that("a per-gene field the row fitter adds is never cut to block 1's genes", {
  skip_g4()
  lns <- asNamespace("limma")
  e <- new.env(parent = lns)
  e$lm.series <- function(M, design = NULL, ndups = 1, spacing = 1, weights = NULL) {
    fit <- get("lm.series", envir = lns)(M, design = design, ndups = ndups, spacing = spacing,
                                          weights = weights)
    fit$rowsum <- rowSums(as.matrix(M))
    fit
  }
  e$gls.series <- function(M, design = NULL, ndups = 2, spacing = 1, block = NULL,
                           correlation = NULL, weights = NULL, ...) {
    fit <- get("gls.series", envir = lns)(M, design = design, ndups = ndups, spacing = spacing,
                                           block = block, correlation = correlation,
                                           weights = weights, ...)
    fit$rowsum <- rowSums(as.matrix(M))
    fit
  }
  lm2 <- limma::lmFit
  environment(lm2) <- e
  set.seed(2)
  y <- matrix(stats::rnorm(200 * 8), 200, 8)
  des <- cbind(1, rep(0:1, 4))
  blk <- rep(1:4, each = 2)
  n <- 0L
  spy <- function(idx, f, workers) { n <<- n + 1L; lapply(idx, f) }

  expect_identical(lmFit_parallel(y, des, workers = 2L, chunks = 4L, parallel_backend = spy,
                                  backend = lm2),
                   lm2(y, des))
  expect_gt(n, 0L)
  expect_identical(lmFit_parallel(y, des, block = blk, correlation = 0.2, workers = 2L,
                                  chunks = 4L, parallel_backend = spy, backend = lm2),
                   lm2(y, des, block = blk, correlation = 0.2))
  withr::local_options(combat.timing = TRUE, combat.timing.min = 0)
  expect_message(lmFit_parallel(y, des, workers = 2L, chunks = 4L, parallel_backend = spy,
                                backend = lm2),
                 "stood down")
})

test_that("a degenerate block design returns the original's early value whatever trim is", {
  skip_g4()
  set.seed(2)
  M <- matrix(stats::rnorm(60), 10, 6)
  for (pb in c("serial", "mclapply")) {
    expect_identical(
      outcome(function() duplicateCorrelation_parallel(M, block = 1:6, trim = NA, workers = 2L,
                                                       parallel_backend = pb)),
      outcome(function() limma::duplicateCorrelation(M, block = 1:6, trim = NA)),
      info = pb)
  }
  blk <- rep(1:3, each = 2)
  des <- cbind(1, stats::model.matrix(~ factor(blk))[, -1])
  expect_identical(
    outcome(function() duplicateCorrelation_parallel(M, design = des, block = blk, trim = "x",
                                                     workers = 2L, parallel_backend = "serial")),
    outcome(function() limma::duplicateCorrelation(M, design = des, block = blk, trim = "x")))
# A usable trim keeps one warning, and a bad trim on a real block design errors the same way.
  expect_identical(
    outcome(function() duplicateCorrelation_parallel(M, block = 1:6, workers = 2L,
                                                     parallel_backend = "serial")),
    outcome(function() limma::duplicateCorrelation(M, block = 1:6)))
  expect_identical(
    outcome(function() duplicateCorrelation_parallel(M, block = rep(1:3, 2), trim = NA,
                                                     workers = 2L, parallel_backend = "serial")),
    outcome(function() limma::duplicateCorrelation(M, block = rep(1:3, 2), trim = NA)))
})

test_that("removeBatchEffect_parallel pairs a backend with the lmFit in its own environment", {
  skip_g4()
  set.seed(7)
  y <- matrix(stats::rnorm(400 * 24, 8, 2), 400, 24)
  b <- factor(rep(1:4, length.out = 24))
  fake <- new.env(parent = asNamespace("limma"))
  lf <- limma::lmFit
  bd <- body(lf)
  bd[[length(bd)]] <- quote({
    fit$coefficients <- fit$coefficients * 2
    new("MArrayLM", fit)
  })
  body(lf) <- bd
  environment(lf) <- fake
  fake$lmFit <- lf
  rbe2 <- limma::removeBatchEffect
  environment(rbe2) <- fake
  direct <- rbe2(y, batch = b)
  expect_false(identical(direct, limma::removeBatchEffect(y, batch = b)))
  expect_identical(removeBatchEffect_parallel(y, batch = b, backend = rbe2, workers = 2L,
                                              chunks = 3L, parallel_backend = "serial"),
                   direct)
  expect_identical(removeBatchEffect_parallel(y, backend = rbe2, workers = 2L),
                   rbe2(y))
})

test_that("a gated call never expands vector weights to a full matrix", {
  skip_g4()
  lns <- asNamespace("limma")
  n_amw <- 0L
  e <- new.env(parent = lns)
  e$asMatrixWeights <- function(...) {
    n_amw <<- n_amw + 1L
    get("asMatrixWeights", envir = lns)(...)
  }
  lm2 <- limma::lmFit
  environment(lm2) <- e
  set.seed(4)
  G <- 300L
  S <- 12L
  M <- matrix(stats::rnorm(G * S), G, S)
  des <- cbind(1, rep(0:1, S / 2))
  aw <- stats::runif(S, 0.5, 2)
  gw <- stats::runif(G, 0.5, 2)
  check <- function(w, expansions, ...) {
    n_amw <<- 0L
    expect_identical(lmFit_parallel(M, des, weights = w, workers = 2L, chunks = 3L,
                                    backend = lm2, ...),
                     limma::lmFit(M, des, weights = w, ...))
    expect_identical(n_amw, expansions)
  }

  withr::local_options(combat.min.ls.cells = 1e9, combat.min.cells = 1e9)
  check(aw, 0L)
  check(matrix(aw, 1L), 0L)
  check(gw, 0L)
  check(aw, 0L, block = rep(1:6, each = 2), correlation = 0.2)

# At 12 x 12 a length-12 vector is gene weights (asMatrixWeights tests that branch first), so the shut weighted gate applies.
  withr::local_options(combat.min.ls.cells = 0)
  n_amw <- 0L
  Ms <- matrix(stats::rnorm(12 * 12), 12, 12)
  w12 <- stats::runif(12, 0.5, 2)
  expect_identical(lmFit_parallel(Ms, des, weights = w12, workers = 2L, backend = lm2),
                   limma::lmFit(Ms, des, weights = w12))
  expect_identical(n_amw, 0L)

  withr::local_options(combat.min.ls.cells = 0, combat.min.cells = 0, combat.min.wt.genes = 0)
  check(aw, 1L)
  check(gw, 1L)
})

test_that("removeBatchEffect_parallel caps workers once and reports the count it used", {
  skip_g4()
  set.seed(7)
  y <- matrix(stats::rnorm(400 * 24, 8, 2), 400, 24)
  b <- factor(rep(1:4, length.out = 24))
  ref <- limma::removeBatchEffect(y, batch = b)
  rss_calls <- 0L
  local_mocked_bindings(
    rp_mem_available = function() 250 / 0.8,
    rp_mem_rss = function() {
      rss_calls <<- rss_calls + 1L
      100 + 30 * (rss_calls - 1L)
    })
  seen <- integer()
  spy <- function(idx, f, workers) { seen <<- c(seen, workers); lapply(idx, f) }
  warns <- character()
  msgs <- character()
  withr::local_options(combat.timing = TRUE, combat.timing.min = 0, combat.mem.guard = TRUE,
                       combat.mem.divergence = 1)
  got <- withCallingHandlers(
    removeBatchEffect_parallel(y, batch = b, workers = 3L, chunks = 3L, parallel_backend = spy),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_identical(got, ref)
  expect_length(grep("workers need", warns), 1L)
  expect_identical(rss_calls, 1L)
  expect_true(length(seen) > 0L && all(seen == 2L))
  expect_true(any(grepl("custom x2", msgs, fixed = TRUE)))
  withr::local_options(combat.mem.guard = "yes")
  expect_error(removeBatchEffect_parallel(y, workers = 2L), "must be TRUE or FALSE")
})

test_that("a NULL correlation is resolved by the parallel duplicateCorrelation, identically", {
  skip_g4()
  set.seed(12)
  G <- 240L
  S <- 12L
  grp <- factor(rep(c("A", "B"), each = S / 2))
  des <- stats::model.matrix(~ grp)
  blk <- factor(rep(seq_len(S / 2), each = 2))
  cts <- matrix(stats::rnbinom(G * S, mu = 2^stats::runif(G, 2, 9), size = 10), G, S,
                dimnames = list(sprintf("g%04d", seq_len(G)), sprintf("s%02d", seq_len(S))))
  v <- limma::voom(cts, des)
  M <- v$E
  Mna <- M
  Mna[7L, 3L] <- NA

  n <- 0L
  spy <- function(idx, f, workers) { n <<- n + 1L; lapply(idx, f) }
  expect_identical(lmFit_parallel(M, des, block = blk, correlation = NULL, workers = 2L,
                                  chunks = 3L, parallel_backend = spy),
                   limma::lmFit(M, des, block = blk, correlation = NULL))
  expect_identical(n, 2L)

  for (k in c(1L, 2L, 3L, 7L)) {
    for (pb in c("serial", "mclapply")) {
      tag <- paste("chunks", k, pb)
      expect_identical(lmFit_parallel(M, des, block = blk, correlation = NULL, workers = 2L,
                                      chunks = k, parallel_backend = pb),
                       limma::lmFit(M, des, block = blk, correlation = NULL), info = tag)
      expect_identical(lmFit_parallel(v, des, block = blk, correlation = NULL, workers = 2L,
                                      chunks = k, parallel_backend = pb),
                       limma::lmFit(v, des, block = blk, correlation = NULL), info = tag)
      expect_identical(lmFit_parallel(Mna, des, block = blk, correlation = NULL, workers = 2L,
                                      chunks = k, parallel_backend = pb),
                       limma::lmFit(Mna, des, block = blk, correlation = NULL), info = tag)
    }
  }
  expect_identical(lmFit_parallel(M, des, block = blk, correlation = NULL, trim = 0.3,
                                  workers = 2L, chunks = 3L),
                   limma::lmFit(M, des, block = blk, correlation = NULL, trim = 0.3))
# The block is encoded in this design, so the consensus is the original's early zero.
  blk3 <- rep(1:2, each = S / 2)
  expect_identical(
    outcome(function() lmFit_parallel(M, des, block = blk3, correlation = NULL, workers = 2L,
                                      chunks = 3L)),
    outcome(function() limma::lmFit(M, des, block = blk3, correlation = NULL)))
})
