skip_no_limma <- function() testthat::skip_if_not_installed("limma")

rbe_fixture <- function(G = 400L, N = 24L, seed = 7L) {
  set.seed(seed)
  y <- matrix(stats::rnorm(G * N, 8, 2), G, N,
              dimnames = list(sprintf("g%04d", seq_len(G)), sprintf("s%02d", seq_len(N))))
  list(y = y,
       batch = factor(rep(1:4, length.out = N)),
       batch2 = factor(rep(1:2, length.out = N)),
       covar = cbind(age = stats::rnorm(N)),
       group = factor(rep(c("A", "B"), length.out = N)))
}

test_that("removeBatchEffect_parallel is identical on every argument path", {
  skip_no_limma()
  withr::local_options(list(combat.min.ls.cells = 0))
  d <- rbe_fixture()
  des <- stats::model.matrix(~ d$group)
  expect_identical(removeBatchEffect_parallel(d$y, batch = d$batch, workers = 2L),
                   limma::removeBatchEffect(d$y, batch = d$batch))
  expect_identical(removeBatchEffect_parallel(d$y, batch = d$batch, design = des, workers = 2L),
                   limma::removeBatchEffect(d$y, batch = d$batch, design = des))
  expect_identical(
    suppressWarnings(removeBatchEffect_parallel(d$y, batch = d$batch, batch2 = d$batch2,
                                                workers = 2L)),
    suppressWarnings(limma::removeBatchEffect(d$y, batch = d$batch, batch2 = d$batch2)))
  expect_identical(removeBatchEffect_parallel(d$y, batch = d$batch, covariates = d$covar,
                                              workers = 2L),
                   limma::removeBatchEffect(d$y, batch = d$batch, covariates = d$covar))
  expect_identical(
    suppressWarnings(removeBatchEffect_parallel(d$y, batch = d$batch, group = d$group,
                                                workers = 2L)),
    suppressWarnings(limma::removeBatchEffect(d$y, batch = d$batch, group = d$group)))
  # the early return that never reaches lmFit at all
  expect_identical(removeBatchEffect_parallel(d$y, workers = 2L),
                   limma::removeBatchEffect(d$y))
})

test_that("removeBatchEffect_parallel is identical across chunk layouts and backends", {
  skip_no_limma()
  withr::local_options(list(combat.min.ls.cells = 0))
  d <- rbe_fixture()
  ref <- limma::removeBatchEffect(d$y, batch = d$batch)
  for (k in c(1L, 2L, 3L, 7L, 64L)) {
    expect_identical(removeBatchEffect_parallel(d$y, batch = d$batch, workers = 2L, chunks = k),
                     ref, info = paste("chunks", k))
  }
  for (b in c("serial", "mclapply")) {
    expect_identical(removeBatchEffect_parallel(d$y, batch = d$batch, workers = 2L,
                                                parallel_backend = b), ref, info = b)
  }
})

test_that("removeBatchEffect_parallel actually dispatches", {
  skip_no_limma()
  withr::local_options(list(combat.min.ls.cells = 0))
  d <- rbe_fixture()
  # identical() alone cannot tell a working parallel layer from a dead one: the original called
  # unchanged returns exactly the same matrix at serial pace
  n <- 0L
  spy <- function(idx, f, w) { n <<- n + 1L; lapply(idx, f) }
  got <- removeBatchEffect_parallel(d$y, batch = d$batch, workers = 2L, chunks = 4L,
                                    parallel_backend = spy)
  expect_identical(got, limma::removeBatchEffect(d$y, batch = d$batch))
  expect_gt(n, 0L)
})

test_that("removeBatchEffect_parallel refuses what lmFit_parallel refuses", {
  skip_no_limma()
  withr::local_options(list(combat.min.ls.cells = 0))
  d <- rbe_fixture()
  # method = "robust" reshapes nothing but is refused upstream because a row split changes the
  # robust weights; ndups >= 2 pairs different genes. Both arrive here through `...`.
  expect_error(removeBatchEffect_parallel(d$y, batch = d$batch, method = "robust", workers = 2L))
})

test_that("the lmFit rebind inside removeBatchEffect is gated for reachability", {
  skip_no_limma()
  d <- rbe_fixture()
  # a limma that wrote limma::lmFit(...) would make the rebind a no-op: correct output, original
  # speed, invisible to every equivalence assertion in this file
  drifted <- limma::removeBatchEffect
  txt <- paste(deparse(body(drifted)), collapse = "\n")
  skip_if(!grepl("lmFit(", txt, fixed = TRUE), "original body is not shaped as expected")
  body(drifted) <- parse(text = sub("lmFit(", "limma::lmFit(", txt, fixed = TRUE))[[1L]]
  expect_error(removeBatchEffect_parallel(d$y, batch = d$batch, workers = 2L,
                                          backend = drifted),
               "no longer calls lmFit")
})
