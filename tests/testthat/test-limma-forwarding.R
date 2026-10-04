## A limma companion forwards only the arguments its caller supplied, so every default it leaves out is the backend's own.

skip_fwd <- function() {
  testthat::skip_if_not_installed("limma")
  testthat::skip_if_not_installed("statmod")
}

fwd_outcome <- function(f) {
  w <- character()
  v <- NULL
  utils::capture.output(v <- withCallingHandlers(
    tryCatch(f(), error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") },
    message = function(m) invokeRestart("muffleMessage")))
  list(value = v, warnings = w)
}

spy_on <- function(fn, args) {
  rec <- new.env()
  rec$seen <- list()
  e <- new.env(parent = environment(fn))
  e$.rec <- rec
  probe <- str2lang(sprintf(".rec$seen[[length(.rec$seen) + 1L]] <- c(%s)",
                            paste0(args, " = missing(", args, ")", collapse = ", ")))
  b <- as.list(body(fn))
  body(fn) <- as.call(c(b[1L], list(probe), b[-1L]))
  environment(fn) <- e
  list(fn = fn, seen = function() rec$seen)
}

fwd_sim <- function(G = 60L, S = 12L, seed = 11L) {
  set.seed(seed)
  list(y = matrix(stats::rnorm(G * S, 8, 2), G, S,
                  dimnames = list(sprintf("g%03d", seq_len(G)), sprintf("s%02d", seq_len(S)))),
       design = cbind(Intercept = 1, Group = rep(0:1, S / 2)),
       block = factor(rep(seq_len(S / 2), each = 2)),
       batch = factor(rep(1:3, length.out = S)))
}

same_defaults <- function(original, companion, skip = character()) {
  vf <- formals(original)
  pf <- formals(companion)
  expect_true(all(names(vf) %in% names(pf)))
  for (a in setdiff(names(vf), skip)) {
    expect_identical(deparse(pf[[a]]), deparse(vf[[a]]), info = paste("default differs for", a))
  }
}


test_that("every duplicateCorrelation argument exists on the companion with the same default", {
  skip_fwd()
  same_defaults(limma::duplicateCorrelation, duplicateCorrelation_parallel)
})

test_that("every lmFit argument exists on the companion with the same default", {
  skip_fwd()
  same_defaults(limma::lmFit, lmFit_parallel)
})

test_that("every removeBatchEffect argument exists on the companion with the same default", {
  skip_fwd()
# The design default is compared only where limma's is NULL, because older limma uses matrix(1, ncol(x), 1) and the companion forwards only supplied arguments.
  same_defaults(limma::removeBatchEffect, removeBatchEffect_parallel,
                skip = if (!is.null(formals(limma::removeBatchEffect)$design)) "design")
})


test_that("duplicateCorrelation_parallel hands a block only the arguments the caller supplied", {
  skip_fwd()
  d <- fwd_sim()
  spy <- spy_on(limma::duplicateCorrelation, c("ndups", "spacing", "trim"))
  got <- duplicateCorrelation_parallel(d$y, d$design, block = d$block, workers = 2L,
                                       chunks = 3L, parallel_backend = "serial",
                                       backend = spy$fn)
  expect_identical(got, limma::duplicateCorrelation(d$y, d$design, block = d$block))
  seen <- spy$seen()
  expect_length(seen, 3L)
  for (s in seen) expect_identical(s, c(ndups = TRUE, spacing = TRUE, trim = TRUE))
})

test_that("duplicateCorrelation_parallel pools with the backend's own trim default", {
  skip_fwd()
  d <- fwd_sim()
  drifted <- limma::duplicateCorrelation
  formals(drifted)$trim <- 0.4
  ref <- drifted(d$y, d$design, block = d$block)
  expect_false(identical(ref, limma::duplicateCorrelation(d$y, d$design, block = d$block)))
  expect_identical(
    duplicateCorrelation_parallel(d$y, d$design, block = d$block, workers = 2L, chunks = 3L,
                                  parallel_backend = "serial", backend = drifted),
    ref)
  expect_identical(
    duplicateCorrelation_parallel(d$y, d$design, block = d$block, trim = 0.15, workers = 2L,
                                  chunks = 3L, parallel_backend = "serial", backend = drifted),
    drifted(d$y, d$design, block = d$block, trim = 0.15))
})


test_that("lmFit_parallel hands the backend only the arguments the caller supplied", {
  skip_fwd()
  d <- fwd_sim()
  args <- c("design", "ndups", "spacing", "block", "weights", "method")
  spy <- spy_on(limma::lmFit, args)
  expect_identical(lmFit_parallel(d$y, d$design, workers = 2L, chunks = 3L,
                                  parallel_backend = "serial", backend = spy$fn),
                   limma::lmFit(d$y, d$design))
  w <- stats::runif(ncol(d$y), 0.5, 2)
  expect_identical(lmFit_parallel(d$y, d$design, weights = w, workers = 2L, chunks = 3L,
                                  parallel_backend = "serial", backend = spy$fn),
                   limma::lmFit(d$y, d$design, weights = w))
  seen <- spy$seen()
  expect_length(seen, 2L)
  expect_identical(seen[[1L]], stats::setNames(args != "design", args))
  expect_identical(seen[[2L]], stats::setNames(!args %in% c("design", "weights"), args))
})

test_that("duplicateCorrelation_parallel follows a backend whose design default changed", {
  skip_fwd()
  d <- fwd_sim()
  drifted <- limma::duplicateCorrelation
  formals(drifted)$design <- quote(cbind(1, rep(0:1, length.out = ncol(object))))
  ref <- drifted(d$y, block = d$block)
  expect_false(identical(ref, limma::duplicateCorrelation(d$y, block = d$block)))
  expect_identical(
    duplicateCorrelation_parallel(d$y, block = d$block, workers = 2L, chunks = 3L,
                                  parallel_backend = "serial", backend = drifted),
    ref)
})

test_that("lmFit_parallel follows a backend whose design default changed", {
  skip_fwd()
  d <- fwd_sim()
  drifted <- limma::lmFit
  formals(drifted)$design <- quote(cbind(1, rep(0:1, length.out = ncol(object))))
  ref <- drifted(d$y)
  expect_false(identical(ref, limma::lmFit(d$y)))
  expect_identical(lmFit_parallel(d$y, workers = 2L, chunks = 3L, parallel_backend = "serial",
                                  backend = drifted),
                   ref)
})

test_that("lmFit_parallel refuses on the method and ndups the backend will actually use", {
  skip_fwd()
  d <- fwd_sim()
  robust <- limma::lmFit
  formals(robust)$method <- "robust"
  expect_error(lmFit_parallel(d$y, d$design, workers = 2L, backend = robust), "robust")
  expect_identical(lmFit_parallel(d$y, d$design, method = "ls", workers = 2L, chunks = 3L,
                                  parallel_backend = "serial", backend = robust),
                   limma::lmFit(d$y, d$design))
  dups <- limma::lmFit
  formals(dups)$ndups <- 2
  expect_error(lmFit_parallel(d$y, d$design, workers = 2L, backend = dups), "ndups")
})

test_that("a trim passed through lmFit_parallel reaches the consensus correlation", {
  skip_fwd()
  d <- fwd_sim()
  ref <- limma::lmFit(d$y, d$design, block = d$block, correlation = NULL, trim = 0.4)
  expect_false(identical(ref, limma::lmFit(d$y, d$design, block = d$block, correlation = NULL)))
  expect_identical(lmFit_parallel(d$y, d$design, block = d$block, correlation = NULL,
                                  trim = 0.4, workers = 2L, chunks = 3L,
                                  parallel_backend = "serial"),
                   ref)
})


test_that("removeBatchEffect_parallel hands the backend only the arguments the caller supplied", {
  skip_fwd()
  d <- fwd_sim()
  args <- c("batch", "batch2", "covariates", "design", "group")
  spy <- spy_on(limma::removeBatchEffect, args)
  expect_identical(removeBatchEffect_parallel(d$y, batch = d$batch, workers = 2L, chunks = 3L,
                                              parallel_backend = "serial", backend = spy$fn),
                   limma::removeBatchEffect(d$y, batch = d$batch))
  expect_identical(spy$seen(), list(stats::setNames(args != "batch", args)))
})

test_that("removeBatchEffect_parallel's lmFit forwards only what the original passed it", {
  skip_fwd()
  d <- fwd_sim()
  e <- new.env(parent = asNamespace("limma"))
  spy <- spy_on(limma::lmFit, c("ndups", "spacing", "block", "weights", "method"))
  e$lmFit <- spy$fn
  rbe <- limma::removeBatchEffect
  environment(rbe) <- e
  expect_identical(removeBatchEffect_parallel(d$y, batch = d$batch, workers = 2L, chunks = 3L,
                                              parallel_backend = "serial", backend = rbe),
                   limma::removeBatchEffect(d$y, batch = d$batch))
  expect_identical(spy$seen(), list(c(ndups = TRUE, spacing = TRUE, block = TRUE,
                                      weights = TRUE, method = TRUE)))
})

test_that("removeBatchEffect_parallel follows a backend whose covariates default changed", {
  skip_fwd()
  d <- fwd_sim()
  drifted <- limma::removeBatchEffect
  formals(drifted)$covariates <- quote(seq_len(ncol(x)))
  ref <- drifted(d$y)
  expect_false(identical(ref, limma::removeBatchEffect(d$y)))
  expect_identical(removeBatchEffect_parallel(d$y, workers = 2L, chunks = 3L,
                                              parallel_backend = "serial", backend = drifted),
                   ref)
})

test_that("removeBatchEffect_parallel matches limma on a vector with the design default", {
  skip_fwd()
  v <- stats::setNames(c(-0.63, 0.18, -0.84, 1.60, 0.33, -0.82), sprintf("g%d", 1:6))
  ref <- fwd_outcome(function() limma::removeBatchEffect(v, covariates = 3))
  expect_identical(fwd_outcome(function() removeBatchEffect_parallel(v, covariates = 3,
                                                                     workers = 1L)),
                   ref)
})


drifted_dupcor_backend <- function() {
  lns <- asNamespace("limma")
  dc <- get("duplicateCorrelation", envir = lns)
  b <- as.list(body(dc))
  n <- length(b)
  body(dc) <- as.call(c(b[seq_len(n - 1L)], list(quote(ngenes <- ngenes)), b[n]))
  e <- new.env(parent = lns)
  e$duplicateCorrelation <- dc
  lm2 <- limma::lmFit
  environment(lm2) <- e
  list(dupcor = dc, lmFit = lm2)
}

test_that("a duplicateCorrelation refusal carries the rnaparallel_refusal class", {
  skip_fwd()
  d <- fwd_sim()
  be <- drifted_dupcor_backend()
  expect_identical(be$dupcor(d$y, d$design, block = d$block),
                   limma::duplicateCorrelation(d$y, d$design, block = d$block))
  expect_error(duplicateCorrelation_parallel(d$y, d$design, block = d$block, workers = 2L,
                                             backend = be$dupcor),
               "ngenes", class = "rnaparallel_refusal")
  expect_error(duplicateCorrelation_parallel(d$y, d$design, workers = 2L),
               class = "rnaparallel_refusal")
})

test_that("lmFit_parallel runs the backend's own duplicateCorrelation whole when the split refuses", {
  skip_fwd()
  d <- fwd_sim()
  be <- drifted_dupcor_backend()
  ref <- limma::lmFit(d$y, d$design, block = d$block, correlation = NULL)
  expect_identical(lmFit_parallel(d$y, d$design, block = d$block, correlation = NULL,
                                  workers = 2L, chunks = 3L, parallel_backend = "serial",
                                  backend = be$lmFit),
                   ref)
})

test_that("an error that is not a refusal still surfaces from the correlation = NULL route", {
  skip_fwd()
  d <- fwd_sim()
  n <- 0L
  first_fails <- function(idx, f, workers) {
    n <<- n + 1L
    if (n == 1L) stop("boom on the first dispatch")
    lapply(idx, f)
  }
  withr::local_options(combat.progress = FALSE)
  expect_error(lmFit_parallel(d$y, d$design, block = d$block, correlation = NULL,
                              workers = 2L, chunks = 3L, parallel_backend = first_fails),
               "boom on the first dispatch")
  expect_identical(n, 1L)
})

test_that("lmFit_parallel with correlation = NULL caps workers once across the nested call", {
  skip_fwd()
  d <- fwd_sim()
  ref <- limma::lmFit(d$y, d$design, block = d$block, correlation = NULL)
  rss_calls <- 0L
  local_mocked_bindings(
    rp_mem_available = function() 250 / 0.8,
    rp_mem_rss = function() {
      rss_calls <<- rss_calls + 1L
      100 + 30 * (rss_calls - 1L)
    },
    .package = "rnaparallel")
  seen <- integer()
  spy <- function(idx, f, workers) { seen <<- c(seen, workers); lapply(idx, f) }
  withr::local_options(combat.mem.guard = TRUE, combat.mem.divergence = 1,
                       combat.progress = FALSE)
  got <- fwd_outcome(function() lmFit_parallel(d$y, d$design, block = d$block,
                                               correlation = NULL, workers = 3L, chunks = 3L,
                                               parallel_backend = spy))
  expect_identical(got$value, ref)
  expect_length(grep("workers need", got$warnings), 1L)
  expect_identical(rss_calls, 1L)
  expect_true(length(seen) >= 2L && all(seen == 2L))
})


test_that("rp_weights_fast agrees with asMatrixWeights on every weight shape", {
  skip_fwd()
  truth <- function(w, dim_full) {
    is.null(w) || !is.null(attr(limma::asMatrixWeights(w, dim_full), "arrayweights"))
  }
  set.seed(3)
  G <- 7L
  S <- 5L
  full_aw <- limma::asMatrixWeights(stats::runif(S), c(G, S))
  cases <- list(
    null = list(NULL, c(G, S)),
    scalar = list(2, c(G, S)),
    gene_vector = list(stats::runif(G), c(G, S)),
    array_vector = list(stats::runif(S), c(G, S)),
    one_row_matrix = list(matrix(stats::runif(S), 1L), c(G, S)),
    one_col_matrix = list(matrix(stats::runif(G), G, 1L), c(G, S)),
    one_by_one = list(matrix(2), c(G, S)),
    full_matrix = list(matrix(stats::runif(G * S), G, S), c(G, S)),
    full_arrayweights = list(full_aw, c(G, S)),
    full_data_frame = list(as.data.frame(matrix(stats::runif(G * S), G, S)), c(G, S)),
    square_vector = list(stats::runif(6L), c(6L, 6L)),
    square_one_row = list(matrix(stats::runif(6L), 1L), c(6L, 6L)),
    square_one_col = list(matrix(stats::runif(6L), 6L, 1L), c(6L, 6L)),
    one_array_vector = list(stats::runif(G), c(G, 1L)),
    one_gene_vector = list(stats::runif(S), c(1L, S)),
    one_gene_scalar = list(2, c(1L, S)))
  for (nm in names(cases)) {
    w <- cases[[nm]][[1L]]
    dim_full <- cases[[nm]][[2L]]
    expect_identical(rnaparallel:::rp_weights_fast(w, dim_full), truth(w, dim_full), info = nm)
  }
  expect_true(rnaparallel:::rp_weights_fast(full_aw, c(G, S)))
  expect_false(rnaparallel:::rp_weights_fast(stats::runif(6L), c(6L, 6L)))
})
