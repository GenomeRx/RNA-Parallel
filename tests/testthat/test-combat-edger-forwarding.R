## How ComBat_seq_parallel and calcNormFactors_parallel hand their arguments to the original.

fwd_in_process <- function(idx, f, workers) lapply(idx, f)

## The value or error message plus every warning and message, with printed output dropped.
fwd_conditions <- function(expr) {
  w <- character()
  m <- character()
  value <- withCallingHandlers(
    tryCatch({
      utils::capture.output(v <- expr)
      v
    }, error = function(e) structure(conditionMessage(e), class = "fwd_error")),
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") },
    message = function(c) { m <<- c(m, conditionMessage(c)); invokeRestart("muffleMessage") })
  list(value = value, warnings = w, messages = m)
}

## A wrapper `function(m, bb, x)` that passes its own `x` on as argument `arg` of `fn`.
fwd_wrapper <- function(fn, arg, extra = list()) {
  w <- function(m, bb, x) NULL
  body(w) <- as.call(c(list(fn, quote(m), quote(bb)), stats::setNames(list(quote(x)), arg),
                       extra))
  w
}

fwd_counts <- function(G = 200L, S = 8L, seed = 1L) {
  set.seed(seed)
  matrix(stats::rnbinom(G * S, mu = 40, size = 4), G, S)
}

test_that("an argument a wrapper passes through but never received fails as in sva", {
  skip_if_not_installed("sva")
  cts <- fwd_counts()
  bat <- rep(1:2, 4)
  for (a in c("group", "covar_mod", "full_mod", "shrink", "shrink.disp", "gene.subset.n")) {
    w_orig <- fwd_wrapper(quote(sva::ComBat_seq), a)
    w_comp <- fwd_wrapper(quote(ComBat_seq_parallel), a,
                          list(workers = 2L, parallel_backend = "serial"))
    ref <- fwd_conditions(w_orig(cts, bat))
    got <- fwd_conditions(w_comp(cts, bat))
    expect_identical(got$value, ref$value, info = a)
    expect_identical(got$warnings, ref$warnings, info = a)
  }
})

test_that("a method argument whose evaluation fails is evaluated once, as in edgeR", {
  skip_if_not_installed("edgeR")
  cts <- fwd_counts()
  for (obj in list(cts, edgeR::DGEList(cts))) {
    n <- 0L
    boom <- function() { n <<- n + 1L; stop("boom") }
    ref <- fwd_conditions(edgeR_norm(obj, method = boom()))
    n_ref <- n
    n <- 0L
    got <- fwd_conditions(calcNormFactors_parallel(obj, method = boom(), workers = 2L))
    expect_true(inherits(ref$value, "fwd_error"), info = class(obj)[1L])
    expect_identical(got, ref, info = class(obj)[1L])
    expect_identical(n, n_ref, info = class(obj)[1L])
  }
})

test_that("a method a wrapper passes through but never received fails as in edgeR", {
  skip_if_not_installed("edgeR")
  cts <- fwd_counts()
  w_orig <- function(m, me) edgeR_norm(m, method = me)
  w_comp <- function(m, me) calcNormFactors_parallel(m, method = me, workers = 2L)
  ref <- fwd_conditions(w_orig(cts))
  expect_true(inherits(ref$value, "fwd_error"))
  expect_identical(fwd_conditions(w_comp(cts)), ref)
})

test_that("the default label names a method passed in a variable", {
  skip_if_not_installed("edgeR")
  withr::local_options(list(combat.timing = TRUE, combat.timing.min = 0))
  cts <- fwd_counts()
  m <- "RLE"
  line <- paste(testthat::capture_messages(suppressWarnings(
    calcNormFactors_parallel(cts, method = m, workers = 2L))), collapse = "")
  expect_match(line, "calcNormFactors RLE", fixed = TRUE)
})

test_that("a common dispersion the across-batch shim cannot assemble is reported as stood down", {
  skip_if_not_installed("sva")
  txt <- paste(deparse(body(sva::ComBat_seq), width.cutoff = 500L), collapse = "\n")
  from <- c("return(estimateGLMCommonDisp(", ", subset = nrow(counts)))")
  to <- c("return(cbind(estimateGLMCommonDisp(", ", subset = nrow(counts)), 1))")
  for (k in seq_along(from)) {
    skip_if_not(grepl(from[k], txt, fixed = TRUE), paste("sva's ComBat_seq no longer contains", from[k]))
    txt <- gsub(from[k], to[k], txt, fixed = TRUE)
  }
  d <- sva::ComBat_seq
  body(d) <- parse(text = txt)[[1L]]
  environment(d) <- environment(sva::ComBat_seq)

  set.seed(4)
  cts <- matrix(stats::rnbinom(300 * 30, mu = 60, size = 4), nrow = 300)
  bat <- rep(1:3, each = 10L)
  ref <- quietly(d(cts, batch = bat, group = NULL))

  rnaparallel:::rp_count_reset()
  got <- quietly(ComBat_seq_parallel(cts, batch = bat, group = NULL, workers = 2L,
                                     parallel_backend = fwd_in_process, backend = d))
  expect_identical(got, ref)
  expect_true("estimateGLMCommonDisp across batches" %in% rnaparallel:::.rp_dispatch$fallback)
})
