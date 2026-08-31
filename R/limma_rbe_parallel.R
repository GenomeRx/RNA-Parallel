## limma_rbe_parallel.R
##
## removeBatchEffect is 26 lines and all but two of them are argument shaping. Its entire cost
## is one `lmFit` call, written as a bare symbol, followed by `x - beta %*% t(X.batch)`. So the
## companion is not a reimplementation and not even a split: it rebinds `lmFit` to
## `lmFit_parallel` in a child of limma's own environment and calls the original unchanged.
## Nothing in the body reduces across genes, which is why a row split inside lmFit is exact
## here for the same reason it is exact in lmFit_parallel.

#' limma's `removeBatchEffect` with its `lmFit` parallelised
#'
#' Runs `limma::removeBatchEffect` itself. The one expensive call inside it, `lmFit`, is rebound
#' to [lmFit_parallel()] in a child of limma's environment; every other symbol in the body still
#' resolves to limma's own code, and the returned matrix is `identical()` to the original's.
#'
#' @section Why this one is worth having:
#' It is called once per cohort in a batch-effect PCA and again on the pooled matrix, and its
#' cost grows worse than linearly in samples. Measured against limma 3.62.2 on 22,000 genes:
#' 8.5 s at 948 samples and 154 batches, 44.2 s at 3,000 samples and 250 batches, and it did not
#' finish inside ten minutes at 9,493 samples. A pipeline that calls it six times a run pays
#' that six times.
#'
#' @section What is not parallelised:
#' The final `x - beta %*% t(X.batch)` is one BLAS call over the whole matrix. It is not split,
#' because a matrix product is not row-associative in floating point once BLAS is threaded, and
#' the whole product costs a fraction of the fit it follows.
#'
#' @section When this is worth reaching for:
#' On large matrices, and not before. This companion IS the `lmFit` call inside the original, so
#' it inherits that function's answer exactly: below the least-squares gate there is nothing to
#' split and the companion is the original plus about a millisecond. Measured on an M3 at the
#' default worker count, companion against original, every arm `identical()`: 1.00x at 2,000
#' genes by 20 samples, 0.92x at 20,000 x 50, 0.98x at 20,000 x 200, then 1.91x at 20,000 x 500
#' once the matrix is large enough for the split to run at all.
#'
#' A wide batch design is what makes it worth having: the original's cost grows worse than
#' linearly in samples, measured 8.5 s at 948 samples and 154 batches against 44.2 s at 3,000
#' samples and 250 batches.
#'
#' @param x,batch,batch2,covariates,design,group Passed to `limma::removeBatchEffect` unchanged.
#' @param ... Passed through to the underlying `lmFit`. `method = "robust"` and `ndups >= 2`
#'   are refused by [lmFit_parallel()] for the reasons given there, so they are refused here.
#' @inheritParams lmFit_parallel
#' @return The batch-corrected matrix `limma::removeBatchEffect` returns, `identical()` to it.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' y <- matrix(rnorm(20000), nrow = 1000)
#' b <- rep(1:4, length.out = ncol(y))
#' identical(removeBatchEffect_parallel(y, batch = b, workers = 4L),
#'           limma::removeBatchEffect(y, batch = b))
#' }
#' @export
removeBatchEffect_parallel <- function(x, batch = NULL, batch2 = NULL, covariates = NULL,
                                       design = matrix(1, ncol(x), 1), group = NULL, ...,
                                       workers = NULL, chunks = NULL,
                                       parallel_backend = getOption("combat.backend", combat_default_backend()),
                                       backend = NULL, label = NULL) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("limma is required: BiocManager::install(\"limma\")", call. = FALSE)
  }
  .spare <- combat_children()
  on.exit(combat_reap(.spare), add = TRUE)

  workers <- rp_prologue(workers)

  .rp <- rp_step_begin(label, "removeBatchEffect", x, parallel_backend, workers)
  on.exit(rp_step_end(.rp), add = TRUE)

  fn <- if (is.null(backend)) limma::removeBatchEffect else backend
  if (!is.function(fn)) stop("`backend` must be a function", call. = FALSE)
  env <- environment(fn)
  if (is.null(env)) {
    stop("the limma backend has no environment, so lmFit cannot be reached.", call. = FALSE)
  }

  # The same gate every other rebind here carries. A limma that wrote `limma::lmFit(...)` would
  # turn this into a pass-through: correct output, original speed, and nothing to see from outside.
  if (!("lmFit" %in% rp_bare_call_heads(body(fn)))) {
    stop("this limma removeBatchEffect no longer calls lmFit as a bare symbol, so rebinding ",
         "cannot reach it. The companion would run the original serially while still returning ",
         "identical() output, which no equivalence test can detect. Refusing to run.",
         call. = FALSE)
  }

  renv <- new.env(parent = env)
  renv$lmFit <- function(object, design = NULL, ndups = NULL, spacing = NULL, ...) {
    lmFit_parallel(object, design = design, ndups = ndups, spacing = spacing, ...,
                   workers = workers, chunks = chunks,
                   parallel_backend = parallel_backend)
  }
  environment(fn) <- renv
  fn(x = x, batch = batch, batch2 = batch2, covariates = covariates,
     design = design, group = group, ...)
}
