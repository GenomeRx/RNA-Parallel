## limma_rbe_parallel.R
##
## removeBatchEffect is 26 lines and all but two of them are argument shaping. Its entire cost
## is one `lmFit` call, written as a bare symbol, followed by `x - beta %*% t(X.batch)`. So the
## companion is not a reimplementation and not even a split: it rebinds `lmFit` to
## `lmFit_parallel` in a child of limma's own environment and calls the original unchanged.
## Nothing in the body reduces across genes, which is why a row split inside lmFit is exact
## here for the same reason it is exact in lmFit_parallel.

#' limma's `removeBatchEffect` with its `lmFit` parallelized
#'
#' Runs `limma::removeBatchEffect` itself. The one expensive call inside it, `lmFit`, is rebound
#' to [lmFit_parallel()] in a child of limma's environment; every other symbol in the body still
#' resolves to limma's own code, and the returned matrix is `identical()` to the original's.
#'
#' @section Why this one is worth having:
#' In a batch-effect PCA it runs once per cohort and once more on the pooled matrix, so a
#' five-cohort analysis calls it six times.
#'
#' @section What is not parallelized:
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
#' @param x,batch,batch2,covariates,design,group Passed to `limma::removeBatchEffect` unchanged.
#'   One left out is not passed at all, so it takes the backend's own default.
#' @param ... Passed through to the underlying `lmFit`. `method = "robust"` and `ndups >= 2`
#'   are refused by [lmFit_parallel()] for the reasons given there, so they are refused here.
#' @param backend Optional `removeBatchEffect` to wrap. Its `lmFit` is the one its own
#'   environment resolves, so the two always come from the same copy of limma. Defaults to
#'   `limma::removeBatchEffect`.
#' @param label Optional name for this call. It is the stage name on the progress bar, which
#'   draws by default on macOS and Linux in an interactive session or on a terminal
#'   (`options(combat.progress = FALSE)` turns it off), and in the timing line when
#'   `options(combat.timing = TRUE)` is set. Defaults to the companion and the matrix shape,
#'   e.g. `removeBatchEffect 22,000 x 948`; pass a cohort name to tell calls apart in a loop.
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
                                       design = NULL, group = NULL, ...,
                                       workers = NULL, chunks = NULL,
                                       parallel_backend = getOption("combat.backend", combat_default_backend()),
                                       backend = NULL, label = NULL) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("limma is required: BiocManager::install(\"limma\")", call. = FALSE)
  }
  .spare <- combat_children()
  on.exit(combat_reap(.spare), add = TRUE)

# Validated here but capped once, below, where memory is read after removeBatchEffect's own copies.
  workers <- rp_uncapped(rp_prologue(workers))

  # timing and quieting are on.exit hooks, so an error unwinds the sink and still reports the
  # elapsed line: a failed run says where it failed instead of vanishing. Placed after the
  # prologue because that is what resolves `workers` from NULL to a number worth printing.
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
  renv$lmFit <- function(...) {
    capped <- rp_mem_cap(workers)
# The lmFit the backend's own body would reach, so the pair never comes from two copies of limma.
    lm_fn <- get("lmFit", envir = env, mode = "function", inherits = TRUE)
    rp_uncapped(lmFit_parallel(..., workers = capped, chunks = chunks,
                               parallel_backend = parallel_backend, backend = lm_fn))
  }
  environment(fn) <- renv
# Only the arguments the caller supplied are forwarded, so every default is the backend's own.
  given <- intersect(c("batch", "batch2", "covariates", "design", "group"), names(match.call()))
  fwd <- lapply(given, as.name)
  names(fwd) <- given
  eval(as.call(c(list(quote(fn), x = quote(x)), fwd, quote(...))))
}
