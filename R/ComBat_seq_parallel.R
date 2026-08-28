## ComBat_seq_parallel.R
##
## Mirrors upstream ComBat_seq.R: the entry point and nothing else. Source this
## plus helper_seq_parallel.R to use the function without installing anything.
##
## The core function is still the original. This file does not contain a copy of
## ComBat-seq; it calls the original ComBat-seq function with six symbols rebound.

#' ComBat-seq with its hot paths parallelised
#'
#' Runs ComBat-seq itself. The algorithm is not reimplemented and not copied: the
#' original ComBat-seq function is called with six symbols rebound in a child of its own
#' environment, so `glmFit`, `glmFit.default`, `match_quantiles` and
#' `estimateGLMTagwiseDisp` resolve to row-parallel versions, `sapply` dispatches the
#' per-batch common-dispersion estimate across batches, and `lapply` dispatches the
#' per-batch tagwise estimate the same way, while every other symbol in the
#' body, `vec2mat`, `estimateGLMCommonDisp`, `getOffset`, `monte_carlo_int_NB`, still
#' resolves to the backend's own. Output is
#' `identical()` to what `sva::ComBat_seq` returns for the same input, not merely
#' close to it.
#'
#' @section Why a rebind and not a rewrite:
#' An earlier version transcribed the algorithm by hand so the loop could be
#' parallelised, and drifted from the original ComBat-seq twice in ways synthetic test data
#' could not expose. The gene filter was written `!all(x <= 1) && var(x) > 0`
#' against the original ComBat-seq's per-batch `all(x == 0)`, which returned 258 genes
#' uncorrected and shifted `lib.size` for every surviving gene as well. And the
#' copy dropped the model intercept the original ComBat-seq keeps, which corrupted every
#' dispersion estimate because `mod` also feeds `estimateGLMCommonDisp`. Both
#' passed on uniform Poisson counts and failed only on real sparse data. Rebinding
#' cannot drift, because there is no second copy of the algorithm to drift.
#'
#' @section Backends for ComBat-seq itself:
#' Works against Bioconductor `sva` or a sourced copy of upstream
#' \url{https://github.com/zhangyuqing/ComBat-seq}, whichever is loaded. Both
#' expose the same 8-argument signature and the same helpers; the helper is taken
#' from the backend's own environment, so the two can never be mixed.
#'
#' @section Backends for the parallelism:
#' `parallel_backend` selects the framework: `"mclapply"` (default, forks),
#' `"future"`, `"BiocParallel"`, `"foreach"`, or `"serial"`. All of them return
#' bit-identical results, because each row chunk is a pure function of its own
#' genes and every backend preserves chunk order. See [combat_backends()].
#'
#' With `"future"` the caller owns the plan; this package will not set one, since
#' that would stamp on a plan established for the whole session. Set e.g.
#' `future::plan(future::multicore, workers = 4)` first.
#'
#' `options(combat.fork = FALSE)` forces a serial run on any backend, with
#' identical output. That is the escape hatch if forking upsets an IDE. Windows
#' cannot fork and falls back to serial automatically, so results are correct there
#' and not faster.
#'
#' @section Dispatches too small to be worth a fork:
#' ComBat-seq calls the quantile match and the dispersion estimate once per batch, not
#' once per run, so batch count sets how thin each dispatch is. A 100-batch design hands
#' over a hundredth of the samples at a time, and on a small matrix that is not worth it:
#' measured 0.50x at 5,000 cells, 0.92x at 10,000 and 1.39x at 20,000.
#' Dispatches below `getOption("combat.min.cells", 20000)` therefore run serially.
#' Without that gate 500 genes by 1,000 samples across 100 batches ran at 0.80x,
#' slower than not parallelising at all. Set the option to 0 to restore the old
#' unconditional behaviour; output is identical either way.
#'
#' The threshold is per path, because the paths do not cost the same per cell.
#' `qnbinom` in the quantile match earns a fork at 20,000 cells; dispersion
#' estimation is cheaper and does not break even below 30,000, measuring
#' 0.65x at 10,000 and 0.79x at 20,000, so it has its own
#' `getOption("combat.min.disp.cells", 30000)`, and the GLM fit its own
#' `getOption("combat.min.glm.cells", 1e5)`. One shared threshold made the
#' dispersion split a net loss on small matrices with many batches.
#'
#' @param counts Raw count matrix, genes in rows, samples in columns.
#' @param batch Batch vector, one entry per column of `counts`. Call
#'   `droplevels()` on a factor first: ComBat-seq counts samples per level, so a
#'   level left over from subsetting reads as a batch with no samples and it stops
#'   with a message that does not mention levels.
#' @param group Optional biological condition, preserved by the correction.
#' @param covar_mod Optional model matrix of other covariates to preserve, e.g.
#'   `cbind(cov1, cov2)`.
#' @param full_mod Fit the full model including `group`.
#' @param shrink Shrink dispersion estimates towards a common value.
#' @param shrink.disp Shrink the dispersion as well as the mean.
#' @param gene.subset.n Genes sampled for shrinkage. Seed before comparing two
#'   runs with this set: the original ComBat-seq's `monte_carlo_int_NB` calls `sample()`, so
#'   unseeded back-to-back runs draw different subsets and differ legitimately.
#' @param workers Maximum concurrent worker processes. `NULL`, the default, resolves to
#'   `min(8, detectCores() - 2)`: 8 is where the measured curves flatten, and the two
#'   held back leave the master and the rest of the machine a core each. Pass a smaller
#'   number when a second R process is also forking; six workers alongside one
#'   kernel-panicked a 24 GB machine, and no core count read here can see that session.
#'
#'   Going past your PERFORMANCE core count can be slower than stopping at it, and the loss is
#'   worst on the largest inputs -- which is exactly where people reach for more workers.
#'   Measured on an 8-core machine with 4 performance cores, TMM on 15,000 genes by 9,000
#'   specimens: 4 workers 10.74 s, 6 workers 7.26 s, 8 workers 11.60 s. Eight was slower than
#'   four. The default resolved to 6 there and was optimal; raising it by hand made it worse.
#' @param chunks Row chunks per parallel step. Defaults to `workers`, which is almost always
#'   what you want -- passing `chunks = workers` explicitly is redundant. Raise it above
#'   `workers` only to cut peak memory per worker, at slightly more fork overhead.
#' @param parallel_backend One of [combat_backends()], or a function
#'   `function(idx, f, workers)` returning a list in the order of `idx`. The
#'   function form plugs in any framework this package does not name, including
#'   `parallel::parLapply` over a FORK or PSOCK cluster, `pbapply::pblapply`,
#'   `furrr::future_map`, `doFuture::%dofuture%`, or any `BiocParallel` BPPARAM.
#'   Defaults to `getOption("combat.backend", "mclapply")`.
#' @param backend Optional ComBat-seq function to wrap. Defaults to
#'   `sva::ComBat_seq`, falling back to a visible top-level `ComBat_seq`.
#'
#' @param label Optional name for this call in the timing line, when
#'   `options(combat.timing = TRUE)` is set. Defaults to the companion and the matrix shape,
#'   e.g. `ComBat-seq 18,270 x 1,500`; pass a cohort name to tell calls apart in a loop.
#' @return A gene-by-sample matrix of adjusted counts, `identical()` to what the
#'   backend's own `ComBat_seq` returns for the same input. `workers = 1` takes a
#'   plain `lapply` and adds no parallelism, so it is the backend function itself
#'   with a thin dispatch wrapper.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("sva", quietly = TRUE)) {
#'   set.seed(1)
#'   count_matrix <- matrix(rnbinom(1600, mu = 50, size = 5), nrow = 200)
#'   batch <- rep(1:2, each = 4)
#'
#'   # upstream's covariate example, with parallelism added
#'   cov1 <- rep(c(0, 1), 4)
#'   cov2 <- c(0, 0, 1, 1, 0, 0, 1, 1)
#'   covar_mat <- cbind(cov1, cov2)
#'   adjusted_counts <- ComBat_seq_parallel(count_matrix, batch = batch, group = NULL,
#'                                          covar_mod = covar_mat, workers = 4L)
#'
#'   identical(adjusted_counts,
#'             sva::ComBat_seq(count_matrix, batch = batch, group = NULL,
#'                             covar_mod = covar_mat))
#' }
#' }
#' @export
ComBat_seq_parallel <- function(counts, batch, group = NULL, covar_mod = NULL,
                                full_mod = TRUE, shrink = FALSE, shrink.disp = FALSE,
                                gene.subset.n = NULL, workers = NULL, chunks = NULL,
                                parallel_backend = getOption("combat.backend", "mclapply"),
                                backend = NULL, label = NULL) {
  # no run leaves workers behind, crashed or not; children that predate this call are
  # someone else's and are spared
  .spare <- combat_children()
  on.exit(combat_reap(.spare), add = TRUE)

  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("edgeR is required: BiocManager::install(\"edgeR\")", call. = FALSE)
  }
  if (is.numeric(workers) && length(workers) == 1L && is.finite(workers) &&
      workers != trunc(workers)) {
    stop("`workers` must be a whole number, not ", workers, call. = FALSE)
  }
  workers <- rp_prologue(workers)

  # timing and quieting are on.exit hooks, so an error unwinds the sink and still reports the
  # elapsed line: a failed run says where it failed instead of vanishing. Placed after the
  # prologue because that is what resolves `workers` from NULL to a number worth printing.
  .rp <- rp_step_begin(label, "ComBat-seq", counts, parallel_backend, workers)
  on.exit(rp_step_end(.rp), add = TRUE)
  if (!is.function(parallel_backend)) {
    parallel_backend <- match.arg(parallel_backend, combat_backends())
  }

  be <- combat_backend(backend)
  mq <- be$match_quantiles

  env <- new.env(parent = be$env)

  # The DGEList branch. offset comes from the WHOLE object before any splitting,
  # which is what edgeR:::glmFit.DGEList does and the one quantity a per-chunk fit
  # could not reconstruct from its own slice.
  env$glmFit <- function(y, design = NULL, dispersion = NULL, prior.count = 0.125,
                         start = NULL, ...) {
    if (inherits(y, "DGEList")) {
      glmFit_rows_parallel(y$counts, design = design, dispersion = dispersion,
                           offset = edgeR::getOffset(y), weights = y$weights,
                           prior.count = prior.count, start = start,
                           workers = workers, chunks = chunks,
                           parallel_backend = parallel_backend)
    } else {
      # An offset in the dots is exactly what makes a matrix call splittable, so honour it.
      # Without one each worker rebuilds lib.size from its own slice of genes, which moved
      # coefficients by 1.7 in testing and raises nothing, so that case is refused rather
      # than guessed. Unreachable from sva::ComBat_seq, which only calls glmFit on a DGEList.
      dots <- list(...)
      if (!is.null(dots$offset)) {
        return(glmFit_rows_parallel(y, design = design, dispersion = dispersion,
                                    offset = dots$offset, weights = dots$weights,
                                    prior.count = prior.count, start = start,
                                    workers = workers, chunks = chunks,
                                    parallel_backend = parallel_backend))
      }
      stop("glmFit was called on a matrix with no offset. rnaparallel cannot split ",
           "that by row: a worker would rebuild library sizes from its own genes and every ",
           "fitted value would change silently. Pass an explicit offset.", call. = FALSE)
    }
  }

  env$glmFit.default <- function(y, design = NULL, dispersion = NULL, offset = NULL,
                                 lib.size = NULL, weights = NULL,
                                 prior.count = 0.125, start = NULL, ...) {
    # accepted by the signature, so refuse it rather than discard it: a row split cannot
    # honour a whole-matrix lib.size, and silently dropping it changes every fitted value
    if (!is.null(lib.size)) {
      stop("pass offset, not lib.size: rnaparallel splits by row and cannot honour ",
           "a whole-matrix lib.size", call. = FALSE)
    }
    glmFit_rows_parallel(y, design = design, dispersion = dispersion, offset = offset,
                         weights = weights, prior.count = prior.count, start = start,
                         workers = workers, chunks = chunks,
                         parallel_backend = parallel_backend)
  }

  env$match_quantiles <- function(counts_sub, old_mu, old_phi, new_mu, new_phi) {
    match_quantiles_parallel(mq, counts_sub, old_mu, old_phi, new_mu, new_phi,
                             workers = workers, chunks = chunks,
                             parallel_backend = parallel_backend)
  }

  # Dispersion is estimated once per batch and used to be the whole serial block.
  # Only exact because ComBat-seq passes prior.df = 0; the helper checks that itself
  # and hands anything else to edgeR unsplit.
  env$estimateGLMTagwiseDisp <- function(y, design = NULL, offset = NULL, dispersion = NULL,
                                         prior.df = 10, trend = TRUE, span = NULL,
                                         AveLogCPM = NULL, weights = NULL, ...) {
    estimateGLMTagwiseDisp_rows_parallel(y, design = design, dispersion = dispersion,
                                         offset = offset, prior.df = prior.df, trend = trend,
                                         span = span, AveLogCPM = AveLogCPM, weights = weights,
                                         workers = workers, chunks = chunks,
                                         parallel_backend = parallel_backend)
  }

  # Common dispersion is estimated once per batch and was the last block still running
  # serially. Profiling a 10,000 by 500 run put it at 13.1% of total time, which by Amdahl
  # caps the achievable speedup near 7x however many workers are added.
  #
  # It cannot be split by row: the objective is optimize() over a SUM of per-gene adjusted
  # profile likelihoods, and splitting that sum changes floating-point accumulation order,
  # which can move the argmax in the last ulp. The batches, however, are independent of each
  # other, contain no RNG, and each returns one scalar. Dispatching across batches is therefore
  # exact where dispatching across rows would not be.
  #
  # The backend body contains exactly two sapply calls: `sapply(batches_ind, length)`, whose
  # FUN is the primitive `length`, and the disp_common call, whose FUN is a closure mentioning
  # estimateGLMCommonDisp. The predicate below selects the second and only the second, and
  # anything it does not recognise falls through to base::sapply untouched. It never sees the
  # lapply that drives monte_carlo_int_NB, whose sample() draws depend on the RNG state the
  # previous batch left behind and which therefore must stay sequential.
  env$sapply <- function(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
    targeted <- tryCatch(
      is.function(FUN) && !is.primitive(FUN) &&
        "estimateGLMCommonDisp" %in% all.names(body(FUN)),
      error = function(e) FALSE)

    # installed under a base-R name, so anything this shim does not reproduce exactly goes
    # back to base: it neither honours simplify = FALSE nor sets names on the result
    if (!targeted || length(X) < 2L || !isTRUE(simplify) ||
        !is.null(names(X)) || is.character(X)) {
      return(base::sapply(X, FUN, ..., simplify = simplify, USE.NAMES = USE.NAMES))
    }

    # No size gate here: each element is a whole common-dispersion estimate over every gene,
    # which is always worth dispatching. Inf says that, where a real cell count would only be
    # computed and then ignored by min_cells = 0.
    parts <- combat_parallel_check(
      combat_parallel_lapply(as.list(X), function(i) FUN(i, ...), workers,
                             parallel_backend, cells = Inf, min_cells = 0,
                             # one scalar per batch, and batch counts exceed worker counts on
                             # real designs, so one fork per worker beats one fork per batch.
                             # Measured at 100 batches: 1428 ms to 1063 ms.
                             preschedule = TRUE),
      "estimateGLMCommonDisp across batches", as.list(X))

    out <- unlist(parts, use.names = FALSE)
    # a batch returning anything but one number means the assumption above is wrong for this
    # input, so hand the whole call back to R rather than assemble something unverified
    if (length(out) != length(X) || !is.numeric(out)) {
      return(base::sapply(X, FUN, ..., simplify = simplify, USE.NAMES = USE.NAMES))
    }
    out
  }

  # The per-batch tagwise dispersion, dispatched ACROSS BATCHES rather than across gene rows.
  #
  # This is the exact route to the work a row split cannot touch. ComBat-seq computes the
  # tagwise dispersion inside `lapply(1:n_batch, ...)`, and each element is a complete
  # estimate on that batch's own columns: independent of every other batch, no RNG, no state
  # carried between them. Splitting gene ROWS inside one batch is what reaches edgeR's
  # one-group kernel and is refused; splitting BATCHES does not go near it.
  #
  # It matters because that stage dominates. Profiled on 6,000 genes by 400 samples across 10
  # batches with group = NULL, estimateGLMTagwiseDisp was 60% of the whole corrected run while
  # match_quantiles was 7%.
  #
  # The predicate selects that one lapply and nothing else. ComBat-seq calls lapply seven
  # times, and one of them drives monte_carlo_int_NB, whose sample() draws depend on the RNG
  # state the previous batch left behind: dispatching it would change results. It is excluded
  # by name as well as by the absence of the tagwise call, because getting this wrong is
  # silent.
  env$lapply <- function(X, FUN, ...) {
    targeted <- tryCatch({
      nm <- all.names(body(FUN))
      is.function(FUN) && !is.primitive(FUN) &&
        "estimateGLMTagwiseDisp" %in% nm &&
        !any(c("mcint_fun", "monte_carlo_int_NB", "sample", "rnorm", "runif") %in% nm)
    }, error = function(e) FALSE)

    # installed under a base-R name, so anything this shim does not reproduce exactly goes
    # straight back to base
    if (!targeted || length(X) < 2L || !is.null(names(X))) {
      return(base::lapply(X, FUN, ...))
    }

    # Each element is one whole-matrix estimate over every gene, so it always earns a
    # dispatch. Nesting is prevented by combat_parallel_lapply, which marks the worker
    # process it dispatches into and runs serially when it finds itself already inside one,
    # so a rebound tagwise in the worker does not open a second pool. That guard used to be
    # `mc.allow.recursive = FALSE`, which is an mclapply argument and covered the fork branch
    # only; on Windows, where foreach/PSOCK is the sole working backend, this nested to
    # workers + workers^2 processes until the flag replaced it.
    # idx is deliberately not passed: its row check compares a chunk's returned rows against
    # the indices it was given, and here one index returns a dispersion per gene. Dead workers
    # and thrown errors are still caught, and the shape is checked below instead.
    parts <- combat_parallel_check(
      combat_parallel_lapply(as.list(X), function(i) FUN(i, ...), workers,
                             parallel_backend, cells = Inf, min_cells = 0),
      "estimateGLMTagwiseDisp across batches")

    # Compared against the matrix FUN actually operates on, not the entry-point argument:
    # ComBat-seq drops genes that are all zero within a batch before this lapply, so on any
    # sparse input the unfiltered row count failed every batch and the whole parallel stage
    # was silently recomputed serially. environment(FUN) is the vendor's own frame.
    expected <- tryCatch(nrow(get("counts", envir = environment(FUN), inherits = FALSE)),
                         error = function(e) NA_integer_)
    ok_shape <- !is.na(expected) && length(parts) == length(X) &&
      all(vapply(parts, function(z) is.numeric(z) && length(z) == expected, logical(1)))
    if (!ok_shape) return(base::lapply(X, FUN, ...))
    parts
  }

  f <- be$fn
  environment(f) <- env
  f(counts = counts, batch = batch, group = group, covar_mod = covar_mod,
    full_mod = full_mod, shrink = shrink, shrink.disp = shrink.disp,
    gene.subset.n = gene.subset.n)
}
