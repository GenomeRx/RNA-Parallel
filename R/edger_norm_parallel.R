## edger_norm_parallel.R
##
## The edgeR normalisation companion. Mirrors the role helper_seq_parallel.R plays for
## ComBat-seq, and reuses its dispatch layer unchanged: combat_row_chunks,
## combat_parallel_lapply, combat_parallel_check, combat_row_order.
##
## Nothing here reimplements edgeR. `calcNormFactors.default` is called ONCE, unchanged, on
## the whole matrix, with four internals rebound in a child of edgeR's own namespace. Two
## separate wins come out of that. Inside `.calcFactorTMM`, `rank` is computed once per
## vector instead of twice. Inside the vendor's per-column loops, the columns are dispatched
## across workers.


# ---- why the vendor is called once, not once per block -----------------------
#
# The parallel axis is SAMPLE COLUMNS, never gene rows: every method here ranks, medians or
# takes a quantile across genes, so a row block is not a smaller version of the problem.
#
# Five quantities are pooled across columns and would change if the vendor were called on a
# column block instead. Calling the vendor once on the full matrix is what hoists all five,
# because edgeR itself computes them, on every column, exactly once:
#
#   1. lib.size <- colSums(x)
#   2. allzero <- .rowSums(x > 0, ...) == 0L, and the row drop that follows it
#   3. refColumn, from .calcFactorQuantile's f75 for TMM or colSums(sqrt(x)) for TMMwsp
#   4. gm <- exp(rowMeans(log(data))) in .calcFactorRLE, a cross-SAMPLE reduction sitting
#      exactly where the column axis cuts. Rebuilt per block it moved RLE by 0.043.
#   5. the tail centering f <- f/exp(mean(log(f)))
#
# Number 1 is also why blocks are the wrong shape for this function at all. colSums over
# FRACTIONAL doubles is not block-associative on arm64, where sizeof.longdouble is 8:
# summing a matrix in pieces and adding the pieces differs from summing it whole by about
# 1.8e-14 relative. Integer counts are exact either way, but a caller normalising
# already-scaled data is not, and identical() does not grade on a curve.
#
# What is left after those five is per column and nothing else, so that is what gets split.


# ---- backend resolution ------------------------------------------------------

#' Resolve an edgeR normalisation backend and the internals it calls
#'
#' Everything is taken from one environment, the backend's own, so the vendor function and
#' the four internals it dispatches to can never be drawn from two different copies of
#' edgeR. `rank`, `quantile` and `apply` are resolved through that same environment rather
#' than named directly, so a rebind always shadows exactly the function the vendor would
#' otherwise have reached.
#'
#' @param fn A `calcNormFactors` default method. Defaults to edgeR's.
#' @return A list with `fn`, `env`, the four `.calcFactor*` internals, the three base
#'   functions the shims stand in for, and the name of each batched loop's index.
#' @noRd
calcnorm_backend <- function(fn = NULL) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("edgeR is required: BiocManager::install(\"edgeR\")", call. = FALSE)
  }
  ns <- asNamespace("edgeR")
  if (is.null(fn)) fn <- get("calcNormFactors.default", envir = ns, inherits = FALSE)
  if (!is.function(fn)) stop("`backend` must be a function", call. = FALSE)

  env <- environment(fn)
  if (is.null(env)) {
    stop("the edgeR backend has no environment, so its internals cannot be resolved.",
         call. = FALSE)
  }

  # gate: the backend must take every argument we forward, or a silent upstream signature
  # change would surface as a wrong normalisation factor rather than an error
  need <- c("object", "lib.size", "method", "refColumn", "logratioTrim", "sumTrim",
            "doWeighting", "Acutoff", "p")
  missing_args <- setdiff(need, names(formals(fn)))
  if (length(missing_args)) {
    stop("this calcNormFactors backend is missing argument(s): ",
         paste(missing_args, collapse = ", "),
         ". rnaparallel was written against the signature edgeR 4.4.2 exposes.",
         call. = FALSE)
  }

  internals <- c(".calcFactorTMM", ".calcFactorTMMwsp", ".calcFactorRLE",
                 ".calcFactorQuantile")
  absent <- internals[!vapply(internals, exists, logical(1), envir = env, inherits = TRUE)]
  if (length(absent)) {
    stop("could not find ", paste(absent, collapse = ", "), " alongside the ",
         "calcNormFactors backend, so the per-column work cannot be reached.", call. = FALSE)
  }
  got <- lapply(internals, get, envir = env, inherits = TRUE)
  names(got) <- internals

  # the two TMM internals are called by name with eight arguments the shim has to reproduce
  # exactly, so a rename or a reorder there must not be absorbed by the shim's own defaults
  tmm_need <- c("obs", "ref", "libsize.obs", "libsize.ref", "logratioTrim", "sumTrim",
                "doWeighting", "Acutoff")
  for (nm in c(".calcFactorTMM", ".calcFactorTMMwsp")) {
    if (!identical(names(formals(got[[nm]])), tmm_need)) {
      stop("`", nm, "` has signature (", paste(names(formals(got[[nm]])), collapse = ", "),
           ") but rnaparallel dispatches it on (", paste(tmm_need, collapse = ", "),
           "). Refusing to run rather than pass arguments into a changed function.",
           call. = FALSE)
    }
  }

  # A rebind only reaches a call written as a bare symbol. If any of these is ever
  # namespace-qualified, the companion degrades to a pass-through: output stays identical(),
  # every equivalence test still passes, and nothing is parallelised or hoisted. Fail loudly.
  # all.names() would flatten `edgeR:::.calcFactorTMM` into its parts and read as reachable,
  # so collect the heads of calls whose head is a bare symbol instead.
  reach <- list(fn = internals, .calcFactorTMM = "rank", .calcFactorRLE = "apply",
                .calcFactorQuantile = "quantile")
  for (nm in names(reach)) {
    host <- if (nm == "fn") fn else got[[nm]]
    unreachable <- setdiff(reach[[nm]], rp_bare_call_heads(body(host)))
    if (length(unreachable)) {
      stop("this edgeR backend no longer calls ", paste(unreachable, collapse = ", "),
           " as a bare symbol inside ", if (nm == "fn") "calcNormFactors.default" else nm,
           ", so rebinding cannot reach it. The package would run serially while still ",
           "returning identical() output, which no equivalence test can detect. ",
           "Refusing to run.", call. = FALSE)
    }
  }

  # The column loops are what the shims batch, and each shim answers by the loop's own index
  # rather than by counting its calls. So the index symbol is read out of the body here
  # instead of being hardcoded, and a body carrying more than one of them is refused: two
  # loops with different index names would leave a shim reading a variable that is not the
  # counter, and every factor would come back off by whatever the two indices disagree on.
  for_index <- function(e, acc = character()) {
    if (!is.call(e)) return(acc)
    if (identical(e[[1L]], as.name("for"))) acc <- c(acc, as.character(e[[2L]]))
    for (i in seq_along(e)) if (is.call(e[[i]])) acc <- for_index(e[[i]], acc)
    unique(acc)
  }
  loops <- list(fn = for_index(body(fn)),
                .calcFactorQuantile = for_index(body(got[[".calcFactorQuantile"]])))
  for (nm in names(loops)) {
    if (length(loops[[nm]]) != 1L) {
      stop("expected exactly one for-loop index in ",
           if (nm == "fn") "calcNormFactors.default" else nm, ", found ",
           if (length(loops[[nm]])) paste(loops[[nm]], collapse = ", ") else "none",
           ". rnaparallel batches that loop and answers by its index, so it refuses to ",
           "guess which variable the index is.", call. = FALSE)
    }
  }

  list(fn = fn, env = env,
       tmm = got[[".calcFactorTMM"]], wsp = got[[".calcFactorTMMwsp"]],
       rle = got[[".calcFactorRLE"]], quant = got[[".calcFactorQuantile"]],
       tmm_loop = loops$fn, quant_loop = loops$.calcFactorQuantile,
       rank = get("rank", envir = env, inherits = TRUE),
       quantile = get("quantile", envir = env, inherits = TRUE),
       apply = get("apply", envir = env, inherits = TRUE))
}


# ---- the rank hoist ----------------------------------------------------------

#' rank(), computed once per vector instead of twice
#'
#' `.calcFactorTMM` builds its trim mask as
#' `(rank(logR) >= loL & rank(logR) <= hiL) & (rank(absE) >= loS & rank(absE) <= hiS)`.
#' R has no common subexpression elimination, so that is four sorts where two would do, and
#' `base::rank` is 40.8% of a plain limma-voom pipeline. The repeat is immediate, so a
#' one-slot cache catches it, and `identical()` on the vector is a memcmp against a sort.
#'
#' Bit-identical by construction, since a cache hit returns the object the first call
#' produced. Measured on 20000 x 200: 1.28s to 0.81s, identical() TRUE, max abs diff 0. No
#' workers involved, so `workers = 1` gets this too.
#'
#' @param rank0 The `rank` the vendor would otherwise have reached.
#' @return A closure standing in for `rank`.
#' @noRd
rp_rank_once <- function(rank0) {
  last <- NULL
  ranked <- NULL
  function(x, na.last = TRUE,
           ties.method = c("average", "first", "last", "random", "max", "min")) {
    # anything but the default call goes straight through: ties.method = "random" is not
    # even a function of its input, so caching it would be wrong as well as pointless
    if (!missing(na.last) || !missing(ties.method)) return(rank0(x, na.last, ties.method))
    if (!is.null(last) && identical(last, x)) return(ranked)
    ranked <<- rank0(x)
    last <<- x
    ranked
  }
}


# ---- the column split --------------------------------------------------------

#' Dispatch a per-column function across column chunks
#'
#' `combat_row_chunks` is a pure index splitter, so the column count goes in where a row
#' count normally would. The axis is COLUMNS here, not rows.
#'
#' @param ncols Number of columns to cover.
#' @param cells Size of the matrix being worked on, for the fork size gate.
#' @param f Function of one index vector, returning one value per column in it.
#' @param min_cells Fork threshold for this path.
#' @param what Label used in error messages.
#' @return One value per column, in column order.
#' @noRd
rp_norm_cols <- function(ncols, cells, f, workers, chunks, parallel_backend, min_cells,
                         what) {
  # the two-row floor is an lm.fit row-axis hazard; columns have no analogue
  idx <- combat_row_chunks(ncols, workers = workers, chunks = chunks)
  parts <- combat_parallel_check(
    combat_parallel_lapply(idx, f, workers, parallel_backend, cells = cells,
                           min_cells = min_cells),
    what, idx)
  out <- unlist(parts)                     # names kept: apply() names RLE from colnames
  if (length(out) != ncols) {
    stop(what, ": bound ", length(out), " value(s) for ", ncols, " column(s). A per-column ",
         "function that no longer returns one value per column would be assembled into the ",
         "wrong columns rather than fail, so this is refused.", call. = FALSE)
  }
  out[combat_row_order(idx)]
}


# ---- how a batched loop stays exact ------------------------------------------
#
# Both loop shims below answer out of a batch computed on the first call, and both prove
# the answer belongs to the call before returning it, on EVERY call:
#
#   1. the calling frame is the same environment the batch was built from,
#   2. the matrix and library sizes in that frame are the same OBJECTS the batch was built
#      from, not merely equal to them,
#   3. the loop's own index, read out of that frame, selects the element.
#
# Checks 1 and 2 are pointer comparisons. identical() returns at once when both arguments
# are the same SEXP, so re-verifying a 20000 x 200 matrix on all 200 calls costs nothing
# measurable, while comparing two equal copies of it costs 0.4 ms each. That is the whole
# reason the loop index is used instead of a call counter: a counter would have to be
# defended by comparing the column contents on every call, which measured 30 ms against a
# 46 ms loop and made upperquartile slower than not doing any of this.
#
# The contents ARE compared once, when the batch is built, which is what ties the frame to
# the arguments the vendor is passing down. After that the three checks above carry it.


#' Batch the vendor's per-column TMM loop
#'
#' The vendor writes `for (i in 1:nsamples) f[i] <- .calcFactorTMM(obs = x[, i], ...)`. The
#' loop is the only thing standing between this and a column split, and it cannot be
#' rebound. So the first call computes every column at once, in parallel, using the vendor's
#' own `.calcFactorTMM`, and the loop's remaining calls are served from that.
#'
#' `x`, `lib.size` and `refColumn` are read from the vendor's own frame rather than rebuilt,
#' which is the point: they are the pooled quantities, already computed once on the full
#' matrix by the code that owns them.
#'
#' @param vfun The vendor internal, with `rank` already hoisted.
#' @param loopvar Name of the vendor's loop index, from [calcnorm_backend()].
#' @param what Label used in error messages.
#' @noRd
rp_factor_shim <- function(vfun, loopvar, workers, chunks, parallel_backend, what) {
  seen <- NULL; fx <- NULL; flib <- NULL; frc <- NULL; targs <- NULL; vals <- NULL
  function(obs, ref, libsize.obs = NULL, libsize.ref = NULL, logratioTrim = 0.3,
           sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e10) {
    fr <- parent.frame()
    got <- mget(c(loopvar, "x", "lib.size", "refColumn", "nsamples"), envir = fr,
                inherits = FALSE, ifnotfound = vector("list", 5L))
    i <- got[[loopvar]]

    if (is.null(vals) || !identical(fr, seen) || !identical(got$x, fx) ||
        !identical(got$lib.size, flib) || !identical(got$refColumn, frc) ||
        !identical(targs, list(logratioTrim, sumTrim, doWeighting, Acutoff))) {
      ok <- is.matrix(got$x) && length(got$refColumn) == 1L &&
        identical(ncol(got$x), as.integer(got$nsamples)) &&
        is.numeric(i) && length(i) == 1L && isTRUE(i >= 1) && isTRUE(i <= ncol(got$x)) &&
        identical(got$x[, got$refColumn], ref) &&
        identical(got$lib.size[got$refColumn], libsize.ref) &&
        identical(got$x[, i], obs) && identical(got$lib.size[i], libsize.obs)
      if (!isTRUE(ok)) {
        stop(what, ": the counts, library sizes and reference column in the backend's frame ",
             "are not the ones it passed down, so the column loop this companion batches is ",
             "not the loop edgeR is running. Refusing to return a factor computed against ",
             "the wrong reference.", call. = FALSE)
      }
      seen <<- fr; fx <<- got$x; flib <<- got$lib.size; frc <<- got$refColumn
      targs <<- list(logratioTrim, sumTrim, doWeighting, Acutoff)
      one <- function(j) vfun(obs = fx[, j], ref = ref, libsize.obs = flib[j],
                              libsize.ref = libsize.ref, logratioTrim = logratioTrim,
                              sumTrim = sumTrim, doWeighting = doWeighting,
                              Acutoff = Acutoff)
      vals <<- rp_norm_cols(ncol(fx), length(fx),
                            function(jj) vapply(jj, one, numeric(1)),
                            workers, chunks, parallel_backend,
                            getOption("combat.min.norm.cells", 2e5), what)
    }

    if (!is.numeric(i) || length(i) != 1L || is.na(i) || i != trunc(i) ||
        i < 1L || i > length(vals)) {
      stop(what, ": the backend's loop index `", loopvar, "` is ", deparse(i),
           ", which selects no column of the batch. Refusing to guess which factor was ",
           "asked for.", call. = FALSE)
    }
    vals[[as.integer(i)]]
  }
}

#' Batch the vendor's per-column quantile loop
#'
#' `.calcFactorQuantile` writes
#' `for (j in seq_len(ncol(data))) f[j] <- quantile(data[, j], probs = p)`. Same treatment
#' and same exactness argument as [rp_factor_shim()].
#'
#' This path runs for `method = "upperquartile"` and for the f75 that picks TMM's reference
#' column, so batching it also parallelises a pooled stage.
#'
#' @param q0 The `quantile` the vendor would otherwise have reached.
#' @param loopvar Name of the vendor's loop index, from [calcnorm_backend()].
#' @noRd
rp_quantile_shim <- function(q0, loopvar, workers, chunks, parallel_backend) {
  what <- "calcNormFactors_parallel quantile columns"
  seen <- NULL; fdata <- NULL; vals <- NULL
  function(x, probs, ...) {
    # a multi-element `p` makes the vendor assign a longer value into f[j], which warns and
    # keeps only the first element. Batching that through vapply(numeric(1)) would error
    # where edgeR warns, so it goes straight through instead.
    if (...length() || !is.numeric(probs) || length(probs) != 1L) {
      return(q0(x, probs = probs, ...))
    }
    fr <- parent.frame()
    got <- mget(c(loopvar, "data", "p"), envir = fr, inherits = FALSE,
                ifnotfound = vector("list", 3L))
    j <- got[[loopvar]]

    if (is.null(vals) || !identical(fr, seen) || !identical(got$data, fdata) ||
        !identical(got$p, probs)) {
      ok <- is.matrix(got$data) && identical(got$p, probs) &&
        is.numeric(j) && length(j) == 1L && isTRUE(j >= 1) && isTRUE(j <= ncol(got$data)) &&
        identical(got$data[, j], x)
      if (!isTRUE(ok)) {
        stop(what, ": the matrix and probability in the backend's frame are not the ones it ",
             "passed down, so the column loop this companion batches is not the loop edgeR ",
             "is running.", call. = FALSE)
      }
      seen <<- fr; fdata <<- got$data
      # unname because the vendor assigns into f[j], which drops the "75%" name anyway
      vals <<- rp_norm_cols(ncol(fdata), length(fdata),
                            function(jj) vapply(jj, function(k)
                              unname(q0(fdata[, k], probs = probs)), numeric(1)),
                            workers, chunks, parallel_backend,
                            getOption("combat.min.order.cells", 4e6), what)
    }

    if (!is.numeric(j) || length(j) != 1L || is.na(j) || j != trunc(j) ||
        j < 1L || j > length(vals)) {
      stop(what, ": the backend's loop index `", loopvar, "` is ", deparse(j),
           ", which selects no column of the batch.", call. = FALSE)
    }
    vals[[as.integer(j)]]
  }
}

#' Column-parallel apply, for the one inside .calcFactorRLE
#'
#' `.calcFactorRLE` is `gm <- exp(rowMeans(log(data)))` and then `apply(data, 2, ...)`.
#' Rebinding the `apply` and nothing else leaves `gm` where it belongs, computed once by
#' edgeR across every sample. That is the trap on this path: `gm` is a cross-SAMPLE
#' reduction sitting exactly where the column axis cuts, and rebuilding it per block moved
#' RLE by 0.043. `FUN` closes over the vendor's frame, so each block is still reading the
#' pooled `gm`.
#'
#' @param apply0 The `apply` the vendor would otherwise have reached.
#' @noRd
rp_apply_shim <- function(apply0, workers, chunks, parallel_backend) {
  function(X, MARGIN, FUN, ..., simplify = TRUE) {
    if (!identical(MARGIN, 2) || !is.matrix(X) || ...length() || !isTRUE(simplify)) {
      return(apply0(X, MARGIN, FUN, ..., simplify = simplify))
    }
    rp_norm_cols(ncol(X), length(X),
                 function(jj) apply0(X[, jj, drop = FALSE], 2, FUN),
                 workers, chunks, parallel_backend,
                 getOption("combat.min.order.cells", 4e6),
                 "calcNormFactors_parallel RLE columns")
  }
}


# ---- entry point -------------------------------------------------------------

#' edgeR normalisation factors with the per-column work parallelised
#'
#' Runs edgeR's own `calcNormFactors`. The method is not reimplemented and not copied:
#' `calcNormFactors.default` is called once, unchanged, on the whole matrix, with
#' `.calcFactorTMM`, `.calcFactorTMMwsp`, `.calcFactorRLE` and `.calcFactorQuantile` rebound
#' in a child of edgeR's own namespace. Every other symbol in the body still resolves to
#' edgeR's own code, so `colSums`, the all-zero row filter, the reference column choice, the
#' RLE geometric mean and the final centering are edgeR's, computed once, on every column.
#' Output is `identical()` to `edgeR::calcNormFactors` for the same input, not merely close.
#'
#' @section The two wins:
#' `.calcFactorTMM` builds its trim mask from `rank(logR)` and `rank(absE)` twice each,
#' which is four sorts where two would do. `rank` is rebound to compute each vector once,
#' which is bit-identical and needs no workers at all: measured 1.28s to 0.81s on
#' 20000 x 200, `identical()` TRUE, max abs diff 0. On top of that the vendor's per-column
#' loops are dispatched across workers, taking the same call to 0.29s at four workers.
#'
#' @section Why the axis is columns:
#' Every method here ranks, medians or takes a quantile across genes, so a gene-row block is
#' not a smaller version of the problem. The per-sample factor is independent only once the
#' pooled quantities are fixed, which is why the vendor is called once on the full matrix
#' rather than once per block. Splitting the matrix would also have to rebuild
#' `lib.size <- colSums(x)`, and `colSums` over fractional doubles is not block-associative
#' on arm64, where `sizeof.longdouble` is 8: about 1.8e-14 relative, on the quantity every
#' method divides by.
#'
#' @section Dispatches too small to be worth a fork:
#' A fork costs about 60 ms here, so small matrices are slower parallelised than not, and
#' each path has its own threshold because they do not cost the same per cell. The TMM and
#' TMMwsp loops dispatch at or above `getOption("combat.min.norm.cells", 2e5)` cells:
#' measured 0.50x at 1e5, 1.00x at 2e5, 1.85x at 5e5. The RLE and upperquartile loops take
#' one order statistic per column and cost roughly a fifteenth as much, so they wait for
#' `getOption("combat.min.order.cells", 4e6)`: measured 0.93x at 4e6 and 1.24x at 1e7. The
#' `rank` hoist is not gated, because it never costs anything, and neither is `workers = 1`.
#'
#' @param object Count matrix, genes in rows and samples in columns, or a `DGEList`.
#' @param lib.size Library sizes. Defaults to `colSums(object)`, computed by edgeR on the
#'   whole matrix.
#' @param method Normalisation method, one of `"TMM"`, `"TMMwsp"`, `"RLE"`,
#'   `"upperquartile"` or `"none"`.
#' @param refColumn Reference column for TMM and TMMwsp. Chosen by edgeR when `NULL`.
#' @param logratioTrim,sumTrim,doWeighting,Acutoff TMM trimming and weighting, passed
#'   through unchanged.
#' @param p Quantile for `method = "upperquartile"`.
#' @param ... Ignored, as in edgeR.
#' @param workers Maximum concurrent worker processes. 4 is the default because 6 has
#'   kernel-panicked a 24 GB machine when a second R process was also forking.
#' @param chunks Column chunks per dispatch. Defaults to `workers`.
#' @param parallel_backend One of [combat_backends()], or a function
#'   `function(idx, f, workers)` returning a list in the order of `idx`. Defaults to
#'   `getOption("combat.backend", "mclapply")`.
#' @param backend Optional `calcNormFactors` default method to wrap. Defaults to edgeR's.
#'
#' @return For a matrix, the named numeric vector `edgeR::calcNormFactors` returns. For a
#'   `DGEList`, the object with `$samples$norm.factors` replaced, exactly as
#'   `edgeR::calcNormFactors` returns it.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' counts <- matrix(rnbinom(20000, mu = 50, size = 5), nrow = 1000)
#' identical(calcNormFactors_parallel(counts, workers = 4L),
#'           edgeR::calcNormFactors(counts))
#' }
#' @export
calcNormFactors_parallel <- function(object, lib.size = NULL,
                                     method = c("TMM", "TMMwsp", "RLE", "upperquartile",
                                                "none"),
                                     refColumn = NULL, logratioTrim = 0.3, sumTrim = 0.05,
                                     doWeighting = TRUE, Acutoff = -1e10, p = 0.75, ...,
                                     workers = 4L, chunks = NULL,
                                     parallel_backend = getOption("combat.backend",
                                                                  "mclapply"),
                                     backend = NULL) {
  # no run leaves workers behind, crashed or not; children that predate this call are
  # someone else's and are spared
  .spare <- combat_children()
  on.exit(combat_reap(.spare), add = TRUE)

  if (is.numeric(workers) && length(workers) == 1L && is.finite(workers) &&
      workers != trunc(workers)) {
    stop("`workers` must be a whole number, not ", workers, call. = FALSE)
  }
  workers <- rp_prologue(workers)
  if (!is.function(parallel_backend)) {
    parallel_backend <- match.arg(parallel_backend, combat_backends())
  }

  be <- calcnorm_backend(backend)
  ns <- be$env

  # anything edgeR gives its own S3 method runs different code from the default one this
  # companion wraps, so it is refused rather than quietly funnelled through as.matrix
  if (!inherits(object, "DGEList")) {
    other <- vapply(class(object), function(cl)
      exists(paste0("calcNormFactors.", cl), envir = ns, inherits = FALSE), logical(1))
    if (any(other)) {
      stop("calcNormFactors_parallel handles a matrix or a DGEList. edgeR has its own ",
           "calcNormFactors method for class ", class(object)[which(other)[1]],
           ", which this companion does not wrap. Call edgeR directly for that class, or ",
           "pass its count matrix.", call. = FALSE)
    }
  }

  # the DGEList method is four lines that end in a bare calcNormFactors() call, so it is
  # rebound the same way everything else here is rather than transcribed
  if (inherits(object, "DGEList")) {
    m <- get("calcNormFactors.DGEList", envir = ns, inherits = FALSE)
    denv <- new.env(parent = environment(m))
    denv$calcNormFactors <- function(object, lib.size = NULL, ...) {
      calcNormFactors_parallel(object, lib.size = lib.size, ..., workers = workers,
                               chunks = chunks, parallel_backend = parallel_backend,
                               backend = backend)
    }
    environment(m) <- denv
    return(m(object, method = method, refColumn = refColumn, logratioTrim = logratioTrim,
             sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff, p = p))
  }

  # one child env per rebind target, so the hoisted rank reaches .calcFactorTMM and nothing
  # else, and the parallel apply reaches .calcFactorRLE and nothing else
  tmm <- be$tmm
  environment(tmm) <- local({ e <- new.env(parent = ns); e$rank <- rp_rank_once(be$rank); e })

  quant <- be$quant
  environment(quant) <- local({
    e <- new.env(parent = ns)
    e$quantile <- rp_quantile_shim(be$quantile, be$quant_loop, workers, chunks,
                                   parallel_backend)
    e
  })

  rle <- be$rle
  environment(rle) <- local({
    e <- new.env(parent = ns)
    e$apply <- rp_apply_shim(be$apply, workers, chunks, parallel_backend)
    e
  })

  env <- new.env(parent = ns)
  env$.calcFactorTMM <- rp_factor_shim(tmm, be$tmm_loop, workers, chunks, parallel_backend,
                                       "calcNormFactors_parallel TMM columns")
  env$.calcFactorTMMwsp <- rp_factor_shim(be$wsp, be$tmm_loop, workers, chunks,
                                          parallel_backend,
                                          "calcNormFactors_parallel TMMwsp columns")
  env$.calcFactorQuantile <- quant
  env$.calcFactorRLE <- rle

  f <- be$fn
  environment(f) <- env
  f(object = object, lib.size = lib.size, method = method, refColumn = refColumn,
    logratioTrim = logratioTrim, sumTrim = sumTrim, doWeighting = doWeighting,
    Acutoff = Acutoff, p = p)
}
