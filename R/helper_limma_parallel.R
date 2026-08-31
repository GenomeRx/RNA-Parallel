## helper_limma_parallel.R
##
## Everything the limma companions need that is not an entry point. Mirrors the role
## helper_seq_parallel.R plays for ComBat-seq, and reuses its dispatch layer unchanged:
## combat_row_chunks, combat_parallel_lapply, combat_parallel_check, combat_row_order.
##
## Nothing here reimplements limma. The original function is called unchanged on each block.
## What this file does is resolve, once and on the full matrix, the arguments a block would
## otherwise reinterpret, and refuse to split when a block would take a different branch
## than the whole matrix would.


# ---- the three ways a limma row split goes silently wrong ---------------------
#
# All three were measured against limma 3.62.2, not reasoned about.
#
# 1. asMatrixWeights dispatches on the BLOCK's row count. Its gene branch
#    (lw == 1 || lw == dim[1]) is tested BEFORE its array branch (lw == dim[2]), so a bare
#    weight vector whose length equals a block's row count changes meaning mid-run. A
#    length-8 array-weight vector read as gene weights in an 8-row block moved coefficients
#    by 0.325 and atanh.correlations by 1.56. Not a last-place difference, a different model.
#
# 2. NoProbeWts is an AND-reduction over every cell that selects between two numerically
#    different algorithms: the whole-matrix lm.fit whose sigma comes from colMeans over the
#    QR effects, and a per-gene loop whose sigma comes from mean() over one gene's effects.
#    Measured: 114 of 400 sigma differ at 4.44e-16 when a block flips. The two branches also
#    return DIFFERENT component sets, since the fast path returns lm.fit's whole object
#    including qr and assign and the slow path returns a bare list without them.
#
# 3. stats::lm.fit contains `if (is.matrix(y) && ny == 1L) y <- drop(y)`. limma feeds it
#    t(M), so a one-gene block presents ny == 1, fit$effects loses its dim, and sigma
#    switches from colMeans (single long-double pass) to mean (two-pass refinement).
#    Measured: 4080 of 16000 one-gene blocks not identical() to serial, and the gene name is
#    lost as well. Every caller here therefore passes min_rows = 2.


#' Expand weights to a full matrix before anything is split
#'
#' A matrix has no shape left to re-dispatch on, which is what makes this the fix for
#' hazard 1 above. `arrayweights` is carried explicitly because [limma::asMatrixWeights()]
#' sets it to mark an expanded per-array vector, and `NoProbeWts` reads it.
#'
#' @param weights `NULL`, a vector, or a matrix, as the user supplied it.
#' @param dim_full `dim()` of the FULL matrix, never a block's.
#' @param env Environment holding limma's `asMatrixWeights`, taken from the backend so the
#'   pair can never be drawn from two different copies of limma.
#' @return `NULL`, or a `dim_full` matrix carrying limma's own `arrayweights` attribute.
#' @noRd
rp_weights_matrix <- function(weights, dim_full, env) {
  if (is.null(weights)) return(NULL)
  amw <- get("asMatrixWeights", envir = env, inherits = TRUE)
  w <- amw(weights, dim_full)
  if (!identical(dim(w), as.integer(dim_full))) {
    stop("asMatrixWeights returned a ", paste(dim(w), collapse = " x "), " matrix for a ",
         paste(dim_full, collapse = " x "), " input. Refusing to split against a weight ",
         "object this package cannot account for.", call. = FALSE)
  }
  w
}


#' Subset an expanded weight matrix without losing what the branch reads
#'
#' `[` drops attributes, and dropping `arrayweights` flips `NoProbeWts` for that block.
#' @noRd
rp_weights_rows <- function(w, ii) {
  if (is.null(w)) return(NULL)
  aw <- attr(w, "arrayweights")
  out <- w[ii, , drop = FALSE]
  if (!is.null(aw)) attr(out, "arrayweights") <- aw
  out
}


#' Decide whether a row split can reproduce the branch the whole matrix takes
#'
#' Returns TRUE when every block provably lands on the same side of `NoProbeWts` as the
#' full matrix, reading the original's own guard expression rather than approximating it.
#' Non-positive weights become NA and are punched into M first, so the finiteness test is
#' taken on the matrix limma will actually branch on.
#'
#' The `gls` punch is deliberately a superset of the original's: gls.series punches on
#' `w < 1e-15` after mapping NA to 0, and this punches non-finite weights as well. A `+Inf`
#' weight therefore reports the split unsafe where the original would have tolerated it, which
#' costs a serial run on an input nobody has. The alternative is a guard that is right by
#' coincidence.
#'
#' @return `TRUE` when splitting is safe, `FALSE` when the caller must run serially.
#' @noRd
rp_branch_stable <- function(M, weights, punch = c("lm", "gls")) {
  punch <- match.arg(punch)

  # Second conjunct FALSE means NoProbeWts is FALSE for the full matrix AND for every
  # block regardless of finiteness, so the branch cannot flip and NAs are harmless.
  if (!is.null(weights) && is.null(attr(weights, "arrayweights"))) return(TRUE)

  # Otherwise NoProbeWts reduces to all(is.finite(M)). All finite means every block is
  # finite too. Any non-finite cell means some block may be wholly finite and flip.
  # The two callers punch weights into M by DIFFERENT rules, and reading the wrong one
  # here reports a split as safe when it is not: a weight strictly inside (0, 1e-15) is
  # punched by gls.series and not by lm.series, which put ten of forty sigma values on the
  # wrong branch with no error and no warning.
  # The punch exists only to make a cell non-finite, and the answer is `all cells finite`, so
  # the weights can be reduced instead of materialised. `w[w <= 0] <- NA` copied the weights,
  # `is.finite(w)` built a second full logical, its negation a third, and the subassignment
  # copied M as well: five full-size allocations to decide one TRUE or FALSE. min() and max()
  # report every failure mode -- NA for any NA, NaN for any NaN, an infinity for either
  # infinity -- so the two rules reduce exactly. lm passes a weight iff it is finite and > 0;
  # gls maps NA to 0 first and punches below 1e-15, so it passes iff finite and >= 1e-15.
  # Measured on array weights: 20,000 x 24 2.13 ms to 1.25 ms, 200,000 x 48 40.85 ms to
  # 27.04 ms on the lm punch and 77.60 ms to 25.76 ms on the gls one.
  if (!length(M)) return(TRUE)     # min() of nothing is +Inf and would flip the answer
  if (!is.null(weights)) {
    lo <- suppressWarnings(min(weights)); hi <- suppressWarnings(max(weights))
    if (!(is.finite(lo) && is.finite(hi))) return(FALSE)
    if (punch == "lm") { if (!(lo > 0)) return(FALSE) } else if (!(lo >= 1e-15)) return(FALSE)
  }
  # rowSums is non-finite for every row holding a non-finite cell, so the cheap scan can
  # only over-report. The exact scan settles the one case it invents, an overflowing sum.
  (is.numeric(M) && all(is.finite(rowSums(M)))) || all(is.finite(M))
}


#' Is an array-weight matrix actually constant down its rows
#'
#' limma's fast branch reads only `weights[1, ]`, so a matrix carrying the `arrayweights`
#' attribute but varying by row gives each block its OWN first row and a different qr. The
#' attribute survives both `*` and `t()`, so `asMatrixWeights(aw, dim(M)) * probe_weights`
#' produces this shape from exported limma calls alone. limma's own voomaLmFit strips the
#' attribute after exactly that multiply. Callers fall back to serial rather than refusing.
#' `rows` is the row each block leads with, the only rows the fast branch ever reads, so a
#' whole-matrix scan would charge every call for rows no block can reach.
#' @noRd
rp_arrayweights_uniform <- function(w, rows) {
  if (is.null(w) || is.null(attr(w, "arrayweights"))) return(TRUE)
  first <- unname(w[1L, ])
  all(vapply(rows, function(i) identical(unname(w[i, ]), first), logical(1)))
}


#' Assert that fields lifted from one block are the same in all of them
#'
#' `qr`, `assign`, `rank`, `pivot` and `cov.coefficients` are functions of the design alone
#' and are therefore taken from the first block. "Therefore" is doing a lot of work in that
#' sentence: a rank-deficient design plus a branch flip produced pivot 1 2 3 against 1 3 2
#' and a 3x3 cov.coefficients against a 2x2 one. Checking is cheap; being wrong is silent.
#' @noRd
rp_invariant <- function(parts, fields, what) {
  if (length(parts) < 2L) return(invisible(TRUE))
  first <- parts[[1L]]

  # Blocks that came back with different component sets ran different algorithms, whatever
  # the surviving fields say. Checked before the field loop because that loop intersects
  # against block 1's names: when block 1 is the slow one, qr and assign are simply absent
  # from the intersection and the mismatch is never compared.
  nm1 <- names(first)
  for (k in seq_along(parts)[-1L]) {
    if (!identical(names(parts[[k]]), nm1)) {
      stop(what, ": block 1 returned (", paste(nm1, collapse = ", "), ") and block ", k,
           " returned (", paste(names(parts[[k]]), collapse = ", "), "). Different component ",
           "sets mean the blocks took different branches inside the original. Refusing to ",
           "assemble one result out of two algorithms.", call. = FALSE)
    }
  }
  for (nm in fields) {
    ref <- first[[nm]]
    for (k in seq_along(parts)[-1L]) {
      if (!identical(parts[[k]][[nm]], ref)) {
        stop(what, ": '", nm, "' differs between blocks, so it is not the design-only ",
             "quantity this package assumes. Block 1 and block ", k, " disagree. ",
             "Refusing to assemble a result from a field that is not block-invariant.",
             call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}


#' Bind row-blocks back together, preserving names
#'
#' `unlist(use.names = FALSE)` on a vector field is the obvious way to do this and it is
#' wrong: limma's fast path names `sigma` from `colnames(t(M))`, so dropping names makes
#' identical() FALSE while every value is bit-equal. Reported as a divergence once already.
#' @noRd
rp_bind_rows <- function(parts, ord, nm, ngenes, what) {
  pieces <- lapply(parts, function(p) p[[nm]])
  if (is.null(pieces[[1L]])) return(NULL)
  if (is.matrix(pieces[[1L]])) {
    dn <- dimnames(pieces[[1L]])
    m <- do.call(rbind, pieces)
    if (nrow(m) != ngenes) {
      stop(what, ": field '", nm, "' bound to ", nrow(m), " rows where ", ngenes,
           " were dispatched", call. = FALSE)
    }
    m <- m[ord, , drop = FALSE]
    # rbind decides dimnames for itself and collapses list(NULL, NULL) to a plain NULL,
    # which is a different object even when every value is bit-equal. gls.series has no
    # coef.names fallback, so an unnamed matrix with an unnamed design reaches exactly that
    # shape, and limma's own implicit design matrix(1, narrays, 1) has no colnames. Measured
    # at default settings in 127 of 2201 random configurations before this was restored.
    if (!is.null(dn)) {
      dimnames(m) <- list(rownames(m), dn[[2L]])
      if (!is.null(names(dn))) names(dimnames(m)) <- names(dn)
    }
    return(m)
  }
  v <- unlist(pieces)                       # names kept on purpose, see above
  if (length(v) != ngenes) {
    stop(what, ": field '", nm, "' bound to ", length(v), " values where ", ngenes,
         " were dispatched", call. = FALSE)
  }
  v[ord]
}


#' Resolve a limma backend and refuse when a rebind target has moved
#'
#' The limma analogue of `combat_backend()`. An unreachable rebind returns correct output
#' at serial speed, and no equivalence test can see the difference, so this refuses to
#' start rather than run a companion that silently does nothing.
#'
#' @param fn A limma function. Defaults to `limma::lmFit`.
#' @param need_args Formals the caller was written against.
#' @param rebound Symbols that must still be reachable as bare calls in `fn`'s body.
#' @return A list with `fn` and `env`.
#' @noRd
limma_backend <- function(fn = NULL, need_args = character(), rebound = character()) {
  if (is.null(fn)) {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("no limma backend found. Install it with BiocManager::install(\"limma\").",
           call. = FALSE)
    }
    fn <- limma::lmFit
  }
  if (!is.function(fn)) stop("`backend` must be a function", call. = FALSE)
  env <- environment(fn)
  if (is.null(env)) {
    stop("the limma backend has no environment, so its helpers cannot be resolved.",
         call. = FALSE)
  }
  missing_args <- setdiff(need_args, names(formals(fn)))
  if (length(missing_args)) {
    stop("this limma backend is missing argument(s): ", paste(missing_args, collapse = ", "),
         ". rnaparallel was written against limma 3.62's signature.", call. = FALSE)
  }
  # Same call-head walk as combat_backend: heads only, and only bare symbols, so a future
  # `limma::lm.series(...)` is caught rather than being read as reachable.
  unreachable <- setdiff(rebound, rp_bare_call_heads(body(fn)))
  if (length(unreachable)) {
    stop("this limma backend no longer calls ", paste(unreachable, collapse = ", "),
         " as a bare symbol, so rebinding cannot reach it. The package would run serially ",
         "while still returning identical() output, which no equivalence test can detect. ",
         "Refusing to run.", call. = FALSE)
  }
  list(fn = fn, env = env)
}
