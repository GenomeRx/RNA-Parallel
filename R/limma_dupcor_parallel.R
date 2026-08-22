## limma_dupcor_parallel.R
##
## Mirrors ComBat_seq_parallel.R: the entry point, plus the one rebind it needs. Sourcing
## this plus helper_limma_parallel.R and helper_seq_parallel.R is enough to use the function.
##
## limma::duplicateCorrelation is called UNCHANGED on each row block, because its per-gene
## loop reads only M[i, ], weights[i, ], design and the block factor. What this file adds is
## resolving, once and on the full matrix, the two arguments a block would otherwise
## re-derive from its own shape; evaluating the vendor's own pooled tail on the concatenated
## result rather than transcribing it; and memoising the design-invariant call statmod repeats
## once per gene. No vendor source is held here: the vendor's own mixedModel2Fit object runs
## unchanged and only the primitives it calls are shadowed. See rp_dupcor_memo().


#' Memoise the design-invariant work inside duplicateCorrelation's per-gene loop
#'
#' `statmod::mixedModel2Fit` calls `La.svd(QtZ, nu = mq, nv = 0)` once per gene, and `QtZ`
#' reads only the `Z` columns of the effects, so every gene presents byte-identical arguments.
#' Measured: 300 genes, 300 calls, every argument list identical to the first, and that one
#' call is 47.5% of the vendor's wall clock. `factor` and `model.matrix` on the block variable
#' are invariant the same way.
#'
#' @section Why this needs no gate:
#' Nothing here holds vendor source. The vendor's own `mixedModel2Fit` object runs unchanged;
#' only the primitives it calls are shadowed, in a child of the vendor's own environment, by
#' memos keyed on a bitwise `identical()` of the actual arguments. A statmod that stops calling
#' `La.svd`, or calls it with per-gene arguments, produces a cache that never hits: correct and
#' slow, never wrong and quiet. The earlier form evaluated the vendor's statements in a stubbed
#' environment, which could return NA for every gene against a drifted statmod, and needed
#' roughly a hundred lines of body-text pinning to notice. None of that is needed to memoise a
#' call, and none of it is here.
#'
#' `num.eq = FALSE` so 0 and -0 are different keys, `single.NA = FALSE` so NaN and NA are
#' different keys. A changed bit can never return a cached value.
#' @noRd
rp_same <- function(a, b) identical(a, b, num.eq = FALSE, single.NA = FALSE)

#' @noRd
rp_dupcor_memo <- function(fn, env) {
  stat_fn <- tryCatch(get("mixedModel2Fit", envir = env, inherits = TRUE),
                      error = function(e) NULL)
  if (!is.function(stat_fn)) return(NULL)
  senv <- environment(stat_fn)
  if (is.null(senv)) return(NULL)
  real_svd <- tryCatch(get("La.svd", envir = senv, inherits = TRUE),
                       error = function(e) NULL)
  if (!is.function(real_svd)) return(NULL)
  fac <- get("factor", envir = env, inherits = TRUE)
  mm  <- get("model.matrix", envir = env, inherits = TRUE)

  menv <- new.env(parent = senv)
  # La.svd is a pure function of its arguments, so a value keyed on all of them is right
  # whichever of the vendor's two call sites produced it.
  menv$La.svd <- local({
    k <- NULL; v <- NULL; has <- FALSE
    function(x, ...) {
      key <- list(x, ...)
      if (has && rp_same(key, k)) return(v)
      vv <- real_svd(x, ...); k <<- key; v <<- vv; has <<- TRUE
      vv
    }
  })

  fast <- stat_fn                     # the installed vendor object, body untouched
  environment(fast) <- menv

  env2 <- new.env(parent = env)
  env2$mixedModel2Fit <- fast
  env2$factor <- local({
    k <- NULL; v <- NULL; has <- FALSE
    function(x, ...) {
      if (nargs() != 1L) return(fac(x, ...))
      if (has && rp_same(x, k)) return(v)
      vv <- fac(x); k <<- x; v <<- vv; has <<- TRUE
      vv
    }
  })
  env2$model.matrix <- local({
    # limma builds the block design from a formula whose environment holds `A`, captured by
    # reference and rebuilt in that frame every gene, so a formula-keyed cache would serve a
    # stale Z the moment the finite-value mask varies.
    kf <- NULL; ka <- NULL; v <- NULL; has <- FALSE
    function(object, ...) {
      if (nargs() != 1L) return(mm(object, ...))
      a <- tryCatch(get("A", envir = environment(object), inherits = FALSE),
                    error = function(e) NULL)
      if (has && identical(object, kf) && rp_same(a, ka)) return(v)
      vv <- mm(object); kf <<- object; ka <<- a; v <<- vv; has <<- TRUE
      vv
    }
  })

  environment(fn) <- env2
  fn
}


#' Do the vendor's post-loop statements survive being run per block?
#'
#' Everything limma runs after its per-gene loop and before the pooled tail executes inside
#' every block, on that block's own `rho`. Against limma 3.62.2 those seven statements
#' decompose: two set block-invariant constants, two clamp `rho` elementwise, two reduce over
#' genes only to decide whether the clamp is needed at all, and the last is `atanh`. A future
#' limma that made any of them depend on which genes share a block would return a quietly
#' different consensus, and the existing gate two lines below reads only the last two.
#'
#' Pinning their text would put vendor source back in this file, which is the thing 0.4.4
#' removed. This tests the property instead: run them on a probe whole, run them on the same
#' probe split in two, and require the concatenation to agree. That catches a break however it
#' is written, and holds nothing.
#' @noRd
rp_tail_decomposes <- function(fn, ndups, max_block) {
  stmts <- as.list(body(fn))
  fori <- which(vapply(stmts, function(s) is.call(s) &&
                         identical(as.character(s[[1]]), "for"), logical(1)))
  if (!length(fori)) return(FALSE)
  lo <- max(fori) + 1L; hi <- length(stmts) - 2L
  if (lo > hi) return(TRUE)
  mid <- stmts[lo:hi]

  # Deliberately asymmetric: every reduction a tail might take (min, max, max of absolute
  # value, mean, sum, count of non-NA) must differ between the pieces, or a value that leaked
  # out of a guard can agree across the split by coincidence. A first version put -0.995 and
  # 0.995 either side and `rho / max(abs(rho))` slipped through because both halves shared an
  # absolute maximum. Straddles both clamps and carries an NA.
  probe <- c(-0.995, -0.40, 0.10, 0.20, 0.9995, 0.30,
             NA_real_, -0.10, 0.55, -0.72, 0.05, 0.99)
  run <- function(rho) {
    e <- new.env(parent = environment(fn))
    e$rho <- rho; e$ndups <- ndups; e$MaxBlockSize <- max_block; e$block <- TRUE
    for (st in mid) eval(st, e)
    e$arho
  }
  whole <- tryCatch(run(probe), error = function(e) NULL)
  if (is.null(whole)) return(FALSE)
  # several split points, because one split can agree by accident where another does not
  for (cut in c(3L, 5L, 6L, 8L, 10L)) {
    got <- tryCatch(c(run(probe[seq_len(cut)]), run(probe[(cut + 1L):length(probe)])),
                    error = function(e) NULL)
    if (is.null(got) || !identical(whole, got)) return(FALSE)
  }
  TRUE
}


#' Intra-block correlation with the per-gene REML fits parallelised
#'
#' Runs [limma::duplicateCorrelation()] itself. The algorithm is not reimplemented and not
#' copied: the vendor function is called on interleaved row blocks and the per-gene `atanh`
#' correlations are concatenated back into gene order. Output is `identical()` to what
#' `limma::duplicateCorrelation` returns for the same input, not merely close to it.
#'
#' @section Why this one:
#' `statmod::mixedModel2Fit` runs one LAPACK SVD per gene, which measured 86.99% of a
#' blocked limma pipeline. The loop is embarrassingly parallel by gene and nothing else in
#' the body is, so a row split moves nearly all of the wall clock.
#'
#' @section The pooled tail:
#' `consensus.correlation` is `tanh(mean(arho, trim = trim, na.rm = TRUE))`, a trimmed mean
#' over every gene, and a trimmed mean does not decompose over blocks. Copying that one line
#' into this package would create a second copy of it to drift. Instead the vendor's last two
#' top-level statements are lifted out of `body()` and evaluated in a child of the vendor's
#' own environment carrying `arho` and `trim`, so a rewritten limma tail runs as the NEW
#' tail. [all.vars()] checks the lifted slice reads nothing else before anything is
#' dispatched, and the function refuses rather than guessing if it does.
#'
#' @section The per-gene memo:
#' `mixedModel2Fit` reaches `La.svd(QtZ, nu = mq, nv = 0)` once per gene, and `QtZ` reads only
#' the `Z` columns of the effects, so every gene presents byte-identical arguments. Measured:
#' 300 genes, 300 calls, every argument list identical to the first, and that one call is 47.5%
#' of the vendor's wall clock. `mixedModel2Fit` is a bare call in the vendor's loop, so it is
#' rebound in a child of limma's environment onto a copy whose `La.svd` is memoised, along with
#' the `factor` and `model.matrix` calls that rebuild the block design per gene.
#'
#' The memo is installed unconditionally, because it is keyed on the arguments rather than on an
#' assumption about them. Where it pays is unweighted input: counted, 300 genes give 300 calls and
#' one distinct argument set, and it reaches 5.12x serial and 13.55x at four workers on 3000 genes
#' by 100 arrays with 50 blocks. Where it does not pay is weighted input, because
#' `mixedModel2Fit` scales `X` by each gene's own weights before `La.svd` sees it, so the argument
#' is not invariant: 300 calls, 300 distinct sets, a 0% hit rate. That costs nothing beyond one
#' failed key comparison per gene, and it is the same fact that made the older body-pinned lift
#' refuse to install on weighted input rather than run and miss.
#'
#' @section What is refused:
#' `block = NULL` errors. That path pairs rows through `unwrapdups`, and
#' `if (spacing == "topbottom") spacing <- nrow(M)/2` is a literal branch on the block's own
#' row count, so an interleaved split would pair different genes and raise nothing.
#'
#' Weights are expanded to a full matrix once, before splitting, because
#' `limma::asMatrixWeights` dispatches on the BLOCK's row count and tests its gene branch
#' before its array branch. A bare length-`narrays` weight vector is therefore read as gene
#' weights by any block whose `nrow` happens to equal `narrays`; measured, that moved
#' `atanh.correlations` by 1.56. See `rp_weights_matrix()`.
#'
#' @section Backends for the parallelism:
#' `parallel_backend` selects the framework, exactly as in [ComBat_seq_parallel()]. All of
#' them return bit-identical results. Below `getOption("combat.min.dupcor.cells", 5000)`
#' gene-by-array cells the call runs serially: one gene is one REML fit, so the fork pays at
#' far smaller sizes than the matrix paths, but not at every size.
#'
#' @param object A matrix, `EList`, `MAList` or `ExpressionSet`, as the vendor takes it.
#' @param design Design matrix. Resolved once from the object when `NULL`, because a block
#'   is handed a bare matrix and would fall back to an intercept instead.
#' @param ndups,spacing Duplicate-spot arguments. Kept so the signature mirrors the vendor;
#'   both are only read on the `block = NULL` path, which this function refuses.
#' @param block Blocking factor, one entry per column. Required.
#' @param trim Fraction trimmed from each end of the pooled mean.
#' @param weights Gene, array or full weight matrix. Expanded once against the full
#'   dimensions before splitting.
#' @param workers Maximum concurrent worker processes. 4 is the default because 6 has
#'   kernel-panicked a 24 GB machine when a second R process was also forking.
#' @param chunks Row chunks. Defaults to `workers`. One-gene chunks are allowed here:
#'   the `lm.fit` one-column demotion behind `lmFit_parallel`'s two-gene floor never
#'   applies, because `duplicateCorrelation` does not reach `lm.fit` per gene.
#' @param parallel_backend One of [combat_backends()], or a function
#'   `function(idx, f, workers)` returning a list in the order of `idx`.
#' @param backend Optional `duplicateCorrelation` function to wrap. Defaults to
#'   `limma::duplicateCorrelation`.
#'
#' @return `list(consensus.correlation, cor, atanh.correlations)`, `identical()` to what the
#'   backend returns for the same input.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("limma", quietly = TRUE) &&
#'     requireNamespace("statmod", quietly = TRUE)) {
#'   set.seed(1)
#'   M <- matrix(rnorm(200 * 12), 200, 12)
#'   subject <- rep(1:6, each = 2)
#'   design <- model.matrix(~ rep(c(0, 1), 6))
#'
#'   dc <- duplicateCorrelation_parallel(M, design, block = subject, workers = 4L)
#'   identical(dc, limma::duplicateCorrelation(M, design, block = subject))
#' }
#' }
#' @export
duplicateCorrelation_parallel <- function(object, design = NULL, ndups = 2L, spacing = 1L,
                                          block = NULL, trim = 0.15, weights = NULL,
                                          workers = 4L, chunks = NULL,
                                          parallel_backend = getOption("combat.backend", "mclapply"),
                                          backend = NULL) {
  # no run leaves workers behind, crashed or not; children that predate this call are
  # someone else's and are spared
  .spare <- combat_children()
  on.exit(combat_reap(.spare), add = TRUE)

  workers <- rp_prologue(workers)
  if (!is.function(parallel_backend)) {
    parallel_backend <- match.arg(parallel_backend, combat_backends())
  }

  if (is.null(block)) {
    stop("`block` is required. With block = NULL the vendor pairs rows through unwrapdups, ",
         "and `if (spacing == \"topbottom\") spacing <- nrow(M)/2` branches on the block's ",
         "own row count, so a row split would pair different genes and raise nothing. ",
         "Call limma::duplicateCorrelation directly for the ndups path.", call. = FALSE)
  }

  if (is.null(backend)) {
    if (!requireNamespace("limma", quietly = TRUE)) {
      stop("no limma backend found. Install it with BiocManager::install(\"limma\").",
           call. = FALSE)
    }
    backend <- limma::duplicateCorrelation
  }
  be <- limma_backend(backend,
                      need_args = c("object", "design", "ndups", "spacing", "block",
                                    "trim", "weights"),
                      # Only the rebinds whose loss would be total. rp_dupcor_memo also
                      # shadows factor and model.matrix, deliberately left out: if a future
                      # limma writes as.factor(block) those two memos stop reaching and the
                      # call loses a little speed, where mixedModel2Fit, getEAWP or
                      # asMatrixWeights going unreachable would leave the companion doing
                      # nothing at all at serial pace, which no equivalence test can see.
                      rebound = c("getEAWP", "asMatrixWeights", "mixedModel2Fit"))

  # Lifted and checked BEFORE any work is dispatched, so a limma this package cannot
  # reproduce costs an error rather than an hour of SVDs and then an error.
  stmts <- as.list(body(be$fn))
  tail_exprs <- stmts[c(length(stmts) - 1L, length(stmts))]
  extra <- setdiff(unlist(lapply(tail_exprs, all.vars)), c("arho", "trim", "mrho"))
  if (length(extra)) {
    stop("this limma backend's last two statements read ", paste(extra, collapse = ", "),
         ", which a pooled tail cannot supply. rnaparallel was written against a tail that ",
         "reads only the concatenated atanh correlations and `trim`. Refusing to run.",
         call. = FALSE)
  }
  # The statements between the loop and that tail run inside every block on its own rho, and
  # the check above never looked at them. Refuse the split unless they still decompose.
  if (!rp_tail_decomposes(be$fn, ndups, max(table(block)))) {
    stop("this limma backend's statements between the per-gene loop and the pooled tail no ",
         "longer give the same result run per block as run over all genes, so a split would ",
         "return a quietly different consensus correlation. Call limma::duplicateCorrelation ",
         "directly.", call. = FALSE)
  }

  eawp <- get("getEAWP", envir = be$env, inherits = TRUE)(object)
  M <- eawp$exprs
  # A block is handed a bare matrix, so getEAWP gives it design NULL and it would silently
  # substitute an intercept. matrix(1, ncol(M), 1) is the vendor's own fallback and is the
  # same in every block, so resolving here changes nothing except the y$design case.
  if (is.null(design)) design <- eawp$design
  design <- if (is.null(design)) matrix(1, ncol(M), 1) else as.matrix(design)
  if (is.null(weights)) weights <- eawp$weights
  # The vendor has degenerate early returns it reaches WITHOUT ever looking at weights, so
  # expanding first turned inputs limma answers into errors here. Fall back rather than
  # refuse: a weight object this package cannot account for is a reason to run serially,
  # not a reason to fail where the vendor succeeds.
  w_raw <- weights
  weights <- tryCatch(rp_weights_matrix(weights, dim(M), be$env),
                      error = function(e) NULL)
  if (!is.null(w_raw) && is.null(weights)) {
    return(be$fn(object = object, design = design, ndups = ndups, spacing = spacing,
                 block = block, trim = trim, weights = w_raw))
  }

  # No input gate. The memo keys on the actual arguments, so a varying finite-value mask
  # misses the cache and runs the real call rather than returning a stale one, and it still
  # pays on the weighted and non-finite inputs the old body-pinned lift refused outright.
  dc <- be$fn
  memo <- rp_dupcor_memo(be$fn, be$env)
  if (!is.null(memo)) dc <- memo

  ngenes <- nrow(M)
  # min_rows = 1: hazard 3 is an lm.fit effect and duplicateCorrelation never reaches
  # lm.fit, so clamping to 2 here only cost parallelism on tiny inputs.
  idx <- combat_row_chunks(ngenes, workers = workers, chunks = chunks, min_rows = 1L)

  fit_block <- function(ii) {
    conds <- list()
    keep <- function(cnd) conds[[length(conds) + 1L]] <<- cnd
    value <- withCallingHandlers(
      dc(object = M[ii, , drop = FALSE], design = design, ndups = ndups,
         spacing = spacing, block = block, trim = trim,
         weights = rp_weights_rows(weights, ii)),
      warning = function(w) { keep(w); invokeRestart("muffleWarning") },
      message = function(m) { keep(m); invokeRestart("muffleMessage") })
    list(value = value, conds = conds)
  }

  # One gene is one REML fit, so this earns a fork at far smaller sizes than the row-split
  # paths do, but not at every size. Counted in cells rather than genes because the per-gene
  # REML scales with the arrays too, and a hundred genes across a thousand arrays is real work
  # a gene-count gate sent serial. Cells bound the decision only loosely below 5,000, where
  # 300x6 gains 2.68x and 40x50 loses at a similar cell count, so the gate sits above the mixed
  # band rather than at a crossover: measured on 4 workers, 0.46x at 25x24, 0.80x at 40x50 and
  # 1.00x at 60x50, against 2.76x at 50x100, 4.39x at 2000x6 and 4.41x at 100x100.
  # length(M) rather than a product, which overflows integer.
  parts <- combat_parallel_check(
    combat_parallel_lapply(idx, fit_block, workers, parallel_backend,
                           cells = length(M),
                           min_cells = getOption("combat.min.dupcor.cells", 5000)),
    "duplicateCorrelation across gene blocks", idx)

  # The rank note and the two degenerate-block warnings are raised once per block, and some
  # backends swallow child output entirely, so what the caller sees would otherwise depend
  # on parallel_backend. Replay each distinct condition once, from here.
  cnds <- unlist(lapply(parts, function(p) p$conds), recursive = FALSE)
  for (cnd in cnds[!duplicated(vapply(cnds, conditionMessage, ""))]) {
    cnd$call <- NULL                     # else the block's own call prints, M[ii, ] and all
    if (inherits(cnd, "warning")) warning(cnd) else message(cnd)
  }

  vals <- lapply(parts, function(p) p$value)
  # combat_parallel_check skips its height check on list results, and the vendor returns a
  # list, so the per-block count is checked here instead. A short block plus a long one can
  # still total ngenes.
  got <- vapply(vals, function(v) length(v$atanh.correlations), integer(1))
  if (!identical(got, lengths(idx))) {
    stop("duplicateCorrelation_parallel: blocks returned ", paste(got, collapse = "/"),
         " correlation(s) where ", paste(lengths(idx), collapse = "/"),
         " gene(s) were dispatched.", call. = FALSE)
  }
  arho <- rp_bind_rows(vals, combat_row_order(idx), "atanh.correlations", ngenes,
                       "duplicateCorrelation_parallel")

  tenv <- new.env(parent = be$env)
  tenv$arho <- arho
  tenv$trim <- trim
  out <- NULL
  for (e in tail_exprs) out <- eval(e, tenv)
  if (!is.list(out)) {
    stop("the lifted limma tail returned a ", class(out)[1],
         " rather than the result list. Refusing to return it.", call. = FALSE)
  }
  out
}
