## limma_lmFit_parallel.R
##
## Mirrors limma's lmFit.R: the entry point and nothing else. Source this plus
## helper_limma_parallel.R and helper_seq_parallel.R to use the function without
## installing anything.
##
## The core function is still limma's. This file does not contain a copy of lmFit; it
## calls limma's own lmFit with two symbols rebound.


#' Split one limma row-fitter across gene blocks
#'
#' The shared body of both rebinds. `block_fn` gets one block's matrix and that block's
#' slice of the once-expanded weights. `serial_fn` is the untouched whole-matrix vendor
#' call, taken whenever a split cannot be proved exact or cannot pay for itself.
#'
#' Only four fields vary by gene. Everything else either comes from the design alone or is
#' an argument echoed back, so it is lifted from the first block and then asserted across
#' all of them.
#'
#' @noRd
rp_row_blocks <- function(M, weights, env, workers, chunks, parallel_backend, what,
                          block_fn, serial_fn, punch = "lm") {
  M <- as.matrix(M)
  w <- rp_weights_matrix(weights, dim(M), env)

  # SIZE GATE FIRST. `fast` reads only is.null(w) and attr(w, "arrayweights"), neither of
  # which rp_branch_stable touches, so it is the same value in either order; and both guards
  # return the same expression, serial_fn(). What changes is that a call destined for the
  # vendor no longer pays a full-matrix scan to get there. It was paying enough to lose:
  # measured on array-weighted input under the gate, the companion was SLOWER than the
  # function it wraps -- 20,000 x 24 vendor 10.8 ms against 14.8 ms, 60,000 x 48 vendor
  # 89.3 ms against 138.4 ms. On a platform where the gate is shut this scan ran on every
  # call and the split never followed it.
  fast <- is.null(w) || !is.null(attr(w, "arrayweights"))
  # Each branch consults its OWN option, and both close where the payload would be copied.
  # Routing both through one option let a raised combat.min.ls.cells silently switch off the
  # weighted branch's split as well; see rp_ls_min_cells().
  min_cells <- if (fast) rp_ls_min_cells("combat.min.ls.cells", 6e6, parallel_backend)
               else      rp_ls_min_cells("combat.min.cells",    2e4, parallel_backend)

  # Under the gate, take the vendor call WHOLE. combat_parallel_lapply honours the gate by
  # walking the blocks serially instead, and on the fast branch that is four lm.fit calls
  # where one would do: measured 0.70x, a companion slower than the function it wraps. An
  # unusable option value falls through to the dispatch, which refuses it there.
  if (isTRUE(length(M) < suppressWarnings(as.numeric(min_cells)))) return(serial_fn())

  # A block landing on the other side of NoProbeWts returns a different component SET, not
  # merely different numbers, so there would be nothing to reassemble. lm.series and
  # gls.series punch weights into M by different rules, so the guard is told which.
  if (!rp_branch_stable(M, w, punch)) return(serial_fn())

  # min_rows = 2 is exactness, not tuning: lm.fit drops a one-column response to a vector.
  idx <- combat_row_chunks(nrow(M), workers, chunks, min_rows = 2L)
  if (length(idx) < 2L) return(serial_fn())

  # The fast branch reads weights[1, ] only. A matrix flagged arrayweights but varying by
  # row therefore gives each block its own first row and a different qr. Reachable from
  # exported calls alone, since the attribute survives `*`: asMatrixWeights(aw, dim(M)) *
  # probe_weights keeps it. limma's own voomaLmFit strips it after exactly that multiply.
  # Checked once the blocks are known, so it reads the leading rows and not the whole matrix.
  if (!rp_arrayweights_uniform(w, vapply(idx, `[[`, integer(1), 1L))) return(serial_fn())

  # A dispatched closure is SERIALISED on a socket backend, and a closure carries its whole
  # defining environment whether the body reads it or not. This frame reaches `M` a second
  # time and, through block_fn, the entry point's own `object`, so the dispatch presented
  # 664.29 MiB of globals for a 55 MiB matrix and `future` refused it outright:
  #   "The total size of the 12 globals exported ... is 664.29 MiB. This exceeds the maximum"
  # Rebuilding against an environment holding only the four objects the body reads leaves the
  # matrix and its weights to travel and nothing else. On a forking backend this is invisible
  # either way, since the child inherits the pages, so the cost was only ever paid where there
  # is no fork() -- which is where the adaptive default sends a Windows caller. The same fix
  # is at edger_norm_parallel.R for the TMM column loop. Nothing about the arithmetic changes.
  lean <- new.env(parent = globalenv())
  lean$M <- M; lean$w <- w; lean$block_fn <- block_fn
  lean$rp_weights_rows <- rp_weights_rows
  per_block <- function(ii) block_fn(M[ii, , drop = FALSE], rp_weights_rows(w, ii))
  environment(per_block) <- lean

  parts <- combat_parallel_check(
    combat_parallel_lapply(
      idx, per_block,
      workers, parallel_backend, cells = length(M), min_cells = min_cells),
    what, idx)

  rp_invariant(parts,
               intersect(c("qr", "assign", "rank", "pivot", "cov.coefficients",
                           "ndups", "spacing", "block", "correlation"),
                         names(parts[[1L]])),
               what)

  ord <- combat_row_order(idx)
  out <- parts[[1L]]
  for (nm in c("coefficients", "stdev.unscaled", "sigma", "df.residual")) {
    out[[nm]] <- rp_bind_rows(parts, ord, nm, nrow(M), what)
  }
  out
}


#' limma's lmFit with its row loops parallelised
#'
#' Runs [limma::lmFit()] itself. The algorithm is not reimplemented and not copied: limma's
#' own `lmFit` is called with `lm.series` and `gls.series` rebound in a child of its own
#' environment, so those two resolve to row-parallel versions while `getEAWP`,
#' `nonEstimable`, `asMatrixWeights`, `duplicateCorrelation` and `new("MArrayLM", ...)`
#' still resolve to limma's. Output is `identical()` to what `limma::lmFit` returns for the
#' same input, on the whole `MArrayLM`, not on selected slots.
#'
#' `lm.series` and `gls.series` are where a voom pipeline actually spends its time:
#' 19.25% of a plain one and 12.80% of a blocked one.
#'
#' @section What a row block reinterprets, and how that is handled:
#' A limma row fitter reads three things off the shape of the matrix it is handed, so a
#' block can quietly answer a different question than the whole matrix would. All three
#' were measured against limma 3.62.2.
#'
#' `asMatrixWeights` tests its gene branch (`lw == 1 || lw == dim[1]`) before its array
#' branch (`lw == dim[2]`), so a bare length-narrays weight vector is read as GENE weights
#' by any block whose row count happens to equal the array count. Splitting a 40 by 8
#' matrix into 8-row blocks moved coefficients by 0.325. Weights are therefore expanded to a
#' full matrix ONCE, against the full dimensions, and row-subset with limma's own
#' `arrayweights` attribute carried along.
#'
#' `NoProbeWts` selects between two numerically different algorithms that also return
#' different component sets: the fast path returns `lm.fit`'s whole object including `qr`
#' and `assign`, the slow path a bare seven-element list without them. One NA cell in a 400
#' gene matrix put the whole matrix on the slow path and every block on the fast one, and
#' 114 of 400 sigma differed. Splitting is refused, in favour of one plain vendor call,
#' unless every block provably lands on the same side as the full matrix.
#'
#' `stats::lm.fit` demotes a one-column response with `if (is.matrix(y) && ny == 1L)
#' y <- drop(y)`, and limma feeds it `t(M)`, so a one-gene block switches sigma from
#' `colMeans` to `mean` and loses the gene name: 4,080 of 16,000 one-gene blocks returned a
#' different sigma and all 2000 lost the name. Every split here holds at least two genes.
#'
#' @section Consensus correlation:
#' `correlation` has no default in `lmFit` and the blocked branch tests
#' `missing(correlation)`, so missingness is forwarded rather than turned into `NULL`. An
#' explicit `correlation = NULL` does reach `gls.series`, which then calls
#' `duplicateCorrelation` itself, a trimmed mean over whatever genes it was given. Per
#' block that is a different consensus for every block. It is resolved once here, on the
#' full matrix, and the resulting scalar is what the blocks see.
#'
#' @section Dispatches too small to be worth a fork:
#' The two branches of `lm.series` are not the same job. With probe weights, or an `EList`
#' from [limma::voom()], limma runs an R loop over genes, and that forks well: measured
#' 2.52x unblocked and 2.85x blocked on 20,000 genes by 24 samples, 2.96x and 3.39x on
#' 60,000 by 48, against 0.95x at 20,000 cells. `getOption("combat.min.cells", 20000)`, the
#' same gate ComBat-seq uses, lands on that break-even and is what this branch takes.
#'
#' Without probe weights limma fits every gene in one vectorised `lm.fit`, which costs
#' milliseconds and is mostly not worth handing to another process. That branch measured
#' 0.24x at 1.2 million cells and 0.56x at 4.8 million, turning over to 1.47x at 12 million
#' and 1.71x at 20 million. It gets its own, much higher
#' `getOption("combat.min.ls.cells", 6e6)`; sharing one gate made it four times slower than
#' calling limma directly. Cells is only a proxy, since the cost rises with array count
#' too, so the gate sits at the smallest size measured to win. Set either option to 0 to
#' dispatch unconditionally; output is identical either way.
#'
#' @section Configurations this refuses:
#' `method = "robust"` goes to `mrlm`, and `ndups >= 2` makes `unwrapdups` reshape rows so
#' that a gene no longer occupies one row. Neither is split, and both stop with an error
#' rather than returning a serial result that looks parallel.
#'
#' @section When this is worth reaching for:
#' It depends entirely on which branch your call takes, and the two are not close.
#'
#' With voom or probe weights limma runs an R loop over genes, which forks well but has a
#' floor. Measured on an M3 at the default worker count, companion against vendor, every arm
#' `identical()`: 0.59x at 1,000 genes by 24 arrays, 0.87x at 2,000 x 24, 1.29x at 4,000 x 24,
#' 1.76x at 8,000 x 24, and 2.79x at 60,000 x 48. The crossover tracks gene count rather than
#' cells, and sits near four thousand genes.
#'
#' Without probe weights limma fits every gene in one vectorised `lm.fit`, which costs
#' milliseconds, so the companion declines to split until the matrix is very large and is
#' parity until then. That is not a missing measurement, it is the gate doing its job.
#'
#' Under either gate the companion is one plain vendor call plus about half a millisecond.
#' Nothing on a single call; worth avoiding in a loop over thousands of small fits, which is
#' the one case a per-call size gate cannot help with.
#'
#' @param object Anything `limma::lmFit` accepts: matrix, `EList` from
#'   [limma::voom()], `MAList`, `ExpressionSet`, `PLMset`, `marrayNorm`, numeric data frame.
#'   `EListRaw` and `RGList` are refused by limma itself with \"please normalize first\".
#' @param design Design matrix, samples in rows.
#' @param ndups Rows containing each gene. Only 1 is supported, see refusals above.
#' @param spacing Spacing between duplicate rows.
#' @param block Blocking vector, one entry per array. Non-`NULL` routes to `gls.series`.
#' @param correlation Inter-duplicate or inter-block correlation. No default, exactly as in
#'   `lmFit`; leave it missing and the blocked branch stops the way `lmFit` does.
#' @param weights Numeric matrix of per-observation weights, a per-gene vector, or a
#'   per-array vector.
#' @param method `"ls"`. `"robust"` is refused, see refusals above.
#' @param ... Passed to `gls.series`, and from there to `duplicateCorrelation`.
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
#' @param chunks Row chunks. Defaults to `workers`; passing `chunks = workers` explicitly is redundant.
#'   Clamped so no chunk holds one gene.
#' @param parallel_backend One of [combat_backends()], or a function
#'   `function(idx, f, workers)` returning a list in the order of `idx`. Defaults to
#'   `getOption("combat.backend", combat_default_backend())`.
#' @param backend Optional `lmFit` to wrap. Defaults to `limma::lmFit`.
#'
#' @param label Optional name for this call in the timing line, when
#'   `options(combat.timing = TRUE)` is set. Defaults to the companion and the matrix shape,
#'   e.g. `ComBat-seq 18,270 x 1,500`; pass a cohort name to tell calls apart in a loop.
#' @return An `MArrayLM`, `identical()` to what the backend's own `lmFit` returns for the
#'   same input. Under either size gate, on a branch a block would flip, or at
#'   `workers = 1`, that is literally what it is: one plain call to the backend with a
#'   dispatch wrapper around it and no parallelism added.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("limma", quietly = TRUE)) {
#'   set.seed(1)
#'   y <- matrix(rnorm(4000), nrow = 500)
#'   design <- cbind(Intercept = 1, Group = rep(0:1, each = 4))
#'
#'   fit <- lmFit_parallel(y, design, workers = 4L)
#'   identical(fit, limma::lmFit(y, design))
#' }
#' }
#' @export
lmFit_parallel <- function(object, design = NULL, ndups = NULL, spacing = NULL,
                           block = NULL, correlation, weights = NULL, method = "ls", ...,
                           workers = NULL, chunks = NULL,
                           parallel_backend = getOption("combat.backend", combat_default_backend()),
                           backend = NULL, label = NULL) {
  # no run leaves workers behind, crashed or not; children that predate this call are
  # someone else's and are spared
  .spare <- combat_children()
  on.exit(combat_reap(.spare), add = TRUE)

  workers <- rp_prologue(workers)

  # timing and quieting are on.exit hooks, so an error unwinds the sink and still reports the
  # elapsed line: a failed run says where it failed instead of vanishing. Placed after the
  # prologue because that is what resolves `workers` from NULL to a number worth printing.
  .rp <- rp_step_begin(label, "lmFit", object, parallel_backend, workers)
  on.exit(rp_step_end(.rp), add = TRUE)
  if (!is.function(parallel_backend)) {
    parallel_backend <- match.arg(parallel_backend, combat_backends())
  }

  method <- match.arg(method, c("ls", "robust"))
  if (method == "robust") {
    stop("method = \"robust\" fits through mrlm, which rnaparallel does not split. Use ",
         "limma::lmFit directly for it.", call. = FALSE)
  }
  if (isTRUE(ndups > 1)) {
    stop("ndups >= 2 makes unwrapdups reshape the matrix, so a gene no longer occupies ",
         "one row and a row split would cut genes in half. Use limma::lmFit directly.",
         call. = FALSE)
  }

  be <- limma_backend(backend,
                      need_args = c("object", "design", "ndups", "spacing", "block",
                                    "correlation", "weights", "method"),
                      rebound = c("lm.series", "gls.series"))
  vendor_lm <- get("lm.series", envir = be$env, inherits = TRUE)
  vendor_gls <- get("gls.series", envir = be$env, inherits = TRUE)

  # ndups is also reachable from y$printer, which lmFit resolves after the check above has
  # already run, so the refusal is repeated where the resolved value arrives.
  refuse_ndups <- function(ndups) {
    if (isTRUE(ndups > 1)) {
      stop("ndups = ", ndups, " came from the data object's printer layout. unwrapdups ",
           "reshapes rows, so rnaparallel cannot split it. Use limma::lmFit directly.",
           call. = FALSE)
    }
  }

  env <- new.env(parent = be$env)

  env$lm.series <- function(M, design = NULL, ndups = 1, spacing = 1, weights = NULL) {
    refuse_ndups(ndups)
    # Leaned for the same reason rp_row_blocks leans its own dispatch closure: this frame's
    # parent is the entry point's, which holds `object`, so a block closure defined here drags
    # a second full copy of the input onto every socket worker.
    lean <- new.env(parent = globalenv())
    lean$vendor_lm <- vendor_lm; lean$design <- design
    lean$ndups <- ndups; lean$spacing <- spacing
    blk <- function(Mi, wi) vendor_lm(Mi, design = design, ndups = ndups, spacing = spacing,
                                      weights = wi)
    environment(blk) <- lean
    rp_row_blocks(
      M, weights, be$env, workers, chunks, parallel_backend, "lmFit_parallel/lm.series",
      blk,
      function() vendor_lm(M, design = design, ndups = ndups, spacing = spacing,
                           weights = weights))
  }

  env$gls.series <- function(M, design = NULL, ndups = 2, spacing = 1, block = NULL,
                             correlation = NULL, weights = NULL, ...) {
    refuse_ndups(ndups)
    # Resolved on the FULL matrix, with the raw weights, which is where the vendor resolves
    # it. Left to the blocks it would be a trimmed mean over each block's genes alone.
    if (is.null(correlation)) {
      dupcor <- get("duplicateCorrelation", envir = be$env, inherits = TRUE)
      correlation <- dupcor(M, design = design, ndups = ndups, spacing = spacing,
                            block = block, weights = weights, ...)$consensus.correlation
    }
    # Same lean rebuild as lm.series. `...` cannot live in a detached environment, so the
    # dots are captured as a list and spliced back in the same position they occupied, which
    # leaves argument matching unchanged.
    dots <- list(...)
    lean <- new.env(parent = globalenv())
    lean$vendor_gls <- vendor_gls; lean$design <- design; lean$ndups <- ndups
    lean$spacing <- spacing; lean$block <- block; lean$correlation <- correlation
    lean$dots <- dots
    blk <- function(Mi, wi) do.call(vendor_gls,
      c(list(Mi, design = design, ndups = ndups, spacing = spacing, block = block,
             correlation = correlation, weights = wi), dots))
    environment(blk) <- lean
    rp_row_blocks(
      M, weights, be$env, workers, chunks, parallel_backend, "lmFit_parallel/gls.series",
      blk,
      function() vendor_gls(M, design = design, ndups = ndups, spacing = spacing,
                            block = block, correlation = correlation, weights = weights, ...),
      punch = "gls")
  }

  f <- be$fn
  environment(f) <- env
  args <- list(object = object, design = design, ndups = ndups, spacing = spacing,
               block = block, weights = weights, method = method, ...)
  # single-bracket with a list(), or an explicit correlation = NULL would add nothing and
  # arrive as missing, which is a different branch of lmFit
  if (!missing(correlation)) args["correlation"] <- list(correlation)
  do.call(f, args, quote = TRUE)
}
