## helper_seq_parallel.R
##
## Mirrors upstream helper_seq.R. Everything ComBat_seq_parallel needs that is not
## the entry point itself. Sourcing this file plus ComBat_seq_parallel.R is enough
## to use the function without installing anything, exactly like upstream.
##
## Nothing here reimplements ComBat-seq. Three hot paths are split by row, the common
## dispersion is dispatched across batches, and each
## slice is handed to the original ComBat-seq function.


# ---- backend resolution ------------------------------------------------------

#' Heads of calls whose head is a bare symbol
#'
#' `all.names()` flattens `edgeR::glmFit` into parts that read as reachable, which is exactly
#' the degradation the rebind gates exist to catch, so only bare heads count. `e[[i]]` is
#' indexed inline: binding it first makes R treat an empty symbol as a missing argument.
#' @noRd
rp_bare_call_heads <- function(e) unique(rp_call_heads_raw(e))

#' The walk itself, without the dedup
#'
#' `unique()` used to run at every node. It keeps first occurrences in input order, and a
#' value's globally-first occurrence is also the first inside its own subtree, so no inner
#' pass could ever delete one: the result is identical, order included, with the dedup done
#' once at the top. Measured on the bodies this package actually walks -- limma::lmFit 483 to
#' 229 us, duplicateCorrelation 612 to 333 us, sva::ComBat_seq 926 to 482 us. That is fixed
#' overhead on every companion call, and most of the wrapper's cost on a small input.
#' @noRd
rp_call_heads_raw <- function(e) {
  if (!is.call(e)) return(character())
  h <- e[[1L]]
  out <- if (is.name(h)) as.character(h) else character()
  for (i in seq_along(e)) {
    if (is.call(e[[i]])) out <- c(out, rp_call_heads_raw(e[[i]]))
  }
  out
}

#' Resolve a ComBat-seq backend and its helper
#'
#' Finds the ComBat-seq function to run and the `match_quantiles` helper it
#' needs, taking both from the same environment so the pair can never be drawn
#' from two different copies of ComBat-seq.
#'
#' Two layouts work with no configuration. Bioconductor `sva` keeps `ComBat_seq`
#' and `match_quantiles` in its namespace. The upstream repository defines both at
#' top level once `ComBat_seq.R` and `helper_seq.R` are sourced. In both cases
#' `environment(fn)` already reaches the helpers, so resolving through it avoids
#' naming a namespace and avoids `:::` altogether.
#'
#' @param fn A ComBat-seq function. Defaults to `sva::ComBat_seq` when `sva` is
#'   installed, otherwise a top-level `ComBat_seq` if one is visible.
#' @return A list with `fn`, `env` and `match_quantiles`.
#' @noRd
combat_backend <- function(fn = NULL) {
  if (is.null(fn)) {
    if (requireNamespace("sva", quietly = TRUE)) {
      fn <- sva::ComBat_seq
    } else if (exists("ComBat_seq", mode = "function")) {
      fn <- get("ComBat_seq", mode = "function")
    } else {
      stop("no ComBat-seq backend found. Install sva (BiocManager::install(\"sva\")) ",
           "or source ComBat_seq.R and helper_seq.R from ",
           "https://github.com/zhangyuqing/ComBat-seq", call. = FALSE)
    }
  }
  if (!is.function(fn)) stop("`fn` must be a function", call. = FALSE)

  env <- environment(fn)
  if (is.null(env)) {
    stop("the ComBat-seq backend has no environment, so its helpers cannot be ",
         "resolved. A primitive or a function stripped by compilation cannot be ",
         "wrapped this way.", call. = FALSE)
  }

  # gate: the backend must take the arguments we forward, and match_quantiles must
  # take the five we split by row. A silent upstream rename or signature change
  # would otherwise surface as wrong numbers rather than an error.
  need <- c("counts", "batch", "group", "covar_mod", "full_mod",
            "shrink", "shrink.disp", "gene.subset.n")
  missing_args <- setdiff(need, names(formals(fn)))
  if (length(missing_args)) {
    stop("this ComBat-seq backend is missing argument(s): ",
         paste(missing_args, collapse = ", "),
         ". rnaparallel was written against the 8-argument signature shared ",
         "by sva 3.54.0 and upstream master.", call. = FALSE)
  }

  if (!exists("match_quantiles", envir = env, inherits = TRUE)) {
    stop("could not find `match_quantiles` alongside the ComBat-seq backend. ",
         "With sva, that means the internal helper was renamed upstream; with a ",
         "sourced copy, that helper_seq.R was not sourced.", call. = FALSE)
  }
  mq <- get("match_quantiles", envir = env, inherits = TRUE)
  mq_need <- c("counts_sub", "old_mu", "old_phi", "new_mu", "new_phi")
  if (!identical(names(formals(mq)), mq_need)) {
    stop("`match_quantiles` has signature (", paste(names(formals(mq)), collapse = ", "),
         ") but rnaparallel splits it by row on (", paste(mq_need, collapse = ", "),
         "). Refusing to run rather than pass arguments positionally into a changed ",
         "function.", call. = FALSE)
  }

  # A rebind only works while the backend calls these as bare symbols. If upstream ever
  # namespace-qualifies one, or swaps sapply for vapply, the rebind becomes unreachable and
  # this package quietly degrades to a pass-through: output stays identical(), every
  # equivalence test still passes, and nothing is parallelised. Fail loudly instead.
  reachable <- rp_bare_call_heads(body(fn))
  rebound <- c("glmFit", "glmFit.default", "match_quantiles",
               "estimateGLMTagwiseDisp", "sapply", "lapply")
  unreachable <- setdiff(rebound, reachable)
  if (length(unreachable)) {
    stop("this ComBat-seq backend no longer calls ", paste(unreachable, collapse = ", "),
         " as a bare symbol, so rebinding cannot reach it. The package would run serially ",
         "while still returning identical() output, which no equivalence test can detect. ",
         "Refusing to run.", call. = FALSE)
  }

  list(fn = fn, env = env, match_quantiles = mq)
}


# ---- row chunking ------------------------------------------------------------

#' Split row indices into interleaved chunks
#'
#' Chunks are INTERLEAVED, not contiguous blocks, because gene order is not random.
#' A count matrix that has been filtered or sorted puts the cheap genes together and
#' the expensive ones together, and `qnbinom` cost rises with count size. Measured on a
#' matrix sorted by expression, contiguous blocks carried 0.36, 2, 17 and 112 million
#' counts: a 313x imbalance, so three workers idled while the fourth did the work, and
#' the speedup fell from 1.66x to 1.47x. Round-robin makes the same split 1.01x even
#' and is immune to whatever order the caller's genes arrive in.
#'
#' The cost is that rows come back permuted, so every caller must restore the original
#' order after binding. `order(unlist(idx))` does that, and is a no-op for a contiguous
#' split, which is why the callers apply it unconditionally.
#'
#' @param ntag Number of rows (genes).
#' @param workers Worker count, used when `chunks` is NULL.
#' @param chunks Chunk count. Clamped to `[1, ntag]`, so it cannot exceed one chunk per gene.
#'   Concurrency is bounded separately by `workers`, so more chunks than workers means more
#'   forks at the same concurrency, not more parallelism.
#' @param interleave Round-robin the rows across chunks. `FALSE` gives the old
#'   contiguous blocks, which is only useful for reproducing the imbalance.
#' @param min_rows Smallest number of rows any chunk may hold. The chunk count is
#'   clamped to `ntag %/% min_rows` so no chunk falls below it. Callers whose original
#'   function branches on block shape pass 2; the default 1 reproduces the old clamp.
#' @return A list of integer vectors covering `seq_len(ntag)` exactly once.
#' @noRd
combat_row_chunks <- function(ntag, workers = 4L, chunks = NULL, interleave = TRUE,
                              min_rows = 1L) {
  ntag <- as.integer(ntag)
  if (ntag < 1L) stop("no rows to split", call. = FALSE)
  # An unusable `chunks` clamps to one chunk rather than erroring: results stay correct and
  # only parallelism is lost. Asserted for 0, -5 and NA in test-parallel.R.
  nch <- suppressWarnings(as.integer(if (is.null(chunks)) workers else chunks))
  if (is.na(nch)) nch <- 1L
  min_rows <- max(1L, suppressWarnings(as.integer(min_rows)))
  if (is.na(min_rows)) min_rows <- 1L
  # min_rows is an exactness constraint, not a tuning knob. limma's lm.series reaches
  # stats::lm.fit, whose `if (is.matrix(y) && ny == 1L) y <- drop(y)` demotes a one-gene
  # block to a vector, flipping sigma from colMeans to mean and dropping the gene name.
  # Measured: 4080 of 16000 one-gene blocks not identical() to serial. At min_rows = 1
  # this reduces to the old clamp exactly, so ComBat-seq callers are unaffected.
  nch <- max(1L, min(nch, ntag %/% min_rows))
  if (nch == 1L) return(list(seq_len(ntag)))
  # `rep_len(seq_len(nch), ntag)` puts chunk k at positions k, k + nch, k + 2*nch, ..., and
  # split() returns those groups in sorted key order, so group k IS seq.int(k, ntag, by = nch).
  # split.default gets there by coercing an ntag-long vector to a factor, which sorts its
  # uniques and builds character levels: measured 0.312 ms against 0.0067 ms at ntag = 18,270.
  # Same integer vectors, same unnamed list, so every downstream row mapping is untouched.
  if (interleave) return(lapply(seq_len(nch), function(k) seq.int(k, ntag, by = nch)))
  unname(split(seq_len(ntag), cut(seq_len(ntag), nch, labels = FALSE)))
}

#' Restore original row order after chunked results are bound together
#'
#' A no-op when the split was contiguous, since `unlist(idx)` is then already sorted.
#' @param idx The chunk list the work was dispatched on.
#' @return Integer permutation to apply to the bound rows.
#' @noRd
combat_row_order <- function(idx) {
  # The inverse permutation, built by scatter rather than by sorting. `unlist(idx)` is a
  # permutation of seq_len(n), so its inverse is what `order()` returns, and computing it
  # directly is O(n) against O(n log n): about 1 ms a bind at 18,000 genes, paid once per
  # dispatch on every companion.
  u <- unlist(idx, use.names = FALSE)
  ord <- integer(length(u))
  ord[u] <- seq_along(u)
  ord
}


#' Where this package's own functions live
#'
#' The parent every lean environment gets. A dispatched closure is rebuilt against an
#' environment holding just what its body reads, and that environment's PARENT decides how the
#' body's remaining calls resolve -- `vapply`, `unname`, `do.call`, `[`. Parented at
#' `globalenv()` they resolved through the user's workspace first, so a user binding named
#' `vapply` changed the companion's numbers while leaving the original untouched: measured, a
#' shadowed `vapply` moved `calcNormFactors_parallel` 8.59e-06 off `edgeR::normLibSizes` with
#' no error and no warning, on the serial backend, on every platform.
#'
#' This resolves to the namespace when the package is installed and to whatever the files were
#' sourced into otherwise, which the package supports and which is exactly where those calls
#' resolved before the closures were leaned. Serialisation is unaffected: a namespace parent is
#' written as a reference, the same as `globalenv()`.
#' @noRd
rp_home <- function() environment(rp_copy_free)


#' Does a dispatch on this backend hand a worker the payload without copying it?
#'
#' The three size gates below were tuned where a worker INHERITS the matrix: a forked child
#' starts almost free and reads the parent's pages copy-on-write, so a split earns its keep as
#' soon as the arithmetic is big enough. Where every chunk is serialised into a worker instead,
#' the cheap-per-cell paths never repay the transfer.
#'
#' The gates used to ask `.Platform$OS.type == "windows"`, which is the wrong question twice
#' over. Windows is a SUFFICIENT condition for copying, not a necessary one: a `multisession`
#' plan is socket-based everywhere, and `options(combat.fork = FALSE)` turns the lot serial by
#' request. And the question is not fork() at all. `foreach` here builds a FORK cluster on
#' Unix, so it does fork, and is still slow, because doParallel's cluster form serialises every
#' TASK over a socket whatever its nodes were made with. Measured on 300,000 x 24, four
#' workers, every arm `identical()`:
#'
#'   limma::lmFit                              0.150 s
#'   mclapply                                  0.111 s   1.35x
#'   foreach, this package's own pool          0.628 s   0.24x
#'   foreach, registerDoParallel(cores = 4)    0.135 s   1.11x
#'
#' So the predicate is about the COPY, not the fork, and `foreach` has to be asked rather than
#' assumed: doParallel reports `doParallelMC` when it is driving `mclapply`, which copies
#' nothing, and `doParallelSNOW` when it is driving a cluster, which copies every task. An
#' unregistered `foreach` gets this package's own cluster and is therefore the copying form.
#'
#' `serial` answers FALSE deliberately. Nothing dispatches there, so the gate decides between
#' one whole original call and the original walked over blocks in one process, and the whole call
#' is the faster of the two: the fast `lm.series` branch measured 0.70x split that way.
#' @param parallel_backend The resolved backend, a name or a function.
#' @return TRUE when a worker reads the payload without a serialised copy.
#' @noRd
rp_copy_free <- function(parallel_backend) {
  if (identical(.Platform$OS.type, "windows")) return(FALSE)
  if (!isTRUE(getOption("combat.fork", TRUE))) return(FALSE)
  # A custom executor keeps the inherited answer where fork() exists. It could be either -- the
  # documented examples span parLapply over a FORK cluster and furrr over a socket plan -- and
  # the two ways to be wrong are not symmetric. Guessing "copies" would silently switch off a
  # split a caller had been getting, which is the class of change this package refuses to make
  # behind someone's back; guessing "inherits" leaves them where they are, and the option is
  # there for anyone who measures otherwise.
  if (is.function(parallel_backend)) return(TRUE)
  switch(as.character(parallel_backend)[1L],
    mclapply = TRUE,
    BiocParallel = TRUE,               # MulticoreParam forks wherever fork() exists
    future = isTRUE(tryCatch(
      future::supportsMulticore() && inherits(future::plan(), "multicore"),
      error = function(e) FALSE)),
    foreach = {
      nm <- tryCatch(if (foreach::getDoParRegistered()) foreach::getDoParName() else NULL,
                     error = function(e) NULL)
      isTRUE(nm %in% c("doParallelMC", "doMC"))
    },
    FALSE)                             # serial, and anything unrecognised
}


#' Size gate for a row or column split, resolved against the backend that will run it
#'
#' 6e6 cells is the least-squares break-even where the payload is INHERITED, and 2e4 the
#' weighted branch's. Measured on
#' Windows over PSOCK with the gate forced open, 1,200 samples, the split never approaches
#' parity and gets worse with every worker added:
#'
#'   18,000 genes  (21.6M cells, original 0.5 s)  0.14x at 2 workers, 0.08x at 4
#'   50,000 genes  (60.0M cells, original 1.4 s)  0.24x at 2 workers, 0.10x at 4
#'
#' On the TCGA cohort the same split measured 0.50x. There is no threshold that rescues that,
#' so without fork the gate CLOSES rather than rising and the fast branch stays whole.
#'
#' Each branch consults its OWN option. That is not a detail: the merged Windows branch routed
#' both branches through one function whose first act was to return `combat.min.ls.cells` when
#' it was set, so a caller who raised the least-squares gate silently raised the voom/weighted
#' one from 2e4 to the same value and switched off a split the docs measure at 2.52x-3.39x.
#' The suite could not see it, because setup-parallel.R sets every gate to 0 and that function
#' returns from its first line for both branches throughout.
#' @param option Name of the option this branch honours.
#' @param fork_default Threshold where a dispatch forks.
#' @param parallel_backend The resolved backend.
#' @return Cell count below which the split does not run.
#' @noRd
rp_ls_min_cells <- function(option, fork_default, parallel_backend) {
  o <- getOption(option)
  if (!is.null(o)) return(o)
  if (rp_copy_free(parallel_backend)) fork_default else Inf
}


#' Gene floor for limma's per-gene least-squares loop
#'
#' The weighted branch of `lm.series`/`gls.series` is an interpreted loop over genes, and its
#' per-gene cost has a large component that does not scale with array count. So what amortises
#' a fork here is GENES per worker, not cells, and a cell gate answers the wrong question: at a
#' fixed 1,000 genes the split behaves the same at 8, 24 and 48 arrays while the cell count
#' moves 8,000 to 48,000, straddling the 2e4 cell gate that used to decide it.
#'
#' Measured on this M3 at the default worker count, gate forced open, median of 9, companion
#' against original, every arm `identical()`, at 8 / 24 / 48 arrays:
#'
#'   500 genes    0.71x  0.57x  0.75x
#'   1,000        1.16x  1.26x  1.25x
#'   2,000        1.39x  1.48x  1.75x
#'   4,000        1.95x  2.02x  2.16x
#'
#' An independent review measured the same shape on the same machine and put 1,000 genes at
#' 0.78x / 0.83x / 0.97x, a loss where this run has a win. The two disagree only there, and
#' both agree from 2,000 up, so the default sits where both agree rather than at either run's
#' own crossover. Being too high forgoes about 1.2x on a thin matrix; being too low returns
#' 0.57x, and a companion slower than the function it wraps is the thing these gates exist to
#' prevent.
#' @return Gene count below which the weighted least-squares split does not run.
#' @noRd
rp_wt_min_genes <- function() {
  o <- getOption("combat.min.wt.genes")
  if (!is.null(o)) return(o)
  2000L
}


#' Size gate for the TMM column split
#'
#' Same shape of problem as [rp_ls_min_cells()] and a different answer, which is why it is
#' measured rather than assumed. 2e5 cells is the fork break-even; over sockets the column loop
#' has to earn a serialised copy per chunk as well. Measured on Windows with the gate forced
#' open:
#'
#'   1.8M cells  (original 1.2 s)   1.05x at 2 workers, 0.89x at 4, 0.64x at 6
#'   21.6M cells (original 15.8 s)  1.58x at 2 workers, 1.53x at 4, 1.36x at 6
#'
#' Unlike the least-squares split this one does pay at scale, so the gate rises instead of
#' closing: break-even sits around two million cells, an order of magnitude above the fork
#' value. On the cohort the companion measured 1.63x with the gate in place.
#' @return Cell count below which the column split does not run.
#' @noRd
rp_norm_min_cells <- function(parallel_backend) {
  o <- getOption("combat.min.norm.cells")
  if (!is.null(o)) return(o)
  if (rp_copy_free(parallel_backend)) 2e5 else 2e6
}


#' Size gate for the order-statistic column loops, RLE and upperquartile
#'
#' These take one order statistic per column and cost roughly a fifteenth as much per cell as
#' the TMM loop, so 4e6 is their fork break-even where TMM's is 2e5. The same reasoning that
#' raises TMM's gate tenfold without fork applies at least as strongly to a loop that cheap,
#' so this rises by the same factor rather than being left at the fork value.
#' @return Cell count below which the column split does not run.
#' @noRd
rp_order_min_cells <- function(parallel_backend) {
  o <- getOption("combat.min.order.cells")
  if (!is.null(o)) return(o)
  if (rp_copy_free(parallel_backend)) 4e6 else 4e7
}


#' The backend to use when the caller has not named one
#'
#' `mclapply` everywhere fork() exists, which is every platform this package was built on. On
#' Windows it cannot fork, so it is correct, serial, and says so once -- a companion that returns
#' the right answer at the speed of the function it wraps.
#'
#' `future` is the only backend that runs real workers there without the caller registering a
#' cluster first, and on the Windows machine that comparison was made on it measured 2.92x on the
#' cohort against `foreach`'s 0.28x. But it needs a
#' plan, and this package will not set one: a caller's plan is theirs, and silently starting
#' sixteen processes inside someone's session is worse than being slow. Defaulting to `future`
#' unconditionally would therefore hand most Windows users a warning on every dispatch and no
#' speedup, which is strictly worse than what they have now.
#'
#' So it is chosen adaptively: taken when a plan is already active and it will actually do
#' something, and left alone otherwise. A user who sets `plan(multisession)` gets parallelism
#' without also having to discover `options(combat.backend=)`; a user who sets nothing gets
#' exactly today's behaviour. Resolved per call, since a plan can be set at any time.
#' @return A backend name from [combat_backends()].
#' @noRd
combat_default_backend <- function() {
  if (!identical(.Platform$OS.type, "windows")) return("mclapply")
  if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("future.apply", quietly = TRUE)) return("mclapply")
  # nbrOfWorkers() is the question that matters: a sequential plan, an unset plan and a
  # one-worker plan all answer 1, and all three mean `future` would resolve in this process.
  n <- tryCatch(future::nbrOfWorkers(), error = function(e) 1L)
  if (isTRUE(n > 1L)) "future" else "mclapply"
}




# ---- memory guard -------------------------------------------------------------

#' This process's own parent PID, or NA when nothing on this build can read it
#'
#' `Sys.getppid()` is not a base R function on every build: it does not exist at all on the
#' R 4.6.1 UCRT Windows build this package is verified against (`exists("Sys.getppid")` is
#' FALSE there, not merely unimplemented for one caller), so calling it unconditionally
#' crashed every dispatch on that build rather than skipping cleanly. Falls back to
#' `ps::ps_ppid()` when that package is installed and `Sys.getppid()` is not present; NA
#' otherwise, and NA is treated by every caller as "cannot tell, skip the check", the same
#' convention as the rest of the memory-guard readers in this file.
#' @noRd
rp_getppid <- function() {
  if (exists("Sys.getppid", where = baseenv(), mode = "function")) {
    return(tryCatch(get("Sys.getppid", envir = baseenv())(), error = function(e) NA_integer_))
  }
  if (requireNamespace("ps", quietly = TRUE)) {
    return(tryCatch(as.integer(ps::ps_ppid()), error = function(e) NA_integer_))
  }
  NA_integer_
}

#' Bytes of RAM the kernel thinks are actually obtainable right now
#'
#' MemAvailable, not MemFree: free excludes reclaimable page cache and reads far
#' lower than what a fork can really have. NA on anything without /proc, which is
#' every non-Linux platform, and every caller treats NA as "cannot tell, proceed".
#' @noRd
rp_mem_available <- function() {
  if (!file.exists("/proc/meminfo")) return(NA_real_)
  l <- tryCatch(readLines("/proc/meminfo", n = 64L), error = function(e) character())
  m <- grep("^MemAvailable:", l, value = TRUE)
  if (!length(m)) return(NA_real_)
  as.numeric(gsub("\\D", "", m[1L])) * 1024
}

#' Total installed bytes of RAM, not what is currently free
#'
#' MemTotal from /proc/meminfo on Linux; `Get-CimInstance Win32_ComputerSystem` via
#' PowerShell on Windows, since R_MAX_VSIZE is meaningful there too and
#' rnaparallel_set_mem_limit() needs a number to halve on any platform a caller runs on.
#' `wmic` was tried first and dropped: it no longer exists on current Windows builds
#' (removed from Windows 11 24H2 onward), so it returned "command not found" rather than a
#' number and rnaparallel_set_mem_limit() read that as "cannot tell" on every affected
#' machine. NA when neither source is readable, which callers treat as "ask the user".
#' @noRd
rp_mem_total <- function() {
  if (file.exists("/proc/meminfo")) {
    l <- tryCatch(readLines("/proc/meminfo", n = 64L), error = function(e) character())
    m <- grep("^MemTotal:", l, value = TRUE)
    if (length(m)) return(as.numeric(gsub("\\D", "", m[1L])) * 1024)
    return(NA_real_)
  }
  if (identical(.Platform$OS.type, "windows")) {
    out <- tryCatch(
      system2("powershell", c("-NoProfile", "-Command",
              "(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory"),
              stdout = TRUE, stderr = FALSE),
      error = function(e) character())
    val <- suppressWarnings(as.numeric(trimws(out)))
    val <- val[!is.na(val) & val > 0]
    if (length(val)) return(val[1L])
  }
  NA_real_
}


#' Resident bytes of this process
#'
#' Field 24 of /proc/self/stat is RSS in pages. Read from stat rather than status
#' because it is one line and one split, and this runs on every dispatch.
#' @noRd
rp_mem_rss <- function() {
  if (!file.exists("/proc/self/stat")) return(NA_real_)
  v <- tryCatch(strsplit(readLines("/proc/self/stat", n = 1L), " ", fixed = TRUE)[[1L]],
                error = function(e) character())
  if (length(v) < 24L) return(NA_real_)
  suppressWarnings(as.numeric(v[24L])) * 4096
}

#' Cap workers at what the machine can actually hold
#'
#' Forking is copy-on-write, so N workers do not cost N times the parent. They cost
#' whatever each one writes to, and for a row-split fit over a large matrix that is a
#' real fraction of it. When the total exceeds what is available the kernel does not
#' return an allocation error to R: it SIGKILLs the process. There is no condition to
#' catch, no traceback, and mclapply reports nothing, so the run simply disappears.
#' Measured on a 40,609 x 9,493 matrix with 125 GB and no swap: 16 workers off a 50 GB
#' parent died, 8 off 100 GB died, 4 off 23 GB died at 111 GB, 2 off 23 GB survived.
#'
#' So refuse at the door instead. Degrading to fewer workers finishes late; being
#' killed loses the whole run and says nothing about why.
#'
#' combat.mem.divergence is the fraction of the parent each worker is assumed to
#' dirty. It is workload-dependent, not a constant, which is why it is an option: a
#' row-split GLM fit dirties far more than a per-column trimmed mean.
#' @noRd
rp_mem_cap <- function(workers) {
  if (!isTRUE(rp_opt_flag("combat.mem.guard", default = TRUE))) return(workers)
  if (workers <= 1L) return(workers)
  avail <- rp_mem_available(); rss <- rp_mem_rss()
  if (is.na(avail) || is.na(rss) || rss <= 0) return(workers)   # cannot tell, do not interfere
  frac <- rp_opt_num("combat.mem.divergence", 0.25)
  if (frac <= 0) return(workers)
  headroom <- avail * 0.8            # leave a fifth for everything that is not this fit
  need <- rss * frac * workers
  if (need <= headroom) return(workers)
  fit <- max(1L, as.integer(floor(headroom / (rss * frac))))
  if (fit >= workers) return(workers)
  warning(sprintf(
    paste0("rnaparallel: %d workers need ~%.0f GB on top of a %.0f GB parent and only ",
           "%.0f GB is available, which on a machine without swap is a kernel kill, not ",
           "an R error. Using %d instead. Set options(combat.mem.divergence=) if this ",
           "workload dirties less, or options(combat.mem.guard=FALSE) to disable."),
    workers, need / 2^30, rss / 2^30, avail / 2^30, fit), call. = FALSE)
  fit
}


# ---- entry-point prologue -----------------------------------------------------

#' Validate the four shared controls once, identically, for every entry point
#'
#' The copies this replaces had already drifted into two textually different variants.
#' Returns the validated worker count; resolves the backend name in the caller's frame.
#' @noRd
rp_prologue <- function(workers) {
  workers <- combat_default_workers(workers)
  if (length(workers) != 1L) stop("`workers` must be a single integer", call. = FALSE)
  if (is.numeric(workers) && is.finite(workers) && workers != trunc(workers)) {
    stop("`workers` must be a whole number, not ", workers, call. = FALSE)
  }
  w <- suppressWarnings(as.integer(workers))
  if (is.na(w) || w < 1L) stop("`workers` must be a positive integer", call. = FALSE)
  # a garbage combat.min.* value should refuse at the door, not one dispatch later
  for (op in c("combat.min.cells", "combat.min.disp.cells", "combat.min.ls.cells",
               "combat.min.norm.cells", "combat.min.order.cells", "combat.min.glm.cells",
               "combat.min.dupcor.cells", "combat.min.batch.cells", "combat.min.wt.genes")) {
    v <- getOption(op)
    if (!is.null(v)) {
      vv <- suppressWarnings(as.numeric(v))
      if (length(vv) != 1L || is.na(vv) || vv < 0) {
        stop("`", op, "` must be a single non-negative number; got ", deparse(v),
             call. = FALSE)
      }
    }
  }
  rp_mem_cap(w)
}


# ---- worker cleanup -----------------------------------------------------------

#' Reap the fork children an entry point created, killing any that are stuck
#'
#' `mclapply` kills its own children when a call ends cleanly, but a child wedged in a
#' signal-unsafe state survives that, and survivors accumulate across the many dispatches
#' one analysis makes. A session that ended holding dozens of them has crashed this
#' machine. Every entry point snapshots the session's children on entry and reaps only
#' what appeared since, so a user's own concurrent workers (a `future` plan, their own
#' `mclapply`) are never touched. `getFromNamespace` rather than `:::`, since `children`
#' is the one handle base R gives to a fork child that no longer responds.
#' @noRd
combat_children <- function() {
  if (.Platform$OS.type == "windows") return(integer())
  ch <- tryCatch(utils::getFromNamespace("children", "parallel")(), error = function(e) NULL)
  vapply(ch, function(p) as.integer(p$pid), integer(1))
}

#' @section What this does NOT cover:
#' Nothing on the default `mclapply` path, and that is not a defect here. `parallel::children()`
#' maps to `mc_children`, which counts only entries that are not DETACHED, and `mclapply`
#' detaches every child inside its own `on.exit(cleanup(mc.cleanup))` -- which runs on error and
#' on interrupt alike, before an entry point's `on.exit(combat_reap(.spare))` is ever reached.
#' Measured: after a normal dispatch, after an error thrown inside a worker, and after an error
#' unwinding through `mclapply` in the caller, `children()` is 0 every time. So on the default
#' backend this function correctly finds nothing, and in a healthy session that is because there
#' is nothing: no live child, no zombie.
#'
#' What it does cover is an ATTACHED child, which is a `future` multicore worker that has not
#' resolved, or a bare `mcparallel()`. Both are reaped. Worth stating because the rail in
#' test-safety-rails.R exercises the `mcparallel` shape, which no companion creates through
#' `mclapply`, so that test passing is not evidence about the default path.
#' @noRd
combat_reap <- function(spare = integer()) {
  if (.Platform$OS.type == "windows") return(invisible(0L))
  children <- utils::getFromNamespace("children", "parallel")
  ours <- function() {
    ch <- tryCatch(children(), error = function(e) NULL)
    Filter(function(p) !(as.integer(p$pid) %in% spare), ch)
  }
  ch <- ours()
  n <- length(ch)
  if (n) {
    try(suppressWarnings(parallel::mccollect(ch, wait = FALSE)), silent = TRUE)
    ch <- ours()
    for (p in ch) try(tools::pskill(p$pid, tools::SIGKILL), silent = TRUE)

    deadline <- proc.time()[["elapsed"]] + 2
    repeat {
      ch <- ours()
      if (!length(ch)) break
      try(suppressWarnings(parallel::mccollect(ch, wait = FALSE)), silent = TRUE)
      ch <- ours()
      if (!length(ch) || proc.time()[["elapsed"]] >= deadline) break
      Sys.sleep(0.01)
    }
    if (length(ch)) {
      warning("could not collect killed worker process(es): ",
              paste(vapply(ch, function(p) p$pid, integer(1)), collapse = ", "),
              call. = FALSE)
    }
  }
  invisible(n)
}


# ---- cluster reuse ------------------------------------------------------------

# Clusters are cached per (type, size) and reused. Building one is expensive and one
# ComBat-seq run dispatches twice for the GLM fits, once per batch for the quantile match,
# once more per batch for tagwise dispersion where that batch has the residual degrees of
# freedom for it, and once across batches for the common dispersion. So 3 + n_batch +
# eligible, up to 2 * n_batch + 3: nine on a 3-batch design and 203 on a 100-batch one. Constructing per dispatch paid the cost every one of
# those times: measured 153 ms for a 4-worker PSOCK cluster, which made socket backends
# look 25x slower than serial when the frameworks themselves were fine.
.combat_clusters <- new.env(parent = emptyenv())

#' Read process identities that survive cleanup retries
#'
#' One `ps` for the whole pool, not one per PID. Probing each worker separately forked an
#' external process per worker per dispatch: measured 9.1 ms per `ps`, 32 ms to answer for a
#' 4-worker pool, against 0.02 ms for the signal-0 liveness test alone. A PID the kernel no
#' longer knows is simply absent from the output, so absence answers liveness too.
#'
#' @param pids Integer PIDs.
#' @return Character vector along `pids`, `NA` where the process is gone or unprobeable.
#' @noRd
combat_pid_identities <- function(pids) {
  out <- rep(NA_character_, length(pids))
  if (!length(pids) || .Platform$OS.type == "windows") return(out)
  got <- tryCatch(
    suppressWarnings(system2(
      "ps", c("-o", "pid=,ppid=,lstart=,command=", "-p", paste(pids, collapse = ",")),
      stdout = TRUE, stderr = FALSE, env = "LC_ALL=C"
    )),
    error = function(e) character()
  )
  got <- trimws(got[nzchar(trimws(got))])
  if (!length(got)) return(out)
  seen <- suppressWarnings(as.integer(sub("^\\s*([0-9]+).*$", "\\1", got)))
  hit  <- match(pids, seen)
  out[!is.na(hit)] <- got[hit[!is.na(hit)]]
  out
}

#' Classify recorded PIDs without treating a reused PID as an owned worker
#'
#' `gone` the process no longer exists; `owned` it is still the worker we recorded;
#' `foreign` the PID has been recycled by something that is not ours, so it must never be
#' signalled; `unknown` it is alive but its identity could not be read. Unknown is treated as
#' ours everywhere downstream: a PID this package recorded is its responsibility until it is
#' shown to belong to somebody else, and refusing to signal it is what made a stuck entry
#' immortal.
#' @noRd
combat_pid_status <- function(pids, identities = NULL) {
  if (!length(pids)) return(character())
  now <- combat_pid_identities(pids)
  alive <- vapply(pids, function(pid) {
    tryCatch(isTRUE(tools::pskill(pid, 0L)), error = function(e) NA)
  }, logical(1))
  vapply(seq_along(pids), function(i) {
    if (identical(alive[i], FALSE) && is.na(now[i])) return("gone")
    if (is.null(identities) || length(identities) < i ||
        is.na(identities[i]) || is.na(now[i])) return("unknown")
    if (identical(now[i], identities[i])) "owned" else "foreign"
  }, character(1))
}

#' Test whether worker processes still exist
#'
#' Liveness only, and deliberately without the identity probe: on a cache hit the question is
#' "are my workers there", and a recycled PID that slipped past this is caught immediately
#' after by the content-checking `clusterCall` probe. Identity is a kill-time guard, paid for
#' only where this package actually signals a process.
#' @noRd
combat_pids_alive <- function(pids) {
  if (!length(pids)) return(logical())
  if (.Platform$OS.type == "windows") return(rep(NA, length(pids)))
  vapply(pids, function(pid) {
    tryCatch(isTRUE(tools::pskill(pid, 0L)), error = function(e) FALSE)
  }, logical(1))
}

#' Wait a bounded time for worker processes to exit
#' @noRd
combat_wait_pids <- function(pids, identities = NULL, timeout = 2) {
  if (!length(pids) || .Platform$OS.type == "windows") return(TRUE)
  deadline <- proc.time()[["elapsed"]] + timeout
  repeat {
    pending <- if (is.null(identities)) {
      combat_pids_alive(pids)
    } else {
      combat_pid_status(pids, identities) %in% c("owned", "unknown")
    }
    if (!any(pending) || proc.time()[["elapsed"]] >= deadline) return(!any(pending))
    Sys.sleep(0.01)
  }
}

#' Resolve the default worker count for a companion entry point
#'
#' `workers = NULL` means "pick one", and the pick is `min(8, detectCores() - 2)`. The upper
#' bound is 8 because that is where the measured curves flatten and beyond which added
#' workers return almost nothing. The `- 2` leaves the machine a core for the master and one
#' for everything else, which is what keeps a small machine from being pushed past what it
#' has: 8 cores gives 6, 4 cores gives 2, a 64-core node still gives 8.
#'
#' Explicit `workers` is never touched, and the ceiling is not a safety guarantee. Forking
#' six workers alongside a SECOND forking R session kernel-panicked a 24 GB machine, and no
#' core count this process can read knows about the other session, so a caller running two
#' at once should pass a smaller number.
#' Without fork() the pick is additionally capped at the PERFORMANCE core count, which is not
#' the core count on a hybrid CPU. The distinction only matters where there is no fork: a forked
#' worker on an efficiency core still adds throughput, which is why the macOS run measures 5.37x
#' at eight workers on a chip with four performance cores and is left alone here. Over sockets
#' every worker also costs a serialised copy, so a slow core stops paying for itself: on an Ultra
#' 185H, 6 performance plus 10 efficiency, the ComBat-seq curve peaks at 6 workers and turns over
#' after. Reading the core count alone would have picked 8.
#' @noRd
combat_default_workers <- function(workers = NULL) {
  if (!is.null(workers)) return(workers)
  ac <- .combat_clusters$allcores
  if (is.null(ac)) {
    ac <- max(1L, suppressWarnings(parallel::detectCores()), na.rm = TRUE)
    .combat_clusters$allcores <- ac
  }
  n <- max(1L, min(8L, ac - 2L))
  if (identical(.Platform$OS.type, "windows")) n <- min(n, rp_perf_cores())
  max(1L, n)
}

#' Performance cores, distinguished from total cores on a hybrid CPU
#'
#' There used to be two of these and they disagreed. This one read the SMT topology and
#' answered 8 on an Apple M3; a second copy inside `combat_parallel_lapply` shelled out to
#' `sysctl hw.perflevel0.physicalcpu` and answered 4, which is right. One routine now, asking
#' each platform the best question it can answer, memoised once.
#'
#' Darwin reports performance cores outright, so it is asked first. Everything else falls back
#' to the topology, because `detectCores(logical = FALSE)` reports 16 on an Ultra 185H that has
#' 6 performance cores and 10 efficiency ones -- not a rounding error, but the difference
#' between a sensible default and one that recruits ten slow workers.
#'
#' The topology answers it without an original table: on Intel hybrid parts only performance cores
#' carry SMT, so the number of logical processors ABOVE the physical count is the number of cores
#' with a second thread. 22 logical against 16 physical gives 6, which is right. The formula
#' degrades correctly everywhere else: on a uniformly hyperthreaded machine every core has a
#' second thread, `logical - physical == physical`, and the guard below returns the physical
#' count; with SMT off or absent the difference is 0 and it does the same. A cgroup-restricted
#' container can report fewer logical than physical, which is also caught.
#' @return Best available performance-core count, never below 1.
#' @noRd
rp_perf_cores <- function() {
  p <- .combat_clusters$perfcores
  if (!is.null(p)) return(p)
  # Darwin reports performance cores outright; nothing else here does.
  p <- suppressWarnings(as.integer(
    if (identical(Sys.info()[["sysname"]], "Darwin"))
      tryCatch(system2("sysctl", c("-n", "hw.perflevel0.physicalcpu"),
                       stdout = TRUE, stderr = FALSE), error = function(e) NA)
    else NA))
  if (length(p) != 1L || is.na(p) || p < 1L) {
    phys <- suppressWarnings(parallel::detectCores(logical = FALSE))
    logi <- suppressWarnings(parallel::detectCores(logical = TRUE))
    if (is.na(phys) || phys < 1L) phys <- if (is.na(logi) || logi < 1L) 1L else logi
    smt <- if (is.na(logi)) 0L else logi - phys
    p <- if (smt > 0L && smt < phys) smt else phys   # hybrid: only P-cores carry SMT
  }
  p <- max(1L, as.integer(p))
  .combat_clusters$perfcores <- p
  p
}

# A retired entry is retried on later calls, so it needs a stop condition. Without one, a
# worker wedged in an uninterruptible syscall -- the case this whole layer exists for -- made
# every later dispatch pay the full timeout, forever, with nothing said.
.combat_retire_tries <- 5L

#' Preserve worker PIDs for signal-only cleanup
#'
#' Stamped with the owning process, exactly as a cached cluster handle is. A forked child
#' inherits this list, and without the stamp it would SIGKILL the workers its parent is still
#' using -- the same hazard the cluster cache guards against by refusing to touch a pool it
#' did not create.
#' @noRd
combat_retire <- function(entry) {
  if (.Platform$OS.type == "windows" || !length(entry$worker_pids)) return(FALSE)
  ids <- entry$worker_identities
  if (length(ids) != length(entry$worker_pids)) ids <- rep(NA_character_, length(entry$worker_pids))
  retired <- .combat_clusters$retired
  retired[[length(retired) + 1L]] <- list(
    pid = Sys.getpid(), pids = entry$worker_pids, identities = ids, tries = 0L
  )
  .combat_clusters$retired <- retired
  TRUE
}

#' Retry retired workers by PID and identity, never through a closed connection
#'
#' @return Number of retired entries cleared. Not a cluster count; see [combat_cluster_stop()].
#' @noRd
combat_retired_reap <- function(timeout = 2) {
  retired <- .combat_clusters$retired
  if (!length(retired)) return(invisible(0L))
  keep <- list()
  reaped <- 0L
  for (entry in retired) {
    # a pool recorded by another process is not ours to signal; forget it without touching it
    if (!identical(entry$pid, Sys.getpid())) next
    status <- combat_pid_status(entry$pids, entry$identities)
    pending <- status %in% c("owned", "unknown")
    if (!any(pending)) {
      reaped <- reaped + 1L
      next
    }
    entry$pids <- entry$pids[pending]
    entry$identities <- entry$identities[pending]
    for (pid in entry$pids) try(tools::pskill(pid, tools::SIGKILL), silent = TRUE)
    if (combat_wait_pids(entry$pids, entry$identities, timeout)) {
      reaped <- reaped + 1L
      next
    }
    entry$tries <- (if (is.null(entry$tries)) 0L else entry$tries) + 1L
    if (entry$tries >= .combat_retire_tries) {
      warning("gave up killing worker process(es) after ", entry$tries, " attempts: ",
              paste(entry$pids, collapse = ", "), call. = FALSE)
      next
    }
    keep[[length(keep) + 1L]] <- entry
  }
  .combat_clusters$retired <- if (length(keep)) keep else NULL
  invisible(reaped)
}

#' Stop an owned cluster without writing after a node is suspected dead
#' @noRd
combat_cluster_teardown <- function(cl, suspected_dead = FALSE) {
  if (is.null(cl)) return(invisible(FALSE))
  close_node <- utils::getFromNamespace("closeNode", "parallel")
  stopped <- FALSE

  if (!suspected_dead) {
    stopped <- tryCatch({
      suppressWarnings(parallel::stopCluster(cl))
      TRUE
    }, error = function(e) FALSE)
  }
  if (!stopped) for (node in cl) try(suppressWarnings(close_node(node)), silent = TRUE)
  # The isOpen sweep is the only evidence-based half of the answer, so it is the whole answer.
  # Gating on `stopped` too returned FALSE for every clean shutdown that went down the
  # fallback path on macOS and Linux, which is every suspected-dead teardown there is.
  invisible(all(vapply(cl, function(node) {
    !tryCatch(isOpen(node$con), error = function(e) FALSE)
  }, logical(1))))
}

#' Get a cached parallel cluster, creating it on first use
#'
#' @param ncore Worker count.
#' @param type Cluster type, `"FORK"` or `"PSOCK"`.
#' @return A cluster object owned by this package. Do not stop it directly; use
#'   [combat_cluster_stop()].
#' @noRd
combat_cluster <- function(ncore, type = if (.Platform$OS.type == "windows") "PSOCK" else "FORK",
                           packages = "edgeR") {
  # One pool per TYPE, not per (type, size). Keying on size too meant trying 2, 4, 6 and 8
  # workers left four pools and 20 worker processes alive until an explicit stop, on a
  # machine whose documented worker ceiling exists because forking too wide panics it.
  key <- type
  entry <- .combat_clusters[[key]]

  # A pool inherited by a forked child is not ours to use or to stop: the child would race
  # the parent on the same sockets, and stopping it would kill the parent's workers. Forget
  # the handle without touching the processes behind it.
  if (!is.null(entry) && !identical(entry$pid, Sys.getpid())) {
    rm(list = key, envir = .combat_clusters)
    entry <- NULL
  }
  # After the ownership check, never before it: the retired list is inherited across a fork
  # exactly as the handle is, and reaping first let a child kill workers its parent still had.
  combat_retired_reap()
  cl <- entry$cl
  # Only when the pool is too SMALL. Rebuilding on any width change is what made foreach
  # collapse: ComBat-seq alternates between batch-level and row-level dispatches and so asks
  # for different widths within one run, and each change tore the pool down and built another.
  # Measured on PSOCK, 20 calls alternating 3/6 workers: 6299 ms against 22 ms at a fixed
  # width, 315 ms per call against 1.1 ms. A wider pool serves a narrower request by handing
  # back the first `ncore` nodes; the surplus sits idle, concurrency is unchanged, and chunk
  # placement is decided by the tag rather than by which node ran what.
  if (!is.null(cl) && length(cl) < ncore) {
    tracked <- combat_retire(entry)
    closed <- combat_cluster_teardown(cl)
    rm(list = key, envir = .combat_clusters)
    if (tracked) combat_retired_reap()
    if (!tracked && .Platform$OS.type == "windows" && !closed) {
      warning("could not close every node in the cached cluster", call. = FALSE)
    }
    entry <- NULL
    cl <- NULL
  }
  # A cached cluster can be dead if the session forked, or the user stopped it by hand,
  # or DESYNCHRONIZED if an interrupt left an unread result sitting in a worker socket.
  # The desynchronized case is the dangerous one: the probe would consume the stale
  # message, conclude the cluster was healthy, and every later dispatch would return the
  # PREVIOUS dispatch's chunk with no error. So check what came back, not just that
  # something did.
  open <- !is.null(cl) && all(vapply(cl, function(node) {
    tryCatch(isOpen(node$con), error = function(e) FALSE)
  }, logical(1)))
  processes_alive <- !is.null(cl) && (
    .Platform$OS.type == "windows" ||
      (length(entry$worker_pids) == length(cl) &&
         all(combat_pids_alive(entry$worker_pids)))
  )
  alive <- open && processes_alive && tryCatch({
    res <- parallel::clusterCall(cl, function() TRUE)
    length(res) == length(cl) && all(vapply(res, isTRUE, logical(1)))
  }, error = function(e) FALSE)
  if (!alive) {
    if (!is.null(cl)) {
      tracked <- combat_retire(entry)
      closed <- combat_cluster_teardown(cl, suspected_dead = TRUE)
      rm(list = key, envir = .combat_clusters)
      if (tracked) combat_retired_reap()
      if (!tracked && .Platform$OS.type == "windows" && !closed) {
        warning("could not close every node in the dead cached cluster", call. = FALSE)
      }
    }
    # Forked siblings inherit the same precomputed port and collide; retry on a socket
    # error with a pid-jittered port rather than dying in one sibling.
    cl <- tryCatch(parallel::makeCluster(ncore, type = type), error = function(e) e)
    tries <- 0L
    while (inherits(cl, "error") && grepl("cannot be opened", conditionMessage(cl)) &&
           tries < 3L) {
      tries <- tries + 1L
      port <- 11000L + (Sys.getpid() * 7L + tries * 131L) %% 20000L
      cl <- tryCatch(parallel::makeCluster(ncore, type = type, port = port),
                     error = function(e) e)
    }
    if (inherits(cl, "error")) stop(cl)
    worker_pids <- tryCatch(
      unlist(parallel::clusterCall(cl, Sys.getpid), use.names = FALSE),
      error = function(e) e
    )
    if (inherits(worker_pids, "error") || length(worker_pids) != length(cl)) {
      combat_cluster_teardown(cl, suspected_dead = TRUE)
      if (inherits(worker_pids, "error")) stop(worker_pids)
      stop("new cluster did not report every worker PID", call. = FALSE)
    }
    # Windows PSOCK workers are still owned by this cluster, but this package has no
    # portable start-time identity probe there.  Keep the placeholder identities for
    # cache shape; Windows never uses PID-only retirement/reaping.
    #
    # A missing identity is NOT fatal. `ps` failing to fork under exactly the memory pressure
    # this package is built for would otherwise abort a caller's analysis on a cluster that
    # had just been built successfully. An unidentified PID is recorded as ours, which
    # combat_pid_status reads as "unknown" and still signals; the cost is that this one
    # worker cannot be told apart from a recycled PID later, not that the run dies now.
    worker_identities <- if (.Platform$OS.type == "windows") {
      rep(NA_character_, length(worker_pids))
    } else {
      combat_pid_identities(worker_pids)
    }
    assign(key, list(cl = cl, pid = Sys.getpid(), worker_pids = worker_pids,
                     worker_identities = worker_identities),
           envir = .combat_clusters)
  }
  # PSOCK workers start empty, so the packages the closures call must be loaded there.
  # FORK workers inherit them. Done on every call rather than only at creation: one pool is
  # cached per TYPE and shared by callers needing different packages, so a cluster built for
  # an edgeR dispatch would otherwise reach a limma closure with limma absent.
  # Send only what this pool has not already loaded. A worker cannot unload a namespace, and
  # every path that replaces the workers writes a fresh cache entry with no record, so the memo
  # resets exactly when they do. Measured 2.78 ms per dispatch of pure socket round trips, paid
  # up to 2 * n_batch + 3 times a run.
  if (type == "PSOCK" && length(packages)) {
    entry <- .combat_clusters[[key]]
    todo <- setdiff(packages, entry$packages_loaded)
    if (length(todo)) {
      for (pkg in todo) {
        invisible(parallel::clusterCall(
          cl, function(p) suppressMessages(requireNamespace(p, quietly = TRUE)), pkg))
      }
      entry$packages_loaded <- c(entry$packages_loaded, todo)
      assign(key, entry, envir = .combat_clusters)
    }
  }
  # A pool wider than asked for serves the request from its first `ncore` nodes.
  if (length(cl) > ncore) cl[seq_len(ncore)] else cl
}

#' Stop and forget every cluster this package cached
#'
#' The `foreach` backend reuses a cached cluster rather than building one per dispatch.
#' Call this when done with a long parallel session, or if a cluster is
#' misbehaving; the next call rebuilds it. Clusters are also stopped when the
#' package namespace unloads.
#'
#' @details
#' Building a cluster is expensive and one ComBat-seq run dispatches 3 + n_batch +
#' eligible tagwise batches, up to 2 * n_batch + 3, so nine on a 3-batch design.
#' Measured on a 4-worker PSOCK cluster: 153 ms to build, and five dispatches
#' went from 631 ms to 5 ms once the cluster was reused.
#'
#' @return Number of clusters stopped, invisibly.
#' @export
combat_cluster_stop <- function() {
  # the cache holds plain flags beside the cluster handles; sweeping them would delete
  # registered_by_us before the foreach reset below reads it, and re-arm the one-time messages
  # perfcores belongs here as much as allcores does. Left off, the `!is.list(entry)` sweep
  # below deleted the memo on every stop and the next call paid another detectCores() pair.
  keys <- setdiff(ls(.combat_clusters),
                  c("warned_windows", "warned_ecores", "registered_by_us", "perf",
                    "perfcores", "allcores", "bpparam", "retired"))
  combat_retired_reap()                 # for its effect; it counts entries, not clusters
  stopped <- 0L
  for (k in keys) {
    entry <- .combat_clusters[[k]]
    # the cache also holds plain flags, and `TRUE$pid` is an error rather than NULL, so a
    # non-list entry is dropped before anything reads a field off it
    if (!is.list(entry)) { rm(list = k, envir = .combat_clusters); next }
    # only stop pools this process created; a handle inherited through a fork belongs to
    # the parent and stopping it would kill workers the parent is still using
    if (!identical(entry$pid, Sys.getpid())) { rm(list = k, envir = .combat_clusters); next }
    tracked <- combat_retire(entry)
    closed <- combat_cluster_teardown(entry$cl)
    # forgetting a handle whose teardown failed, with no PIDs recorded either, orphans the
    # workers with nothing left that can reach them. Windows never retires, so that is the
    # platform this protects; keep the handle and let the next call try again.
    if (!tracked && !closed) {
      warning("could not close every node in the cached cluster; it will be retried",
              call. = FALSE)
      next
    }
    rm(list = k, envir = .combat_clusters)
    if (tracked) combat_retired_reap()
    stopped <- stopped + 1L
  }
  # Only reset foreach if WE registered it. Doing it unconditionally replaced the backend
  # of a caller who had registered their own and never used the foreach path here at all.
  if (isTRUE(.combat_clusters$registered_by_us) &&
      requireNamespace("foreach", quietly = TRUE)) {
    try(foreach::registerDoSEQ(), silent = TRUE)
    .combat_clusters$registered_by_us <- NULL
  }
  invisible(stopped)
}

.onUnload <- function(libpath) combat_cluster_stop()


# ---- the one dispatch point --------------------------------------------------

#' Backends this package can dispatch to
#'
#' All parallelism in this package funnels through one internal dispatch point, so
#' supporting another framework is a branch there and nothing else.
#'
#' @details
#' Every backend returns bit-identical results, because each row chunk is a pure
#' function of its own genes and every backend listed preserves chunk order.
#'
#' \describe{
#'   \item{`mclapply`}{Default. Forks via `parallel::mclapply`. Unix only; falls
#'     back to serial on Windows.}
#'   \item{`future`}{`future.apply::future_lapply`. The caller owns the plan, this
#'     package will not set one. Warns if the plan resolves in one process.}
#'   \item{`BiocParallel`}{`BiocParallel::bplapply` with `MulticoreParam`. No new
#'     dependency in practice, since \pkg{sva} depends on it.}
#'   \item{`foreach`}{`foreach::%dopar%` over a cached cluster from `doParallel`,
#'     FORK on Unix and PSOCK on Windows. Slower than the original on the limma and
#'     edgeR paths, measured 0.24x to 0.41x against `limma::lmFit` at four workers,
#'     because `doParallel` serialises each task's closure and the closure captures
#'     the matrix. Correct, and worth choosing only where a fork is unavailable.}
#'   \item{`serial`}{Plain `lapply`. Same output, no workers.}
#' }
#'
#' @return Character vector of accepted `parallel_backend` values.
#' @examples
#' combat_backends()
#' @export
combat_backends <- function() {
  c("mclapply", "future", "BiocParallel", "foreach", "serial")
}

#' Apply a function over row chunks, on a chosen backend
#'
#' Every backend must satisfy three properties or the identical() promise breaks:
#' results come back in the order the chunks went in, each chunk is evaluated
#' exactly once, and no backend rewrites the returned values. Order is the one
#' most easily lost, which is why nothing here uses an unordered map.
#'
#' The chunk bodies contain no RNG, so the choice of backend cannot shift results
#' through the random stream. ComBat-seq's own `sample()` call lives in
#' `monte_carlo_int_NB`, on the serial side of the companion.
#'
#' @param idx List of index vectors from `combat_row_chunks()`.
#' @param f Function applied to one index vector.
#' @param workers Maximum concurrent workers.
#' @param parallel_backend One of [combat_backends()], or a function
#'   `function(idx, f, workers)` returning a list in the order of `idx`. The
#'   function form is the extension point: any framework can be plugged in without
#'   this package growing a branch for it. `options(combat.fork = FALSE)` still
#'   forces serial, so a custom executor cannot defeat the escape hatch.
#' @param cells Size of the work being dispatched, in matrix cells. Below
#'   `getOption("combat.min.cells", 20000)` the dispatch runs serially, because
#'   the fork costs more than the work it saves. `Inf`, the default, means a
#'   caller that has not measured its own work size always dispatches.
#' @return A list, one element per chunk, in order.
#' @noRd
combat_parallel_lapply <- function(idx, f, workers,
                                   parallel_backend = getOption("combat.backend", combat_default_backend()),
                                   cells = Inf,
                                   min_cells = getOption("combat.min.cells", 2e4),
                                   preschedule = FALSE) {
  # `f` is evaluated on every path this function has, but not until a worker touches it, and
  # until then it is a promise. serialize() writes a promise together with its PRENV, so each
  # dispatched task carried the caller's evaluation frame, that frame's own unforced promises,
  # and through them the original's frames and the entry point's raw inputs. Measured on the
  # upperquartile dispatch: 14,387,651 B per task before, 2,026,039 B after. `idx` needs no
  # such treatment; the tagging below forces it before any branch.
  force(f)
  workers <- as.integer(workers)
  if (is.na(workers)) stop("`workers` must be a positive integer", call. = FALSE)

  # Validate BEFORE the size gate, or validation becomes size-dependent: a misspelled
  # backend name would be refused on a big matrix and silently accepted on a small one.
  custom <- is.function(parallel_backend)
  if (!custom) parallel_backend <- match.arg(parallel_backend, combat_backends())

  # Same reasoning for the packages a named backend needs. Checked below the gate, an
  # uninstalled framework succeeded quietly on a small dispatch and errored on an
  # otherwise identical large one, which is the worst possible time to find out.
  if (!custom) {
    needs <- switch(parallel_backend, future = "future.apply", BiocParallel = "BiocParallel",
                    foreach = c("foreach", "doParallel"), character(0))
    absent <- needs[!vapply(needs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(absent)) {
      stop("parallel_backend = \"", parallel_backend, "\" needs ",
           paste(absent, collapse = " and "), call. = FALSE)
    }
  }

  # ComBat-seq calls the hot paths once per BATCH, so a 100-batch run dispatches
  # thin slices rather than one fat one, and each fork has to earn its cost against
  # a fraction of a percent of the matrix. Measured on 10-column slices: 0.50x at
  # 5,000 cells, 0.92x at 10,000, 1.39x at 20,000. Below the threshold this was not
  # merely slower than it could be, it was slower than not forking at all:
  # 500 genes x 1000 samples x 100 batches ran at 0.80x against plain sva.
  #
  # The threshold is per path, not global, because the paths cost different amounts per
  # cell. `qnbinom` in match_quantiles is dear enough to pay for a fork at 20,000 cells;
  # dispersion estimation is cheaper per cell and does not break even until about 30,000,
  # so one shared threshold made the dispersion split a net loss on small matrices with
  # many batches. See `combat.min.disp.cells`.
  #
  # as.numeric because a character option would make this a lexicographic comparison,
  # which silently switches the gate off for exactly the sizes it exists to catch.
  # Refuse a bad threshold rather than ignore it. Coercing and then skipping on NA meant a
  # typo like "2000O" silently switched the safety gate off, which is the opposite of what
  # someone setting the option wants.
  fk <- getOption("combat.fork", TRUE)
  if (!(is.logical(fk) && length(fk) == 1L && !is.na(fk))) {
    stop("`combat.fork` must be TRUE or FALSE; got ", deparse(fk),
         ". A garbage value used to force every dispatch serial with no signal.",
         call. = FALSE)
  }

  mc <- suppressWarnings(as.numeric(min_cells))
  if (length(mc) != 1L || is.na(mc) || mc < 0) {
    stop("`combat.min.cells` and `combat.min.disp.cells` must each be a single ",
         "non-negative number; got ", deparse(min_cells), call. = FALSE)
  }
  # counted, not just taken: a call that ran entirely under the gates returns identical()
  # output at serial pace and is indistinguishable from one that forked, unless something says so
  if (isTRUE(cells < mc)) { rp_count(FALSE); return(lapply(idx, f)) }

  # Every job carries its chunk number, and the number comes back with the result. Checking
  # only the LENGTH of what a backend returns cannot tell a correct answer from the same
  # chunks in the wrong order, and binding a reordered list scrambles genes silently rather
  # than failing. With the tag, a permutation is detected and put back; a duplicate or a
  # missing chunk is refused.
  tag <- function(k) structure(idx[[k]], combat_chunk = k)
  idx_tagged <- lapply(seq_along(idx), tag)
  # f_tagged wraps every job, so its own frame travels on every task of every backend. Left on
  # this frame it carried `idx`, the full duplicate `idx_tagged`, `tag` and `untag`: measured
  # 22,295 B per task at 6,000 rows, about 150 KB at 18,270 genes, times chunks, times up to
  # 2 * n_batch + 3 dispatches a run.
  .lean_tag <- new.env(parent = parent.env(environment()))
  .lean_tag$f <- f
  # Captured as plain values in the MASTER, not read from .rp_dispatch inside the worker: a
  # PSOCK worker is a fresh process that does not inherit the master's global env or its
  # options() unless told to, and rp_progress_dir()/the stage label would silently read as
  # unset there. mclapply's forked children do inherit a copy, but capturing here keeps one
  # code path for all four backends instead of a fork-only shortcut that quietly breaks PSOCK.
  .lean_tag$master_pid     <- Sys.getpid()
  .lean_tag$progress_dir   <- rp_progress_dir()
  .lean_tag$progress_stage <- .rp_dispatch$progress_label %||% "dispatch"
  f_tagged <- function(ii) {
    k <- attr(ii, "combat_chunk")
    attr(ii, "combat_chunk") <- NULL
    # Mark THIS process as running a dispatch, so a rebound hot path that dispatches again
    # from inside the worker runs serially instead of opening a second pool. The flag is set
    # here rather than in the master because the worker is the only process the nested call
    # can consult: a PSOCK worker is a fresh R process that inherits nothing set after it
    # started, and a fork child gets its own copy of the environment.
    prev <- Sys.getenv("RNAPARALLEL_IN_WORKER", unset = NA_character_)
    Sys.setenv(RNAPARALLEL_IN_WORKER = "1")
    on.exit(if (is.na(prev)) Sys.unsetenv("RNAPARALLEL_IN_WORKER")
            else Sys.setenv(RNAPARALLEL_IN_WORKER = prev), add = TRUE)
    # File progress covers the stretch the console tick cannot: once the master calls into
    # mclapply/future/BiocParallel/foreach it blocks until every chunk returns, so this is
    # the only place inside that block anything can be reported from. No-op when
    # combat.progress.dir was unset in the master (rp_progress_file_write() checks `dir`).
    # Bare symbols, not `.lean_tag$progress_dir`: `environment(f_tagged) <- .lean_tag` makes
    # THIS function's own scope .lean_tag, so a name bound in it resolves directly, the same
    # way `f` above is called bare rather than as `.lean_tag$f`. `.lean_tag$progress_dir`
    # would look up a variable named `.lean_tag` INSIDE .lean_tag itself, which does not
    # exist, and future's globals scan does not pull the enclosing binding in either --
    # measured as `object '.lean_tag' not found` under a real multisession plan.
    if (!is.null(progress_dir)) {
      rp_progress_file_write(progress_dir, progress_stage, k, "start")
      on.exit(rp_progress_file_write(progress_dir, progress_stage, k, "done"), add = TRUE)
    }
    # A fork whose master was killed is reparented to init and keeps running, holding its
    # share of the matrix for as long as the machine is up. Nothing collects its result and
    # nothing reaps it, so the memory the master died for stays gone: measured at 111 GB and
    # 116 GB held by such orphans on two separate runs. SIGTERM does not clear them either,
    # because R installs a handler and the worker is blocked, so only SIGKILL works and only
    # if somebody notices. getppid() == 1 is the one signal a fork can read for itself.
    #
    # Sys.getppid() is not a base R function on every build: it does not exist at all on the
    # R 4.6.1 UCRT Windows build this was verified against (`exists("Sys.getppid")` is FALSE,
    # not merely unavailable to a worker), and the original code called it unconditionally,
    # which crashed every single dispatch under a real future::multisession run with "could
    # not find function Sys.getppid" -- not a worker-only failure, the whole mechanism assumed
    # a base function that is not universally present. rp_getppid() below fails to NA rather
    # than erroring, and NA skips the check exactly like the platforms that never had
    # Sys.getppid to begin with: on Windows there is no fork at all for this to protect
    # against, so skipping is correct there regardless.
    #
    # Guarded on the pid differing from the master's: f_tagged also runs IN the master under a
    # serial or custom backend, and quitting there would take the caller's whole session down.
    # A master launched with setsid (any nohup/detached render) legitimately has ppid 1.
    if (Sys.getpid() != master_pid && isTRUE(rp_getppid() == 1L))
      quit(save = "no", status = 0L, runLast = FALSE)
    list(combat_chunk = k, value = f(ii))
  }
  environment(f_tagged) <- .lean_tag

  # Put results back in dispatch order and strip the wrapper. Anything that is not a tagged
  # result (a try-error, a NULL from a killed worker, a condition from foreach) is passed
  # through untouched at its own position for combat_parallel_check() to report on.
  untag <- function(out) {
    ids <- vapply(out, function(o) {
      if (is.list(o) && !is.null(o$combat_chunk)) as.integer(o$combat_chunk) else NA_integer_
    }, integer(1))
    ok <- !is.na(ids)
    if (any(ok)) {
      if (anyDuplicated(ids[ok]) || any(ids[ok] < 1L) || any(ids[ok] > length(idx))) {
        stop("the parallel backend returned duplicated or out-of-range chunks, so the ",
             "result cannot be reassembled. Genes would be scrambled rather than an error ",
             "raised, which is why this is checked.", call. = FALSE)
      }
    }
    res <- vector("list", length(idx))
    filled <- logical(length(idx))
    # Single-bracket assignment with a list() wrapper, NEVER res[[i]] <- value. A killed
    # mclapply child returns a plain NULL, and `x[[i]] <- NULL` DELETES the element and
    # shrinks the list rather than storing it. That silently dropped the dead chunk, shifted
    # every later chunk up one slot, and scrambled genes: measured 100 genes NA and 200
    # carrying another gene's dispersion, with combat_parallel_check reporting success.
    for (j in seq_along(out)) {
      if (ok[j]) { res[ids[j]] <- list(out[[j]]$value); filled[ids[j]] <- TRUE }
    }
    # Untagged elements are only ever error placeholders. A backend that strips the tag
    # wrapper and returns a real value cannot be placed: with equal-sized chunks it would pass
    # every later check while genes sat in another chunk's rows.
    for (j in which(!ok)) {
      o <- out[[j]]
      if (!(is.null(o) || inherits(o, "try-error") || inherits(o, "condition"))) {
        stop("the parallel backend returned a value without the chunk tag this package ",
             "attached. Its position cannot be recovered, and placing it by order would bind ",
             "genes into the wrong rows rather than raise anything. Refusing.", call. = FALSE)
      }
    }
    spare <- which(!filled)
    for (j in which(!ok)) if (length(spare)) { res[spare[1]] <- list(out[[j]]); spare <- spare[-1] }
    res
  }

  # A dispatch already running inside one of this package's workers must not open a second
  # pool. ComBat-seq dispatches the tagwise loop ACROSS BATCHES and ships the original closure,
  # whose environment still carries the rebound `estimateGLMTagwiseDisp`; inside the worker
  # that symbol dispatches AGAIN over gene rows. The result is workers + workers^2 processes:
  # measured 2 outer and 4 inner for `workers = 2L` on Windows, which extrapolates to 272 at
  # the 16-worker arm, each one a fresh R process with edgeR and limma loaded.
  #
  # This used to be handled by `mc.allow.recursive = FALSE`, but that is an argument to
  # parallel::mclapply and guards the fork branch alone. On Windows mclapply is serial and
  # foreach/PSOCK is the only backend that runs workers at all, so the one platform that
  # needed the guard was the one platform without it. The flag is process-local and set by
  # `f_tagged` above, so this covers every backend, custom executors included.
  #
  # A caller's OWN parallel loop is unaffected: nothing marks their workers, so a companion
  # called once per cohort inside their loop still parallelises, which is the nesting pattern
  # the documentation recommends.
  if (nzchar(Sys.getenv("RNAPARALLEL_IN_WORKER"))) {
    rp_count(FALSE)
    return(lapply(idx, f))
  }

  if (custom) {
    if (!fk) { rp_count(FALSE); return(lapply(idx, f)) }
    rp_count(TRUE)
    out <- parallel_backend(idx_tagged, f_tagged, workers)
    if (!is.list(out) || length(out) != length(idx)) {
      stop("a custom parallel_backend must return a list of length ", length(idx),
           ", in the order of `idx`; got ", class(out)[1], " of length ", length(out),
           call. = FALSE)
    }
    return(untag(out))
  }

  # options(combat.fork = FALSE) forces serial regardless of backend: slower,
  # identical output, no worker processes. The escape hatch when an IDE wedges,
  # since forking inside RStudio is not officially supported.
  forced_serial <- !fk
  degenerate <- workers <= 1L || length(idx) <= 1L
  if (parallel_backend == "serial" || forced_serial || degenerate) {
    rp_count(FALSE)
    return(lapply(idx, f))
  }
  rp_count(TRUE)

  # Clamped to real cores as well as chunks: `workers` is a user number and nothing
  # else bounded it, so workers = 64 on an 8-core box forked 64 processes at once.
  #
  # R CMD check sets _R_CHECK_LIMIT_CORES_ and then errors on more than two, which is
  # why this cap has to exist rather than trusting detectCores(). It only started
  # mattering once the test suite began actually forking; before that no test reached
  # a worker at all, so check passed while proving nothing.
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  # Past the performance-core count each added worker returns less throughput than the one
  # before it, so more workers can be slower rather than faster. Not a core-type effect: forks
  # migrate across the performance and efficiency clusters, and eight children given identical
  # work finish within 1.08x of each other. That is reported once and
  # left to the caller: `workers` is an explicit request, and cutting it behind the caller's
  # back makes the argument mean something other than what it says. Capping it silently also
  # made `workers = 8` run 8 chunks at 4-way concurrency, which reads as a worker-count result
  # when it is really a chunking result.
  # Memoised: both branches spawn a process on macOS, about 10 ms each, and one ComBat-seq
  # run dispatches up to 2 * n_batch + 3 times. The core count does not change mid-session.
  perf <- rp_perf_cores()

  # Memoised for the same reason as perf: detectCores() spawns a process on macOS, measured
  # 13 ms, and this ran on every forking dispatch on every backend.
  cores_cap <- if (nzchar(chk) && !identical(tolower(chk), "false")) 2L else {
    ac <- .combat_clusters$allcores
    if (is.null(ac)) {
      ac <- max(1L, suppressWarnings(parallel::detectCores()), na.rm = TRUE)
      .combat_clusters$allcores <- ac
    }
    ac
  }
  ncore <- min(workers, length(idx), cores_cap)

  # Only where the dispatch does NOT fork. A forked worker on an efficiency core still adds
  # throughput, because it shares the parent's pages rather than being handed a copy: the
  # macOS run measures 5.37x at eight workers on a chip with four performance cores. Ungated,
  # this fired on a stock M3 at the package's OWN default -- min(8, 8 - 2) = 6 against 4
  # performance cores -- so the package warned that its default might be slower than a number
  # it had declined to pick, and contradicted its own published measurement. Over sockets the
  # warning is real, because there each added worker also costs a serialised copy.
  if (!rp_copy_free(parallel_backend) && ncore > perf &&
      is.null(.combat_clusters$warned_ecores)) {
    .combat_clusters$warned_ecores <- TRUE
    message("workers = ", workers, " exceeds the ", perf, " performance core(s) this machine ",
            "reports. Without fork() each added worker also costs a serialised copy, so past ",
            "that point this can be slower than workers = ", perf, ".")
  }

  switch(parallel_backend,
    mclapply = {
      # fork only. Windows has no fork, so fall back rather than error.
      if (.Platform$OS.type == "windows") {
    # a silent serial run looks identical to a parallel one until you time it
    if (is.null(.combat_clusters$warned_windows)) {
      # Not foreach. Measured on Windows, foreach fell 1.18x, 0.94x, 0.57x, 0.28x at 2, 4, 8
      # and 16 workers while future held 2.92x on the cohort of the day, and combat_default_backend()
      # picks future itself once a plan exists. Sending people to the slower one at the exact
      # moment they discover the problem is the opposite of helping.
      message("mclapply cannot fork on Windows, so this ran serially. ",
              "Set future::plan(future::multisession, workers = N) and this will use the ",
              "future backend automatically.")
      .combat_clusters$warned_windows <- TRUE
    }
    rp_count_serial_after_all()   # it did not fork; the line must not claim it did
    return(lapply(idx, f))
  }
      # mc.allow.recursive = FALSE, or a caller who wraps this in their own
      # mclapply/future_lapply over cohorts multiplies the worker count instead of
      # reusing it: 3 cohorts x 4 workers measured 12 concurrent grandchildren.
      # R degrades to lapply inside an already-forked child, which is what we want.
      untag(parallel::mclapply(idx_tagged, f_tagged, mc.cores = ncore,
                               mc.preschedule = preschedule, mc.allow.recursive = FALSE))
    },

    future = {
      # The caller owns the plan. Setting one here would stamp on a plan the user
      # established for the whole session, which is the usual complaint about
      # packages that touch future's global state.
      # Test behaviour, not class. Under RStudio `supportsMulticore()` is FALSE and
      # `plan(multicore)` falls back per future rather than rewriting the plan object,
      # so its class never says "sequential" and the old check never fired: measured
      # every chunk resolving in the parent PID with no warning at all.
      if (future::nbrOfWorkers() < 2L ||
          (inherits(future::plan(), "multicore") && !future::supportsMulticore())) {
        warning("parallel_backend = \"future\" but the active future plan resolves ",
                "in one process, so this will run serially. Set e.g. ",
                "future::plan(future::multicore, workers = ", ncore, ").",
                call. = FALSE)
        rp_count_serial_after_all()   # same reason as Windows: report what ran
      }
      # `ncore` was computed and then used only in the warning above, so a six-worker plan
      # ran six workers however small `workers` was. Chunking the jobs caps the number of
      # futures in flight at `ncore`, which is the only lever this backend gives us without
      # rewriting the caller's plan.
      untag(future.apply::future_lapply(
        idx_tagged, f_tagged, future.seed = NULL,
        future.chunk.size = ceiling(length(idx) / ncore)))
    },

    BiocParallel = {
      # MulticoreParam forks and is unavailable on Windows, where BiocParallel
      # itself substitutes a serial param, so this stays correct there.
      bp <- .combat_clusters$bpparam
      if (is.null(bp) || BiocParallel::bpnworkers(bp) != ncore) {
        bp <- BiocParallel::MulticoreParam(workers = ncore, stop.on.error = FALSE,
                                           RNGseed = NULL)
        .combat_clusters$bpparam <- bp
      }
      untag(BiocParallel::bplapply(idx_tagged, f_tagged, BPPARAM = bp))
    },

    foreach = {
      # foreach's backend is process-global and there is no exported way to read it back,
      # so do not overwrite one that already exists. If the caller registered a backend,
      # theirs is used and nothing global is touched. Only when none is registered do we
      # register our cached pool, and then we reset afterwards.
      #
      # An earlier attempt restored with `if (!prev) registerDoSEQ()`, which is wrong twice:
      # registerDoSEQ() itself makes getDoParRegistered() TRUE, so the restore worked only
      # on the first of a run's many dispatches, and a caller who HAD a backend had it
      # silently replaced by ours and never given back.
      # Testing getDoParRegistered() alone is not enough: registerDoSEQ() in the restore
      # below makes it TRUE, so from the second dispatch on, the branch was skipped, nothing
      # was registered, and %dopar% ran sequentially in the parent while the cached cluster
      # sat idle. Measured before this fix: 4 worker PIDs on dispatch 1, then 1 (the parent)
      # on every dispatch after. doSEQ is the sequential backend, so it counts as no backend.
      prev_backend <- if (foreach::getDoParRegistered()) foreach::getDoParName() else NULL
      own_backend <- is.null(prev_backend) || identical(prev_backend, "doSEQ")
      if (own_backend) {
        cl <- combat_cluster(ncore, packages = c("edgeR", "limma"))
        doParallel::registerDoParallel(cl)
        .combat_clusters$registered_by_us <- TRUE
        on.exit({
          try(foreach::registerDoSEQ(), silent = TRUE)
          .combat_clusters$registered_by_us <- NULL
        }, add = TRUE)
      }
      # the operator has to be bound locally: foreach is in Suggests, so it is not
      # imported, and `%dopar%` is not available by qualification alone
      `%dopar%` <- foreach::`%dopar%`
      i <- NULL  # keeps R CMD check quiet about the foreach index

      # A caller's own `doRNG` registration rewrites the MASTER's .Random.seed on every
      # %dopar%. ComBat-seq's own sample() lives in monte_carlo_int_NB on the serial side, so
      # a leaked stream is a different answer from the original under shrink = TRUE. Nothing
      # here consumes the stream, so restoring it costs nothing and makes the claim above
      # -- that dispatch never moves the caller through the random stream -- actually true.
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        seed0 <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
        on.exit(assign(".Random.seed", seed0, envir = globalenv()), add = TRUE)
      }

      # One task per chunk, whoever owns the backend. Grouping chunks to hold a caller's
      # backend down to `workers` was measured 2.08x slower on a wide registered backend,
      # made each worker retain its whole group's results, and clamped a remote backend to
      # this machine's core count. A backend the caller registered is theirs: its width
      # governs, `workers` does not, and that is documented rather than silently enforced.
      out <- foreach::foreach(i = idx_tagged, .packages = "edgeR",
                              .errorhandling = "pass") %dopar% f_tagged(i)
      untag(out)
    },

    stop("unhandled parallel_backend: ", parallel_backend,
         ". combat_backends() names it but combat_parallel_lapply() has no branch for it.",
         call. = FALSE)
  )
}


# ---- worker failure ----------------------------------------------------------

#' Check worker results and fail with the real cause
#'
#' Backends report failure in different shapes, and the difference matters.
#' `mclapply` returns a `try-error` carrying a condition for a thrown error but a
#' plain `NULL` for a *killed* child, the out-of-memory case. Calling
#' `conditionMessage()` on that NULL raises its own dispatch error and buries the
#' actual cause, which is exactly what an earlier version of this code did.
#' `future` and `BiocParallel` normally throw before reaching here; `foreach` with
#' `.errorhandling = "pass"` returns the condition object itself.
#'
#' @param parts Result list from `combat_parallel_lapply()`.
#' @param idx The chunk list the work was dispatched on. When supplied, each
#'   returned chunk is checked against the number of rows it was given.
#' @param what Label used in the error message.
#' @return `parts` unchanged when every element is usable.
#' @noRd
combat_parallel_check <- function(parts, what, idx = NULL) {
  is_err <- function(p) inherits(p, "try-error") || inherits(p, "condition")
  errored <- vapply(parts, is_err, logical(1))
  died <- vapply(parts, is.null, logical(1))

  # A chunk can also come back the wrong SIZE, which neither of the checks above sees. A
  # custom executor returning the right number of results with one of them short passed
  # validation and then introduced NA rows during the bind, silently. Rows are checked
  # where the result is a matrix or a vector; the GLM path returns a list of fit fields
  # and is skipped rather than guessed at.
  # Count first. `got != want` RECYCLES when the lengths differ, so three results against
  # four chunks compared FALSE FALSE FALSE FALSE and passed validation with only a warning.
  if (!is.null(idx) && length(parts) != length(idx)) {
    stop(what, ": the backend returned ", length(parts), " result(s) for ", length(idx),
         " chunk(s). A missing or extra chunk would be bound into the wrong rows.",
         call. = FALSE)
  }

  if (!is.null(idx) && !any(errored) && !any(died)) {
    want <- lengths(idx)
    got <- vapply(parts, function(p) {
      if (is.null(p) || (is.list(p) && !is.data.frame(p))) NA_integer_ else as.integer(NROW(p))
    }, integer(1))
    wrong <- !is.na(got) & got != want
    if (any(wrong)) {
      k <- which(wrong)[1]
      stop(what, ": chunk ", k, " came back with ", got[k], " row(s) where ", want[k],
           " were dispatched. A backend that returns a short, padded or reordered chunk ",
           "would corrupt the result rather than fail, so this is refused.", call. = FALSE)
    }
  }

  if (!any(errored) && !any(died)) return(parts)

  msg <- character(0)
  if (any(died)) {
    msg <- c(msg, sprintf("%d chunk(s) returned NULL, meaning the worker process died. %s",
                          sum(died),
                          "That is almost always the kernel killing it for memory: lower `workers`, raise `chunks`, or set options(combat.fork = FALSE)."))
  }
  if (any(errored)) {
    texts <- vapply(parts[errored], function(p) {
      cond <- if (inherits(p, "condition")) p else attr(p, "condition")
      if (is.null(cond)) "no condition attached" else conditionMessage(cond)
    }, character(1))
    detail <- texts[1]

    # EVERY chunk failing with the SAME text is definitive evidence of an environment problem
    # rather than a data one, and saying so is the difference between a two-minute fix and a
    # diagnostic pass over the cohort. A real consumer read "estimateGLMTagwiseDisp across
    # batches: 197 chunk(s) raised an error" as a data fault and went looking at the counts;
    # all 197 had failed identically because the package had been reinstalled underneath the
    # session. The message was accurate and still sent them the wrong way.
    same <- length(unique(texts)) == 1L && sum(errored) > 1L
    msg <- c(msg, if (same) {
      sprintf(paste0("all %d chunk(s) failed with the identical error, which points at the ",
                     "environment rather than the data: %s"), sum(errored), detail)
    } else {
      sprintf("%d chunk(s) raised an error, first was: %s", sum(errored), detail)
    })

    # and if the reason is the one this package can actually diagnose, say it outright
    if (isTRUE(rnaparallel_stale()) ||
        any(grepl("lazy-load database", texts, fixed = TRUE))) {
      msg <- c(msg, paste("rnaparallel was reinstalled while this session was running, so its",
                          "lazy-load database no longer matches what the session loaded.",
                          "Restart R; nothing about the data is wrong."))
    }
  }

  # A condition class lets a caller branch on structure instead of regexing prose, which is
  # what a caller deciding whether to fall back to uncorrected counts actually needs.
  stop(structure(
    class = c("rnaparallel_dispatch_error", "error", "condition"),
    list(message = paste0(what, ": ", paste(msg, collapse = " ")), call = NULL,
         stage = what, n_chunks = length(parts),
         n_errored = sum(errored), n_died = sum(died))))
}


# ---- the two parallelised paths ---------------------------------------------

#' Does this design send edgeR down its one-group kernel
#'
#' edgeR's glmFit dispatches on the design: a layout whose columns are exactly its factor
#' levels goes to mglmOneWay -> mglmOneGroup, and anything else to mglmLevenberg. The one-group
#' kernel is NOT a pure function of the gene it is fitting when that gene's fit does not
#' converge. Measured on edgeR 4.4.2: with low counts and one library over-sequenced 1000x, the
#' same gene's coefficient changes with which other rows share the matrix, and a solo fit
#' returns a subnormal like 2.37e-314. Splitting rows therefore changes the answer.
#'
#' There is no cheap per-gene detector. The bad values are finite, ordinary-looking, and carry
#' no convergence flag, so this is a design-level refusal rather than a result-level filter.
#' NULL, a vector, and an intercept-only matrix all count as one-group.
#' @noRd
combat_design_oneway <- function(design) {
  if (is.null(design)) return(TRUE)
  d <- as.matrix(design)
  if (ncol(d) < 1L) return(TRUE)
  ok <- tryCatch(nlevels(edgeR::designAsFactor(d)) == ncol(d), error = function(e) TRUE)
  isTRUE(ok)
}


#' @noRd
rp_rows <- function(x, ii) if (is.null(x) || is.null(dim(x))) x else x[ii, , drop = FALSE]

#' @noRd
rp_per_gene <- function(x, ii) {
  if (is.null(x) || (is.null(dim(x)) && length(x) > 1L)) x[ii] else rp_rows(x, ii)
}

#' Row-parallel edgeR GLM fit
#'
#' Separable because `edgeR:::.compressOffsets` uses the offset it is handed
#' rather than recomputing `log(colSums(y))`, and `addPriorCount`,
#' `mglmLevenberg` and `mglmOneWay` are per gene from there on. Offset and
#' dispersion must therefore always arrive explicitly: a NULL offset inside a
#' worker would silently rebuild library sizes from that worker's slice of genes,
#' which changes every fitted value without raising anything.
#'
#' @section Return value is not a DGEGLM:
#' This returns a plain list carrying the fields ComBat-seq reads and
#' nothing else. It is correct at that one call site and a trap anywhere else:
#' `glmLRT()` or `glmQLFTest()` on it would fail, or worse half-work. The
#' function stays internal for that reason.
#'
#' @param y Count matrix, genes in rows.
#' @param design Design matrix.
#' @param dispersion Scalar, or a per-gene vector or matrix.
#' @param offset Offset matrix or vector. Never NULL when called from a DGEList.
#' @param weights Optional observation weights.
#' @param prior.count Prior count passed through to edgeR.
#' @param start Optional starting coefficients.
#' @param workers Maximum concurrent workers.
#' @param chunks Row chunks. Defaults to `workers`.
#' @param parallel_backend One of [combat_backends()].
#' @return A list of fit components.
#' @noRd
glmFit_rows_parallel <- function(y, design, dispersion, offset, weights = NULL,
                                 prior.count = 0.125, start = NULL,
                                 workers = 4L, chunks = NULL,
                                 parallel_backend = getOption("combat.backend", combat_default_backend())) {
  y <- as.matrix(y)
  idx <- combat_row_chunks(nrow(y), workers = workers, chunks = chunks)

  # Without an explicit offset each worker rebuilds library sizes from its own genes, which
  # moved coefficients by 1.8 in testing and raises nothing. Refuse rather than guess.
  if (is.null(offset) && length(idx) > 1L) {
    stop("glmFit_rows_parallel needs an explicit offset: a worker would rebuild library ",
         "sizes from its own slice of genes and every fitted value would change silently.",
         call. = FALSE)
  }

  # matrices are row-subset; per-sample vectors and scalars pass through untouched
  # A per-gene dispersion is a bare vector with no dim, so rows() handed every chunk the
  # full-length vector and edgeR rejected it: worked at chunks = 1, failed at chunks = 4.
  # The documented "scalar, or a per-gene vector or matrix" only held for a single chunk.
  # Serial edgeR refuses a per-gene vector whose length is neither 1 nor nrow(y); the split
  # path used to slice it anyway and hand each chunk a correctly sized piece of the wrong
  # thing. Refuse the same inputs the serial call refuses.
  check_len <- function(x, nm) {
    if (!is.null(x) && is.null(dim(x)) && length(x) > 1L && length(x) != nrow(y)) {
      stop(nm, " has wrong length: ", length(x), " for ", nrow(y), " rows", call. = FALSE)
    }
    invisible(TRUE)
  }
  check_len(dispersion, "dispersion")


  # A one-group design reaches a kernel whose output depends on the block it was fitted in, so
  # this path is not row-splittable at any chunk count. Called whole instead. Measured before
  # the gate: ComBat_seq_parallel returned 1, 2 and 3 differing cells at chunks 2, 4 and 8 on
  # a 240x12 matrix with one library over-sequenced 1000x.
  if (combat_design_oneway(design)) {
    return(edgeR::glmFit.default(y, design = design, dispersion = dispersion, offset = offset,
                                 lib.size = NULL, weights = weights,
                                 prior.count = prior.count, start = start))
  }

  # Rebuilt against an environment holding only what the body reads. A closure is serialised
  # WITH its defining environment, so on a socket backend this frame's live bindings and its
  # unforced promises travel with every task, and the promises reach back through the original's
  # frames into the entry point's raw inputs. Invisible on a forking backend, where the child
  # inherits the pages, which is why it survived this long.
  # Measured on 1,200 x 60: the two glmFit dispatches went 5,560,242 B to 1,201,690 B and
  # 8,013,222 B to 1,777,238 B.
  .lean <- new.env(parent = parent.env(environment()))
  .lean$y <- y; .lean$design <- design; .lean$dispersion <- dispersion
  .lean$offset <- offset; .lean$weights <- weights
  .lean$prior.count <- prior.count; .lean$start <- start
  fit_rows <- function(ii) {
    f <- edgeR::glmFit.default(
      y[ii, , drop = FALSE], design = design, dispersion = rp_per_gene(dispersion, ii),
      offset = rp_rows(offset, ii), lib.size = NULL, weights = rp_rows(weights, ii),
      prior.count = prior.count, start = rp_rows(start, ii))
    # only the six fields the parent reads back. edgeR also returns full counts, dispersion
    # and offset slices that are discarded on arrival, about 540 MB a dispatch on a cohort
    # this size, serialised through the pipe and held in the parent for nothing.
    list(coefficients = f$coefficients, fitted.values = f$fitted.values,
         df.residual = f$df.residual, unshrunk.coefficients = f$unshrunk.coefficients,
         method = f$method, failed = f$failed)
  }
  environment(fit_rows) <- .lean

  fits <- combat_parallel_check(
    combat_parallel_lapply(idx, fit_rows, workers, parallel_backend, cells = length(y),
                           min_cells = getOption("combat.min.glm.cells", 1e5)),
    "glmFit_rows_parallel", idx)

  # A gene whose fit failed carries state from the block it sat in, so one failure anywhere
  # invalidates the split. edgeR reports it, so this is checked rather than assumed.
  if (any(vapply(fits, function(f) any(f$failed != 0), logical(1)))) {
    return(edgeR::glmFit.default(y, design = design, dispersion = dispersion, offset = offset,
                                 lib.size = NULL, weights = weights,
                                 prior.count = prior.count, start = start))
  }

  # chunks are interleaved, so every gene-indexed result comes back permuted
  ord <- combat_row_order(idx)
  # `ord` is a permutation of seq_len(nrow(y)), and the only sorted permutation of 1:n is 1:n,
  # so an unsorted test decides exactly whether the gather moves anything. At one chunk both
  # the rbind and the gather are no-ops on the values and pure copies in cost: measured on an
  # 18,270 x 1,500 field, 136 ms for the rbind and 53 ms for the permute, both removed. That is
  # the documented `workers = 1` and `chunks = 1` path, which paid for two full copies of every
  # bound field to return what it was given. is.matrix guards the single-chunk shortcut so a
  # vector field still goes through rbind, which would legitimately promote it to one row.
  bind <- function(nm) {
    pieces <- lapply(fits, function(f) f[[nm]])
    m <- if (length(pieces) == 1L && is.matrix(pieces[[1L]])) pieces[[1L]]
         else do.call(rbind, pieces)
    if (nrow(m) != nrow(y)) {
      stop("glmFit_rows_parallel: field '", nm, "' bound to ", nrow(m), " rows where ",
           nrow(y), " were dispatched", call. = FALSE)
    }
    if (is.unsorted(ord)) m[ord, , drop = FALSE] else m
  }
  vec <- function(nm) {
    v <- unlist(lapply(fits, function(f) f[[nm]]), use.names = FALSE)
    if (length(v) != nrow(y)) {
      stop("glmFit_rows_parallel: field '", nm, "' bound to ", length(v), " values where ",
           nrow(y), " were dispatched", call. = FALSE)
    }
    if (is.unsorted(ord)) v[ord] else v
  }
  # `deviance` is deliberately absent. It is the ONE field edgeR returns that is not a pure
  # function of its own gene: in the Levenberg branch `mglmLevenberg` passes its raw slot
  # through, and for a gene whose fit aborts that slot holds a neighbour's value. Measured on
  # an 803 x 27 fit, gene 34 came back 109.70 whole, 75.60 at 3 chunks, 14.52 at 5, and
  # 2.47e-323 fitted alone. Paired with gene 1 it returned gene 1's deviance. Serial and
  # forked agree, so it is edgeR's own behaviour and not a parallel artefact.
  #
  # ComBat-seq reads coefficients and fitted.values only, so nothing is lost by omitting it,
  # and omitting it keeps the promise this file makes elsewhere: every field returned here is
  # chunk-count independent. Returning a value that changes with `chunks` would break that
  # promise quietly. If a future caller needs it, recompute from the bound fitted values with
  # edgeR::nbinomDeviance() rather than binding the per-chunk slots.
  out <- list(coefficients = bind("coefficients"),
              fitted.values = bind("fitted.values"),
              df.residual = vec("df.residual"),
              method = fits[[1]]$method,
              design = design,
              prior.count = prior.count)
  if (prior.count > 0) out$unshrunk.coefficients <- bind("unshrunk.coefficients")
  out
}

# Derived from sva (Zhang, Parmigiani and Johnson), Artistic-2.0.
#
# The sva 3.54.0 `match_quantiles` body, deparsed at width.cutoff = 500. The only place
# this package holds original code, so it is pinned:
# `match_quantiles_rows` below is a transcription of exactly this text and runs only while the
# backend still deparses to it.
.match_quantiles_pinned <- c(
  "{",
  "    new_counts_sub <- matrix(NA, nrow = nrow(counts_sub), ncol = ncol(counts_sub))",
  "    for (a in 1:nrow(counts_sub)) {",
  "        for (b in 1:ncol(counts_sub)) {",
  "            if (counts_sub[a, b] <= 1) {",
  "                new_counts_sub[a, b] <- counts_sub[a, b]",
  "            }",
  "            else {",
  "                tmp_p <- pnbinom(counts_sub[a, b] - 1, mu = old_mu[a, b], size = 1/old_phi[a])",
  "                if (abs(tmp_p - 1) < 1e-04) {",
  "                  new_counts_sub[a, b] <- counts_sub[a, b]",
  "                }",
  "                else {",
  "                  new_counts_sub[a, b] <- 1 + qnbinom(tmp_p, mu = new_mu[a, b], size = 1/new_phi[a])",
  "                }",
  "            }",
  "        }",
  "    }",
  "    return(new_counts_sub)",
  "}")

#' Whole-slice transcription of the pinned `match_quantiles` body
#'
#' The same three branches in the same order, taken over every cell at once instead of one
#' cell at a time. `pnbinom` and `qnbinom` are vectorised C and each cell reads only its own
#' arguments, so one call over the selected cells computes exactly what the cell loop computes.
#' The `1 +` on the `qnbinom` branch is part of the original body and is easy to drop when
#' transcribing from a description of it rather than from the body itself.
#'
#' `old_phi` and `new_phi` are per GENE, so the row each selected cell belongs to has to be
#' recovered from its linear index to look the dispersion up. That is the only arithmetic here
#' the cell loop does not do, and it indexes rather than computes.
#'
#' The logical `NA` start matters: the original's result type is whatever its assignments promote
#' that matrix to, integer while every cell keeps its count and double once one cell takes the
#' `qnbinom` branch. `out[] <- counts_sub` promotes identically and, unlike `out <- counts_sub`,
#' keeps the original's absent dimnames rather than carrying the input's.
#'
#' This replaced a per-row loop, which is why the type and dimnames notes above are stated
#' rather than assumed: measured 0.508 s to 0.370 s on an 18,270 x 28 slice, and `identical()`
#' to `sva::match_quantiles` on 600 randomised cases spanning integer and double storage,
#' dimnames present and absent, all-counts-<=1 inputs and dispersions down to 1e-8.
#'
#' Reached only through `combat_mq_dispatch()`.
#' @noRd
match_quantiles_rows <- function(counts_sub, old_mu, old_phi, new_mu, new_phi) {
  out <- matrix(NA, nrow = nrow(counts_sub), ncol = ncol(counts_sub))
  out[] <- counts_sub
  big <- which(counts_sub > 1)
  if (length(big)) {
    r <- ((big - 1L) %% nrow(counts_sub)) + 1L        # gene each selected cell belongs to
    tmp_p <- stats::pnbinom(counts_sub[big] - 1, mu = old_mu[big], size = 1/old_phi[r])
    q <- which(abs(tmp_p - 1) >= 1e-04)
    if (length(q)) {
      out[big[q]] <- 1 + stats::qnbinom(tmp_p[q], mu = new_mu[big[q]],
                                        size = 1/new_phi[r[q]])
    }
  }
  out
}

#' Choose the quantile matcher for one call
#'
#' Two refusals, both falling back to the backend's own function, which is always correct
#' because it is the thing being matched.
#'
#' The body gate is the price of holding a copy at all. A transcription cannot follow an
#' upstream edit, so the backend's deparsed body must still equal `.match_quantiles_pinned`
#' byte for byte; one changed character and the slice goes back to the original.
#'
#' The NA gate exists because the original ERRORS on a missing count: `if (counts_sub[a, b] <=
#' 1)` is `if (NA)`. `which()` drops NA instead, so the row form would return an NA cell and
#' no error at all. `old_mu` and `old_phi` reach that same condition through `tmp_p`, so they
#' are checked too. Falling back preserves the original's own error message.
#' @noRd
combat_mq_dispatch <- function(mq, counts_sub, old_mu, old_phi) {
  # NULL means call the original WHOLE: these inputs reach the original's own `if (NA)` error or
  # a degenerate shape, and running it unsliced preserves that behaviour byte for byte,
  # error message included. `old_phi <= 0` is included because `size = 1/old_phi` turns
  # non-positive dispersions into NaN probabilities the row form would silently keep.
  # is.finite has no data.frame method, and the original accepts data.frame counts, so a
  # non-matrix goes to the original whole before anything here can touch it.
  #
  # Three conditions below are not reachable from ComBat-seq, whose `mu_hat` is a matrix of
  # fitted values and whose `phi` always has one entry per gene. They are here because this is
  # the gate's contract, not ComBat-seq's: a NEGATIVE finite `old_mu`, an `old_phi` shorter
  # than the matrix, or an `old_mu` whose dim differs all make `pnbinom` return NaN or read the
  # wrong cell, and the original ERRORS on that (`if (abs(NaN - 1) < 1e-04)` is `if (NA)`) while
  # the vectorised form's `which()` drops the NA and returns a plausible half-matched matrix
  # with no signal. The whole point of falling back is to preserve the original's own behaviour,
  # so the gate has to be complete rather than complete-for-today's-only-caller.
  #
  # The finiteness tests are reductions rather than `all(is.finite(x))`, which allocates a
  # logical the size of the matrix on the path that is 66% of ComBat-seq's serial time. A
  # numeric x is all-finite exactly when it holds no NA and its min and max are both finite.
  finite_all <- function(x) !anyNA(x) && is.finite(min(x)) && is.finite(max(x))
  if (!is.matrix(counts_sub) || !is.matrix(old_mu) ||
      nrow(counts_sub) == 0L || ncol(counts_sub) == 0L ||
      !identical(dim(old_mu), dim(counts_sub)) ||
      length(old_phi) != nrow(counts_sub) ||
      !finite_all(counts_sub) || !finite_all(old_mu) ||
      !finite_all(old_phi) || any(old_phi <= 0) || any(old_mu < 0)) {
    return(NULL)
  }
  # deparse honours options(scipen), so an analyst's .Rprofile could silently shut this gate
  # and hand every slice back to the cell loop with no signal. Pinned to the setting the
  # text was captured under; the comparison is over the body, not the print options.
  op <- options(scipen = 0L); on.exit(options(op), add = TRUE)
  if (!identical(deparse(body(mq), width.cutoff = 500L), .match_quantiles_pinned)) {
    # the gate did its job, and that is exactly why it has to be visible: correct numbers at
    # original speed is indistinguishable from correct numbers at our speed
    rp_note_fallback("match_quantiles")
    return(mq)
  }
  match_quantiles_rows
}

#' Row-parallel quantile matching
#'
#' The dominant cost in ComBat-seq. Profiling a 10,000 by 500 run with `group = NULL` puts it
#' at 66.5% of serial time, against 14.0% for tagwise dispersion and 13.1% for common
#' dispersion.
#'
#' Splitting by row is exact rather than approximate. The original ComBat-seq body is a cell
#' loop where `new_counts_sub[a, b]` reads only `counts_sub[a, b]`,
#' `old_mu[a, b]`, `new_mu[a, b]`, `old_phi[a]` and `new_phi[a]`, so every cell
#' depends on its own gene and nothing else. No cross-row term exists to lose.
#'
#' Each slice is matched by `match_quantiles_rows()`, a row-vectorised transcription of the
#' sva 3.54.0 body and the one copy of original code this package holds. Drift is gated, not
#' assumed away: `combat_mq_dispatch()` compares the backend's body against the pinned text
#' byte for byte and hands the slice back to the backend's own function the moment they
#' differ, or the moment an input carries an NA.
#'
#' @param mq The backend's `match_quantiles`, from `combat_backend()`.
#' @param counts_sub Count matrix for the genes being adjusted.
#' @param old_mu,new_mu Fitted means, same shape as `counts_sub`.
#' @param old_phi,new_phi Per-gene dispersions.
#' @param workers Maximum concurrent workers.
#' @param chunks Row chunks. Defaults to `workers`.
#' @param parallel_backend One of [combat_backends()].
#' @return A matrix the same shape as `counts_sub`.
#' @noRd
match_quantiles_parallel <- function(mq, counts_sub, old_mu, old_phi, new_mu, new_phi,
                                     workers = 4L, chunks = NULL,
                                     parallel_backend = getOption("combat.backend", combat_default_backend())) {
  idx <- combat_row_chunks(nrow(counts_sub), workers = workers, chunks = chunks)

  # gate resolved once, before the fork, so a worker inherits the decision rather than
  # re-deparsing the backend body once per chunk
  mq_fun <- combat_mq_dispatch(mq, counts_sub, old_mu, old_phi)
  if (is.null(mq_fun)) return(mq(counts_sub, old_mu, old_phi, new_mu, new_phi))

  # Slicing inside the worker, not before dispatch. Shipping each chunk a pre-sliced payload
  # was tried and measured on Windows/PSOCK, where the closure below is serialised rather than
  # inherited through copy-on-write: +3.9%, -6.6%, +4.3%, +8.0% at 2, 4, 6 and 8 workers,
  # against 6.4% drift in the serial reference over the same run. That is noise, and one arm
  # was slower. `future.apply` resolves the globals a closure actually reads rather than
  # shipping its whole frame, so the cost this was aimed at was not being paid to begin with.
  # Rebuilt against an environment holding only what the body reads. A closure is serialised
  # WITH its defining environment, so on a socket backend this frame's live bindings and its
  # unforced promises travel with every task, and the promises reach back through the original's
  # frames into the entry point's raw inputs. Invisible on a forking backend, where the child
  # inherits the pages, which is why it survived this long.
  # Measured on 1,200 x 60 counts: 8,603,855 B shipped per task against 662,975 B after.
  .lean <- new.env(parent = parent.env(environment()))
  .lean$mq_fun <- mq_fun; .lean$counts_sub <- counts_sub; .lean$old_mu <- old_mu
  .lean$old_phi <- old_phi; .lean$new_mu <- new_mu; .lean$new_phi <- new_phi
  do_rows <- function(ii) mq_fun(counts_sub = counts_sub[ii, , drop = FALSE],
                                 old_mu = old_mu[ii, , drop = FALSE],
                                 old_phi = old_phi[ii],
                                 new_mu = new_mu[ii, , drop = FALSE],
                                 new_phi = new_phi[ii])
  environment(do_rows) <- .lean

  parts <- combat_parallel_check(
    combat_parallel_lapply(idx, do_rows, workers, parallel_backend, cells = length(counts_sub)),
    "match_quantiles_parallel", idx)
  # Same argument as glmFit_rows_parallel's bind(): one chunk makes both the rbind and the
  # gather no-ops on the values, and a sorted `ord` makes the gather one too.
  m <- if (length(parts) == 1L && is.matrix(parts[[1L]])) parts[[1L]]
       else do.call(rbind, parts)
  ord <- combat_row_order(idx)
  if (is.unsorted(ord)) m[ord, , drop = FALSE] else m
}

#' Row-parallel tagwise dispersion estimation
#'
#' The third hot path, and the one that used to be called unparallelisable. ComBat-seq
#' estimates dispersion once per batch, which is 14.0% of serial time on the same
#' 10,000 by 500 profile and grows with
#' batch count, so on a 100-batch design it is a large serial block sitting in the middle
#' of an otherwise parallel run.
#'
#' @section Why this is exact, and only at `prior.df = 0`:
#' `dispCoxReidInterpolateTagwise` moderates each gene's adjusted profile likelihood
#' towards a smoothed curve computed across genes, which is a genuine cross-row term:
#' `(apl + prior.n * apl.smooth) / (1 + prior.n)`. ComBat-seq calls this with
#' `prior.df = 0`, and `prior.n` is `prior.df / (nlibs - ncoefs)`, so the moderation
#' weight is exactly zero and the expression collapses to `apl`, so each gene's dispersion
#' depends on its own counts alone. That holds only while every gene's `apl` is FINITE:
#' `0 * NaN` and `0 * -Inf` are both `NaN`, and `apl.smooth` is built across genes, so one
#' poisoned gene would contaminate a different neighbour set in each arm. Counts at the top
#' of the double range can do it, which is why the gate below also checks magnitude. Verified rather than assumed: chunked and whole
#' matrix results were `identical()` on 2,000 and 8,000 genes, with and without dead,
#' constant and near-empty genes present. Any other `prior.df` is handed straight to
#' edgeR unsplit.
#'
#' `trend` is left alone rather than forced off. It defaults to `TRUE` and does drive a
#' moving average across genes, but that feeds `apl.smooth`, which the zero weight
#' discards, so splitting stays exact either way.
#'
#' @section What must be passed explicitly:
#' Three quantities are functions of the whole gene set and would otherwise be rebuilt
#' from a worker's own slice, silently: `offset`, which defaults to `log(colSums(y))`;
#' `span`, which defaults to `(10/ntags)^0.23` and therefore changes with chunk size; and
#' `AveLogCPM`. Computing all three once up front and passing them down is the same
#' discipline `glmFit_rows_parallel()` applies to offsets. A split that omits them returns
#' different numbers on every matrix tested.
#'
#' @param y Count matrix for one batch, genes in rows.
#' @param design Design matrix for that batch's samples.
#' @param dispersion Common dispersion, scalar or per gene.
#' @param offset Offset vector or matrix. Computed once here when NULL.
#' @param prior.df Prior degrees of freedom. Anything other than 0 falls back to edgeR.
#' @param trend Passed through to edgeR unchanged.
#' @param span Smoothing span. Computed once here when NULL.
#' @param AveLogCPM Average log CPM per gene. Computed once here when NULL.
#' @param weights Optional observation weights.
#' @param workers Maximum concurrent workers.
#' @param chunks Row chunks. Defaults to `workers`.
#' @param parallel_backend One of [combat_backends()].
#' @return Numeric vector of per-gene dispersions, in the input row order.
#' @noRd
estimateGLMTagwiseDisp_rows_parallel <- function(y, design = NULL, dispersion = NULL,
                                                 offset = NULL, prior.df = 10, trend = TRUE,
                                                 span = NULL, AveLogCPM = NULL, weights = NULL,
                                                 workers = 4L, chunks = NULL,
                                                 parallel_backend = getOption("combat.backend", combat_default_backend())) {
  y <- as.matrix(y)
  ntag <- nrow(y)
  # serial edgeR stops on a dispersion whose length is neither 1 nor nrow(y). Splitting by row
  # would hand each chunk a valid-looking slice of an invalid vector and return numbers.
  if (!is.null(dispersion) && is.null(dim(dispersion)) &&
      length(dispersion) > 1L && length(dispersion) != ntag) {
    stop("dispersion has wrong length: ", length(dispersion), " for ", ntag, " rows",
         call. = FALSE)
  }

  # gate: only the zero-moderation case is provably row-separable, so anything else
  # goes to edgeR whole rather than being split on an assumption
  # A one-group design is excluded for the same reason glmFit_rows_parallel excludes it:
  # adjustedProfileLik fits with glmFit inside every chunk, so the one-group kernel's
  # block dependence reaches the returned dispersions as finite, plausible, wrong numbers.
  # The existing non-finite post-check cannot see them.
  # Design test FIRST. `&&` is order-independent for operands with no side effects, and these
  # have none, so the decision is unchanged. What changes is that the two full-slice scans are
  # no longer paid to reach a veto: ComBat-seq hands this a per-batch design of mod[batch, ],
  # which is one-way on every batch of every run, so `separable` is always FALSE here and both
  # scans were computed and discarded. Measured on an 18,270 x 28 slice, the shape a 54-batch
  # cohort passes: 1.66 ms to 0.033 ms, plus about 6 MB of transient allocation per batch that
  # was charged to the worker's RSS.
  #
  # The surviving scan is one pass and no allocation. `all(is.finite(y))` builds an n x m
  # logical and `abs(y)` an n x m double; anyNA plus max plus min answer the same question
  # about a numeric y, since a non-finite value is NA, NaN, Inf or -Inf and the three tests
  # between them catch all four.
  separable <- is.numeric(prior.df) && length(prior.df) == 1L &&
    !is.na(prior.df) && prior.df == 0 &&
    !combat_design_oneway(design) &&
    !anyNA(y) && max(y) < 1e150 && min(y) > -1e150
  if (!separable || ntag < 2L) {
    return(edgeR::estimateGLMTagwiseDisp(y, design = design, offset = offset,
                                         dispersion = dispersion, prior.df = prior.df,
                                         trend = trend, span = span,
                                         AveLogCPM = AveLogCPM, weights = weights))
  }

  # the three whole-matrix quantities, computed once, exactly as edgeR would
  if (is.null(offset)) offset <- log(colSums(y))
  if (is.null(span)) span <- if (ntag > 10) (10 / ntag)^0.23 else 1
  if (is.null(AveLogCPM)) AveLogCPM <- edgeR::aveLogCPM(y, offset = offset, weights = weights)

  idx <- combat_row_chunks(ntag, workers = workers, chunks = chunks)

  # Rebuilt against an environment holding only what the body reads. A closure is serialised
  # WITH its defining environment, so on a socket backend this frame's live bindings and its
  # unforced promises travel with every task, and the promises reach back through the original's
  # frames into the entry point's raw inputs. Invisible on a forking backend, where the child
  # inherits the pages, which is why it survived this long.
  .lean <- new.env(parent = parent.env(environment()))
  .lean$y <- y; .lean$design <- design; .lean$offset <- offset
  .lean$dispersion <- dispersion; .lean$trend <- trend; .lean$span <- span
  .lean$AveLogCPM <- AveLogCPM; .lean$weights <- weights
  disp_rows <- function(ii) edgeR::estimateGLMTagwiseDisp(
    y[ii, , drop = FALSE], design = design, offset = rp_rows(offset, ii),
    dispersion = rp_per_gene(dispersion, ii), prior.df = 0, trend = trend,
    span = span, AveLogCPM = AveLogCPM[ii], weights = rp_rows(weights, ii))
  environment(disp_rows) <- .lean

  # higher threshold than the other two paths: measured on 10-column slices this path
  # returns 0.65x at 10,000 cells, 0.79x at 20,000, 1.08x at 30,000 and 1.44x at 50,000
  parts <- combat_parallel_check(
    combat_parallel_lapply(idx, disp_rows, workers, parallel_backend, cells = length(y),
                           min_cells = getOption("combat.min.disp.cells", 3e4)),
    "estimateGLMTagwiseDisp_rows_parallel", idx)
  out <- unlist(parts, use.names = FALSE)[combat_row_order(idx)]

  # The separability gate above checks `y`, not the adjusted profile likelihood computed from
  # it. A degenerate fit could in principle return a non-finite apl for finite counts, and
  # `0 * NaN` is NaN, whose spread depends on AveLogCPM neighbours and therefore on the chunk
  # layout. Recompute unsplit rather than return a layout-dependent answer.
  if (!all(is.finite(out))) {
    return(edgeR::estimateGLMTagwiseDisp(y, design = design, offset = offset,
                                         dispersion = dispersion, prior.df = prior.df,
                                         trend = trend, span = span,
                                         AveLogCPM = AveLogCPM, weights = weights))
  }
  out
}
