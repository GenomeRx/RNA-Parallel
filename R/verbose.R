## verbose.R
##
## One elapsed line per companion call, so a long run says what it is doing and a stuck one
## says where it stopped. Nothing here touches a returned value: the timer is an on.exit hook
## and the quieting is a stdout sink, so both unwind on an error exactly as they do on success.


# ---- options -----------------------------------------------------------------

#' Read a single logical option, refusing a garbage value rather than ignoring it
#'
#' The `combat.min.*` options already refuse a typo instead of silently switching a gate off.
#' Same reasoning: someone who sets `combat.timing = "yes"` wants timing, and quietly giving
#' them silence is the opposite of what they asked for.
#' @noRd
rp_opt_flag <- function(name, default = FALSE) {
  v <- getOption(name, default)
  if (!(is.logical(v) && length(v) == 1L && !is.na(v))) {
    stop("`", name, "` must be TRUE or FALSE; got ", deparse(v), call. = FALSE)
  }
  v
}

#' @noRd
rp_opt_secs <- function(name, default = 0) {
  v <- suppressWarnings(as.numeric(getOption(name, default)))
  if (length(v) != 1L || is.na(v) || v < 0) {
    stop("`", name, "` must be a single non-negative number of seconds; got ",
         deparse(getOption(name)), call. = FALSE)
  }
  v
}


# ---- silencing the original ----------------------------------------------------

# sva::ComBat_seq announces itself with cat(), not message(): "Found 154 batches", "Estimating
# dispersions", "Fitting the GLM model" and five more, once per call. suppressMessages() cannot
# touch any of it, because none of it is a condition. Only a stdout sink can.
#
# stdout ONLY, deliberately. Sinking stderr too would swallow warnings, and the one thing worse
# than a chatty run is a silent one that dropped a warning about the data.

#' @noRd
rp_quiet_begin <- function() {
  con <- file(nullfile(), open = "w")
  sink(con, type = "output")
  con
}

#' @noRd
rp_quiet_end <- function(con) {
  # Pop only our own level. A caller who had their own sink open keeps it, and a call that
  # errored between sink() and here still unwinds, because this runs from on.exit().
  if (sink.number() > 0L) try(sink(type = "output"), silent = TRUE)
  try(close(con), silent = TRUE)
  invisible(NULL)
}


# ---- what actually ran -------------------------------------------------------

# "mclapply x6" is what was REQUESTED. Whether the work reached a worker is a different
# question and the more useful one: seven size gates can send a call serial, and a companion
# that returned identical() output at serial pace looks exactly like one that forked. So count
# both and report what happened rather than what was asked for.
.rp_dispatch <- new.env(parent = emptyenv())

#' @noRd
rp_count_reset <- function() {
  .rp_dispatch$par <- 0L
  .rp_dispatch$ser <- 0L
  .rp_dispatch$fallback <- character()
  .rp_dispatch$progress_last <- NULL
  invisible(NULL)
}


#' @noRd
rp_count <- function(parallel) {
  if (isTRUE(parallel)) .rp_dispatch$par <- rp_or0(.rp_dispatch$par) + 1L
  else                  .rp_dispatch$ser <- rp_or0(.rp_dispatch$ser) + 1L
  rp_progress_tick()
  invisible(NULL)
}

#' Reclassify the dispatch just counted as parallel: it fell back to serial after all
#'
#' The count is taken before the backend switch, because that is the one place every backend
#' passes through. Two arms then discover they cannot actually fork -- mclapply on Windows, and
#' a `future` plan that resolves in one process -- and both already say so in a warning. Without
#' this the timing line printed "future x2" beside the package's own warning that it had just
#' run serially, which is the exact claim the engine column exists to make impossible.
#' @noRd
rp_count_serial_after_all <- function() {
  .rp_dispatch$par <- max(0L, rp_or0(.rp_dispatch$par) - 1L)
  .rp_dispatch$ser <- rp_or0(.rp_dispatch$ser) + 1L
  invisible(NULL)
}

#' Record that a pinned original excerpt stood down for this call
#'
#' `match_quantiles` is the one piece of original text this package still holds, and its gate is a
#' byte-exact body match against sva. An sva release that reformats that body -- never mind
#' rewrites it -- silently hands every slice back to the original's own cell loop. The numbers stay
#' correct and the run gets 1.4-1.7x slower forever, which is the failure nobody reports because
#' nobody can see it. So say it in the same line that already reports what ran.
#' @noRd
rp_note_fallback <- function(what) {
  cur <- .rp_dispatch$fallback
  if (!(what %in% cur)) .rp_dispatch$fallback <- c(cur, what)
  invisible(NULL)
}

#' @noRd
rp_or0 <- function(x) if (is.null(x)) 0L else x


# ---- single-line progress (default on) ----------------------------------------

# combat.timing prints one line at the END of a call. ComBat-seq alone dispatches its hot
# paths up to 2*n_batch + 3 times per call, so on a large cohort (hundreds of batches) there
# is nothing on screen between "computing" and the final line, and a stuck run looks exactly
# like a slow one. This is the fix: one line, overwritten in place with a carriage return, so
# it never scrolls and never floods a log. On by default so every parallel call is visible
# without opting in; set options(combat.progress = FALSE) to silence it.
#
# This tick fires in the MASTER process between dispatches -- see the file-based mechanism
# below for the part that survives the master blocking inside one big parallel call.

#' @noRd
rp_progress_tick <- function() {
  if (!rp_opt_flag("combat.progress", default = TRUE)) return(invisible(NULL))
  # Throttled to 4/sec: enough to prove the run is alive, not enough to slow it down or
  # flood a log file that doesn't understand a bare carriage return (each tick still costs
  # a Sys.time() read).
  now <- Sys.time()
  last <- .rp_dispatch$progress_last
  if (!is.null(last) && as.numeric(now - last, units = "secs") < 0.25) return(invisible(NULL))
  .rp_dispatch$progress_last <- now

  n <- rp_or0(.rp_dispatch$par) + rp_or0(.rp_dispatch$ser)
  label <- .rp_dispatch$progress_label %||% "dispatch"
  # No total is known in advance (batch count varies by call site), so this counts up rather
  # than filling a bar to a percentage that would have to guess. A rising count is still
  # unambiguous evidence of life; a percentage that never appears is worse than none.
  cr <- "\r"
  cat(cr, sprintf("  %s: %d dispatched ", label, n), sep = "", file = stderr())
  utils::flush.console()
  invisible(NULL)
}

#' @noRd
rp_progress_done <- function() {
  if (!rp_opt_flag("combat.progress", default = TRUE)) return(invisible(NULL))
  # Clear the line rather than leaving a stale count sitting there once the step's own
  # combat.timing line (if any) prints below it.
  cr <- "\r"
  cat(cr, strrep(" ", 60L), cr, sep = "", file = stderr())
  utils::flush.console()
  invisible(NULL)
}


# ---- file progress (opt-in, for a blocking parallel call) --------------------

# The console tick above only fires in the MASTER process, and only between dispatches: the
# moment mclapply/future/BiocParallel/foreach blocks for the actual parallel work, the master
# is synchronously waiting and cannot print anything until it returns. On a single large
# dispatch (one ComBat-seq stage split into a few hundred chunks across 16 workers) that block
# can run for hours, and the tick above goes silent for the whole stretch: exactly the gap
# that made a live run indistinguishable from a hung one on a real 14h53m + 13h pooled run.
#
# Workers cannot fix this by writing to the master's console either: a forked child's stdout
# is not reliably multiplexed back to an RStudio Server session, and PSOCK/BiocParallel
# workers do not share a console at all. A file each worker can append to, read from a
# SEPARATE session while the master blocks, is the only channel that survives all four
# backends. One file per worker PID avoids write contention between workers.
#
# Off unless the caller sets a directory. Every write is one line; the cost is a
# file-append syscall per chunk, not per gene, so at hundreds of chunks over hours it is
# immaterial next to the compute itself.

#' @noRd
rp_progress_dir <- function() {
  d <- getOption("combat.progress.dir", NA_character_)
  if (is.na(d) || !nzchar(d)) return(NULL)
  d
}

#' Append one chunk-progress line to this worker's own file
#'
#' TSV so it parses without guessing a delimiter: unix time, stage label, chunk index,
#' event ("start" or "done"). One file per worker PID means every writer only ever appends
#' to a file nothing else touches, so no locking is needed on any of the four backends,
#' including PSOCK workers that share nothing with each other. `dir` is passed explicitly
#' rather than read via `rp_progress_dir()`/`getOption()`, because the one caller that matters
#' (the dispatch wrapper in `combat_parallel_lapply()`) already resolved it once in the master
#' and captured it as a plain value in the worker's closure -- a PSOCK worker does not inherit
#' the master's `options()`, so reading the option again inside the worker would silently see
#' the default instead.
#' @noRd
rp_progress_file_write <- function(dir, stage, chunk, event) {
  if (is.null(dir)) return(invisible(NULL))
  pid <- Sys.getpid()
  path <- file.path(dir, sprintf("rnaparallel-%d.tsv", pid))
  line <- sprintf("%.0f\t%s\t%d\t%s\n", as.numeric(Sys.time()), stage, chunk, event)
  # append = TRUE, one write per line: a worker that dies mid-chunk leaves a "start" with no
  # matching "done", which is itself useful (rnaparallel_progress() reports it as stalled)
  # rather than losing the row a buffered/batched write would risk on a killed process.
  try(cat(line, file = path, append = TRUE), silent = TRUE)
  invisible(NULL)
}

#' Summarise chunk progress from a `combat.progress.dir`, with an ETA
#'
#' Reads every `rnaparallel-*.tsv` file in `dir`, pairs each chunk's start/done rows, and
#' reports completed chunks, a mean seconds-per-chunk from the ones that finished, and a
#' projected finish time. Meant to be called from a SEPARATE R session while the run that is
#' writing the files is still blocked inside its parallel call: that is the whole point of
#' writing to a file rather than a console the blocked session cannot flush anyway.
#'
#' @param dir Directory passed as `options(combat.progress.dir = ...)` in the running session.
#' @return Invisibly, a list with `done`, `started`, `stalled` (started, never finished) and
#'   `eta` (a `POSIXct`, or `NA` if fewer than two chunks have finished). Also prints one line.
#' @examples
#' \dontrun{
#' # in the running session:
#' options(combat.progress.dir = "/tmp/rnaparallel-progress")
#' ComBat_seq_parallel(counts, batch = batch, group = group, workers = 16L)
#'
#' # from a second session, while the first is still running:
#' rnaparallel_progress("/tmp/rnaparallel-progress")
#' }
#' @export
rnaparallel_progress <- function(dir) {
  files <- list.files(dir, pattern = "^rnaparallel-.*\\.tsv$", full.names = TRUE)
  if (!length(files)) {
    message("no rnaparallel-*.tsv files in ", dir, " yet")
    return(invisible(list(done = 0L, started = 0L, stalled = 0L, eta = as.POSIXct(NA))))
  }
  rows <- do.call(rbind, lapply(files, function(f) {
    tryCatch(
      utils::read.delim(f, header = FALSE, sep = "\t",
                        col.names = c("ts", "stage", "chunk", "event"),
                        colClasses = c("numeric", "character", "integer", "character")),
      error = function(e) NULL)
  }))
  if (is.null(rows) || !nrow(rows)) {
    message("no readable progress rows in ", dir, " yet")
    return(invisible(list(done = 0L, started = 0L, stalled = 0L, eta = as.POSIXct(NA))))
  }
  starts <- rows[rows$event == "start", ]
  dones  <- rows[rows$event == "done", ]
  key <- function(d) paste(d$stage, d$chunk)
  done_key <- key(dones)
  stalled <- starts[!(key(starts) %in% done_key), , drop = FALSE]

  n_done <- nrow(dones)
  n_started <- nrow(starts)
  secs_per_chunk <- NA_real_
  eta <- as.POSIXct(NA)
  if (n_done >= 2L) {
    # matched by (stage, chunk): a chunk's own start row, not the run's earliest start, so a
    # straggler chunk started late does not get credited with an unfairly long duration.
    m <- merge(starts[c("stage", "chunk", "ts")], dones[c("stage", "chunk", "ts")],
              by = c("stage", "chunk"), suffixes = c("_start", "_done"))
    durs <- m$ts_done - m$ts_start
    secs_per_chunk <- mean(durs[durs >= 0])
    remaining <- n_started - n_done
    if (!is.na(secs_per_chunk) && remaining > 0L) {
      eta <- Sys.time() + secs_per_chunk * remaining
    }
  }
  msg <- sprintf("%d done, %d started, %d stalled%s%s",
                n_done, n_started, nrow(stalled),
                if (!is.na(secs_per_chunk)) sprintf(", %s/chunk", rp_secs(secs_per_chunk)) else "",
                if (!is.na(eta)) sprintf(", ETA %s", format(eta, "%H:%M")) else "")
  message(msg)
  invisible(list(done = n_done, started = n_started, stalled = nrow(stalled), eta = eta))
}


# ---- the step wrapper --------------------------------------------------------

#' Default label for a companion call
#'
#' `ComBat-seq 18,270 x 1,500`. Enough to tell two calls apart in a loop over cohorts without
#' the caller naming every one; `label` overrides it when they want to.
#' @noRd
rp_label <- function(what, x) {
  d <- tryCatch(dim(x), error = function(e) NULL)
  if (length(d) != 2L) return(what)
  sprintf("%s %s x %s", what, format(d[1L], big.mark = ","), format(d[2L], big.mark = ","))
}

#' Begin a timed, optionally quiet, progress-ticking companion step
#'
#' Returns a handle for `rp_step_end()`, or NULL when timing, quiet, and progress are all off
#' (progress defaults on, so this is the explicit `combat.progress = FALSE` case). Registering
#' the teardown with `on.exit()` in the caller is what makes this exception-safe: the sink
#' unwinds and the elapsed line still prints when the original throws, so a failed run reports
#' where it failed rather than vanishing.
#' @noRd
rp_step_begin <- function(label, what, x, backend, workers) {
  timing <- rp_opt_flag("combat.timing")
  quiet  <- rp_opt_flag("combat.quiet")
  progress <- rp_opt_flag("combat.progress", default = TRUE)
  if (!timing && !quiet && !progress) return(NULL)
  # NOT reentrant, on purpose. calcNormFactors_parallel on a DGEList reaches the original's
  # DGEList method, which calls the companion again on the counts matrix, so one user-facing
  # call is two nested calls here and printed itself twice. Only the outermost reports, and
  # only it resets the counters, so the inner dispatches are still attributed to it.
  if (rp_or0(.rp_dispatch$depth) > 0L) {
    .rp_dispatch$depth <- .rp_dispatch$depth + 1L
    return(structure(list(nested = TRUE), class = "rp_step"))
  }
  # Build the handle FIRST and claim the depth last. Opening the sink can fail, and a depth
  # raised before a throw would never be lowered -- every later call would read as nested and
  # print nothing, which is a silent loss of exactly the output this exists to produce.
  h <- list(t0 = proc.time()[["elapsed"]],
            label = if (is.null(label)) rp_label(what, x) else as.character(label)[1L],
            backend = if (is.function(backend)) "custom" else as.character(backend)[1L],
            workers = workers,
            timing = timing,
            con = if (quiet) rp_quiet_begin() else NULL)
  rp_count_reset()
  .rp_dispatch$progress_label <- h$label
  .rp_dispatch$depth <- 1L
  h
}

#' @noRd
rp_step_end <- function(h) {
  if (is.null(h)) return(invisible(NULL))
  .rp_dispatch$depth <- max(0L, rp_or0(.rp_dispatch$depth) - 1L)
  if (isTRUE(h$nested)) return(invisible(NULL))
  rp_progress_done()
  if (!is.null(h$con)) rp_quiet_end(h$con)
  if (!isTRUE(h$timing)) return(invisible(NULL))
  secs <- proc.time()[["elapsed"]] - h$t0
  if (secs < rp_opt_secs("combat.timing.min")) return(invisible(NULL))

  par <- rp_or0(.rp_dispatch$par)
  ser <- rp_or0(.rp_dispatch$ser)
  engine <- if (par > 0L) sprintf("%s x%d", h$backend, h$workers) else "serial"
  note <- if (par > 0L && ser > 0L) sprintf("  %d par / %d gated", par, ser)
          else if (par == 0L && ser > 0L) sprintf("  %d gated", ser)
          else ""
  fb <- .rp_dispatch$fallback
  if (length(fb)) note <- paste0(note, sprintf("  [%s stood down]", paste(fb, collapse = ", ")))

  message(sprintf("  %-34s %-16s %8s%s",
                  substr(h$label, 1L, 34L), engine, rp_secs(secs), note))
  invisible(NULL)
}

#' Elapsed time a human reads at a glance
#' @noRd
rp_secs <- function(s) {
  if (s >= 3600) sprintf("%.1fh", s / 3600)
  else if (s >= 60) sprintf("%.1fm", s / 60)
  else sprintf("%.1fs", s)
}

#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
