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

#' Read a single non-negative numeric option, refusing a garbage value
#'
#' `what` only changes the wording: `combat.timing.min` reads naturally as "a number of
#' seconds", `combat.mem.divergence` does not (it errored with that exact wording for a
#' value that is a fraction, not a duration, before this had a `what` argument at all).
#' `rp_opt_secs()`/`rp_opt_num()` below are the two names every call site already uses;
#' both are one-line wrappers so neither existing caller needed to change.
#' @noRd
rp_opt_number <- function(name, default, what) {
  v <- suppressWarnings(as.numeric(getOption(name, default)))
  if (length(v) != 1L || is.na(v) || v < 0) {
    stop("`", name, "` must be a single non-negative ", what, "; got ",
         deparse(getOption(name)), call. = FALSE)
  }
  v
}

#' @noRd
rp_opt_secs <- function(name, default = 0) rp_opt_number(name, default, "number of seconds")

#' @noRd
rp_opt_num <- function(name, default = 0) rp_opt_number(name, default, "number")


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
  # Only the MASTER should ever draw this tick. The two size/nested-worker gates in
  # combat_parallel_lapply() that fall through to a plain lapply() call rp_count(FALSE),
  # which calls this, from INSIDE a forked child on a recursive dispatch (ComBat-seq's
  # tagwise loop across batches, dispatched again per-batch inside an outer worker) --
  # printing that worker's own stale "N dispatched" to the shared stderr, contradicting
  # whatever the real master is drawing at the same moment. RNAPARALLEL_IN_WORKER is set
  # for the whole lifetime of a worker process (see f_tagged in helper_seq_parallel.R), so
  # this is a correct "am I the master" check, not a proxy for it.
  if (nzchar(Sys.getenv("RNAPARALLEL_IN_WORKER"))) return(invisible(NULL))
  # A reporter is rendering the bar on this same line with its own carriage returns. Two writers
  # on one line interleave into unreadable output, so the counter stands down while it runs.
  if (isTRUE(.rp_dispatch$reporter)) return(invisible(NULL))
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
#
# One file per worker PID, and the file is APPEND-ONLY for the life of that PID: a worker
# pool reused across multiple companion calls in one long pipeline script (five companions,
# one shared `future::multisession` plan, `combat.progress.dir` set once at the top) keeps
# writing to the SAME files across every call, since the PID does not change. That is
# correct within any one call (rp_progress_summarise() needs every row to tell "started,
# never finished" from "not started yet"), but nothing here truncates a file BETWEEN calls,
# so the recommendation for a long pipeline is a fresh subdirectory per stage
# (`options(combat.progress.dir = file.path(base, "lmfit"))`, then `file.path(base, "combat")`)
# rather than one directory for the whole script. Truncating mid-run is not the fix: that
# would corrupt exactly the "stalled" and "done" counts this mechanism exists to report.

#' @noRd
rp_progress_dir <- function() {
  d <- getOption("combat.progress.dir", NA_character_)
  if (is.na(d) || !nzchar(d)) {
    # Default on: with no dir the run is silent for its whole duration, which is the state this
    # package exists to fix. Under tempdir(), so it never lands in a repo.
    d <- file.path(tempdir(), "rnaparallel-progress")
  }
  if (!dir.exists(d)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(d)) return(NULL)
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
#' A single continuously-updating bar DURING one blocking call is not something the running
#' session can print: once it calls into `mclapply`/`future`/`BiocParallel`/`foreach` it is
#' synchronously waiting and cannot redraw a console until the call returns, which is the whole
#' reason this writes to files instead. `watch = TRUE` gets the live-bar behaviour anyway, from
#' the side that CAN keep drawing: this function's own process, polling the files and
#' redrawing a real `[#####-----] 47%` bar with the current stage name, once a second, until
#' every chunk in the last dispatch is done.
#'
#' @param dir Directory passed as `options(combat.progress.dir = ...)` in the running session.
#' @param watch If `TRUE`, poll and redraw a live bar every `interval` seconds instead of
#'   returning once. Meant for a SEPARATE terminal/session next to the one running the actual
#'   computation; stop it with Ctrl-C or `interval` reaching a stall (see below).
#' @param interval Seconds between redraws in watch mode.
#' @param stall_after Seconds with no new "done" OR "start" row before watch mode gives up
#'   and returns, so a finished or crashed run does not poll forever with nobody watching.
#'   Default 10 minutes: long enough to survive one very slow chunk, short enough to actually
#'   stop.
#' @return Invisibly, a list with `done`, `started`, `stalled` (started, never finished) and
#'   `eta` (a `POSIXct`, or `NA` if fewer than two chunks have finished) -- on every return
#'   path, including a stall: the last summary read is still the most useful answer, and
#'   was already computed, so there is no reason to throw it away for `NULL`.
#' @examples
#' \dontrun{
#' # in the running session:
#' options(combat.progress.dir = "/tmp/rnaparallel-progress")
#' ComBat_seq_parallel(counts, batch = batch, group = group, workers = 16L)
#'
#' # from a second session, while the first is still running:
#' rnaparallel_progress("/tmp/rnaparallel-progress")
#'
#' # or, for a live-updating bar in that second session:
#' rnaparallel_progress("/tmp/rnaparallel-progress", watch = TRUE)
#' }
#' @export
rnaparallel_progress <- function(dir, watch = FALSE, interval = 1, stall_after = 600) {
  if (isTRUE(watch)) {
    # The package refuses garbage everywhere else (rp_opt_flag, rp_opt_secs, every
    # combat.min.* gate); this was the one public entry point that did not. interval <= 0
    # either busy-polls (0) or reaches Sys.sleep() as an outright error (negative), and both
    # are silent until they happen rather than refused at the door like everything else.
    if (!(is.numeric(interval) && length(interval) == 1L && !is.na(interval) && interval > 0)) {
      stop("`interval` must be a single positive number of seconds; got ",
           deparse(interval), call. = FALSE)
    }
    if (!(is.numeric(stall_after) && length(stall_after) == 1L && !is.na(stall_after) &&
          stall_after > 0)) {
      stop("`stall_after` must be a single positive number of seconds; got ",
           deparse(stall_after), call. = FALSE)
    }
    return(rp_progress_watch(dir, interval, stall_after))
  }
  rp_progress_once(dir)
}

#' Read every worker's progress file, optionally tailing rather than re-reading whole
#'
#' Each worker only ever APPENDS a line (`rp_progress_file_write()`, `append = TRUE`,
#' one `cat()` per chunk); no line already on disk is ever rewritten. That makes byte-offset
#' tailing safe: bytes read on a previous call can never have changed, so re-reading them
#' again buys nothing. Pass `cache` (a mutable environment, one entry per file path holding
#' the last-seen file size and the parsed rows through that point) to skip re-reading and
#' re-parsing everything already seen; omit it for a one-shot call, where there is only ever
#' one read and caching has nothing to save.
#'
#' Matters specifically in `watch = TRUE` mode: a redraw every `interval` seconds over a
#' multi-hour run means the same worker files get read from byte 0 every single poll without
#' this, so the read cost grows with elapsed time even though each poll only cares about the
#' handful of lines written since the last one.
#' @noRd
rp_progress_read <- function(dir, cache = NULL) {
  files <- list.files(dir, pattern = "^rnaparallel-.*\\.tsv$", full.names = TRUE)
  if (!length(files)) return(NULL)
  rows <- do.call(rbind, lapply(files, rp_progress_read_file, cache = cache))
  if (is.null(rows) || !nrow(rows)) return(NULL)
  rows
}

#' @noRd
rp_progress_parse_lines <- function(lines) {
  if (!length(lines)) return(NULL)
  tryCatch(
    utils::read.delim(text = paste(lines, collapse = "\n"), header = FALSE, sep = "\t",
                      col.names = c("ts", "stage", "chunk", "event"),
                      colClasses = c("numeric", "character", "integer", "character")),
    error = function(e) NULL)
}

#' @noRd
rp_progress_read_file <- function(path, cache) {
  if (is.null(cache)) {
    lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
    return(rp_progress_parse_lines(lines))
  }
  sz <- tryCatch(file.size(path), error = function(e) NA_real_)
  prev <- cache[[path]]
  # File shrank (rare: cleared between polls) means the byte offset from last time no
  # longer means anything on this file, so start over rather than seek past the end.
  if (!is.null(prev) && !is.na(sz) && sz >= prev$size) {
    if (sz == prev$size) return(prev$rows)   # nothing new: skip the read entirely
    con <- tryCatch(file(path, open = "rb"), error = function(e) NULL)
    if (is.null(con)) return(prev$rows)
    on.exit(close(con), add = TRUE)
    if (prev$size > 0) seek(con, where = prev$size, origin = "start")
    new_lines <- tryCatch(readLines(con, warn = FALSE), error = function(e) character())
    new_rows <- rp_progress_parse_lines(new_lines)
    combined <- if (is.null(new_rows)) prev$rows
                else if (is.null(prev$rows)) new_rows
                else rbind(prev$rows, new_rows)
    cache[[path]] <- list(size = sz, rows = combined)
    return(combined)
  }
  lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
  rows <- rp_progress_parse_lines(lines)
  cache[[path]] <- list(size = if (is.na(sz)) 0 else sz, rows = rows)
  rows
}

#' @noRd
rp_progress_summarise <- function(rows) {
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
  # The stage shown is whichever one has the most recent activity, so a multi-stage run (five
  # companions in one script, one progress.dir) shows the stage actually running right now
  # rather than the alphabetically first one.
  cur_stage <- if (nrow(rows)) rows$stage[which.max(rows$ts)] else NA_character_
  list(done = n_done, started = n_started, stalled = nrow(stalled), eta = eta,
      secs_per_chunk = secs_per_chunk, stage = cur_stage)
}

#' @noRd
rp_progress_once <- function(dir) {
  rows <- rp_progress_read(dir)
  if (is.null(rows)) {
    message("no rnaparallel-*.tsv files in ", dir, " yet")
    return(invisible(list(done = 0L, started = 0L, stalled = 0L, eta = as.POSIXct(NA))))
  }
  s <- rp_progress_summarise(rows)
  msg <- sprintf("%d done, %d started, %d stalled%s%s",
                s$done, s$started, s$stalled,
                if (!is.na(s$secs_per_chunk)) sprintf(", %s/chunk", rp_secs(s$secs_per_chunk)) else "",
                if (!is.na(s$eta)) sprintf(", ETA %s", format(s$eta, "%H:%M")) else "")
  message(msg)
  invisible(s[c("done", "started", "stalled", "eta")])
}

#' A live `|====------|` bar, this package's own format, drawn in the WATCHING process, not
#' the blocked one
#'
#' This is the process that can actually keep redrawing: the running session is synchronously
#' blocked inside its parallel call and cannot. Stops on its own once `started` chunks stop
#' growing for `stall_after` seconds (the run finished, or nobody is writing to `dir` at all)
#' so a call left running does not poll an abandoned directory forever.
#'
#' Single line, overwritten in place with a carriage return, same mechanism as the console
#' tick: a redraw needs a live process polling a SEPARATE process's progress files, which is
#' a different problem from printing your own progress once. A two-line redraw would need
#' ANSI cursor-up escapes that not every terminal honours identically, so one line keeps the
#' `|===---|` look without that risk.
#' Render the progress bar from a process that is not blocked
#'
#' The console tick is emitted by the master and silenced in workers on purpose, since sixteen
#' processes interleaving carriage returns on one line is unreadable. That leaves nobody able to
#' report during a dispatch: the master is inside mclapply until the last chunk returns. Forking a
#' reporter gives the bar a writer that is neither blocked nor competing, reading the same per
#' worker files rnaparallel_progress(watch=TRUE) reads.
#'
#' Started and stopped inside combat_parallel_lapply so it never reaches combat_reap's snapshot.
#' @noRd
rp_reporter_start <- function(dir) {
  if (is.null(dir) || .Platform$OS.type != "unix") return(NULL)
  if (nzchar(Sys.getenv("RNAPARALLEL_IN_WORKER"))) return(NULL)   # nested dispatch, master already has one
  if (!isTRUE(getOption("combat.fork", TRUE))) return(NULL)
  # Quarter second before the first frame, then twice a second. Long enough that a dispatch
  # finishing instantly never paints and erases, short enough that every call doing real work
  # shows a bar.
  h <- tryCatch(parallel::mcparallel({ Sys.sleep(0.25); rp_progress_watch(dir, 0.5, .Machine$integer.max) }),
                error = function(e) NULL)
  if (!is.null(h)) .rp_dispatch$reporter <- TRUE
  h
}

#' @noRd
rp_reporter_stop <- function(h) {
  .rp_dispatch$reporter <- FALSE
  if (is.null(h)) return(invisible(NULL))
  try(tools::pskill(h$pid, tools::SIGKILL), silent = TRUE)
  try(suppressWarnings(parallel::mccollect(h, wait = FALSE)), silent = TRUE)
  # Erase the bar rather than newline past it. The reporter dies mid-frame, so whatever it had
  # painted stays on the terminal and reads as belonging to whichever call prints next: a TMM bar
  # sitting under an lmFit heading is worse than no bar. rp_step_end prints the real summary.
  cat("\r", strrep(" ", 120L), "\r", sep = "", file = stderr())
  invisible(NULL)
}

#' @noRd
rp_progress_watch <- function(dir, interval, stall_after) {
  last_activity <- Sys.time()
  last_started <- -1L
  last_done <- -1L
  last_summary <- NULL
  cr <- "\r"
  width <- 50L    # this package's own bar width
  # One cache across every poll in this watch call: the whole reason to tail is to stop
  # re-reading and re-parsing lines already seen, which only pays off across REPEATED reads
  # of the same files. A fresh cache per poll would tail nothing and cost the same as none.
  cache <- new.env(parent = emptyenv())
  repeat {
    rows <- rp_progress_read(dir, cache = cache)
    if (is.null(rows)) {
      cat(cr, strrep(" ", 70L), cr, "  waiting for ", dir, " ...", sep = "")
      utils::flush.console()
    } else {
      s <- rp_progress_summarise(rows)
      last_summary <- s
      # Activity is either a chunk STARTING or a chunk FINISHING: tracking only `started`
      # (as this used to) means a run whose chunks were all dispatched early (preschedule,
      # or simply a fast dispatcher) but whose done rows are still trickling in hits the
      # stall exit mid-run, with the misleading "no new chunks" message on a run that is
      # very much still producing chunks -- just not new START rows.
      if (s$started != last_started || s$done != last_done) {
        last_activity <- Sys.time(); last_started <- s$started; last_done <- s$done
      }
      pct <- if (s$started > 0L) s$done / s$started else 0
      filled <- round(pct * width)
      bar <- paste0("|", strrep("=", filled), strrep("-", width - filled), "|")
      eta_txt <- if (!is.na(s$eta)) sprintf("  ETA %s", format(s$eta, "%H:%M")) else ""
      stage_txt <- substr(s$stage %||% "", 1L, RP_LABEL_WIDTH)
      line <- sprintf("%s %3.0f%%  %-*s %d/%d%s",
                      bar, pct * 100, RP_LABEL_WIDTH, stage_txt, s$done, s$started, eta_txt)
      cat(cr, strrep(" ", width + 60L), cr, line, sep = "")
      utils::flush.console()
      if (s$started > 0L && s$done >= s$started &&
          as.numeric(Sys.time() - last_activity, units = "secs") > 2) {
        cat("\n")
        return(invisible(s[c("done", "started", "stalled", "eta")]))
      }
    }
    if (as.numeric(Sys.time() - last_activity, units = "secs") > stall_after) {
      cat("\n  no new chunks in ", stall_after, "s, stopping watch", sep = "")
      if (!is.null(last_summary) && last_summary$stalled > 0L) {
        cat(sprintf(" (%d chunk%s never finished)", last_summary$stalled,
                    if (last_summary$stalled == 1L) "" else "s"))
      }
      cat("\n")
      return(invisible(
        if (is.null(last_summary)) list(done = 0L, started = 0L, stalled = 0L, eta = as.POSIXct(NA))
        else last_summary[c("done", "started", "stalled", "eta")]))
    }
    Sys.sleep(interval)
  }
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
#' Returns a handle for `rp_step_end()`, or NULL when timing, quiet, console progress, and file
#' progress (`combat.progress.dir`) are all off (progress defaults on, so this is the explicit
#' `combat.progress = FALSE` case with no directory set). Registering the teardown with
#' `on.exit()` in the caller is what makes this exception-safe: the sink unwinds and the
#' elapsed line still prints when the original throws, so a failed run reports where it failed
#' rather than vanishing.
#' @noRd
rp_step_begin <- function(label, what, x, backend, workers) {
  timing <- rp_opt_flag("combat.timing")
  quiet  <- rp_opt_flag("combat.quiet")
  progress <- rp_opt_flag("combat.progress", default = TRUE)
  # File progress needs the real stage label even when the console tick is off: a caller who
  # wants only combat.progress.dir (say, on a headless RStudio Server run where the console
  # tick is pointless) still needs .rp_dispatch$progress_label set below to something other
  # than the generic "dispatch" fallback, or every stage's TSV rows read identically and
  # rnaparallel_progress() cannot tell a ComBat-seq run from an lmFit run in the same dir.
  pdir <- rp_progress_dir()
  file_progress <- !is.null(pdir)
  # Created here, once, in the MASTER, rather than left to the first worker's own
  # rp_progress_file_write(): that function's `cat(..., append = TRUE)` does not create
  # missing directories and wraps everything in `try(silent = TRUE)`, so a typo'd or
  # unwritable combat.progress.dir used to fail completely silently -- every worker's write
  # dropped on the floor, and rnaparallel_progress() reported "no files yet" for the entire
  # run with no indication anything was ever wrong. Failing loudly here, once, before any
  # work starts, is strictly better than a run that "worked" with a progress feature quietly
  # doing nothing for its whole duration.
  if (file_progress) {
    ok <- dir.create(pdir, recursive = TRUE, showWarnings = FALSE) || dir.exists(pdir)
    if (!ok) {
      stop("combat.progress.dir (", pdir, ") could not be created. Check the path and ",
           "permissions, or unset the option to run without file progress.", call. = FALSE)
    }
    probe <- file.path(pdir, sprintf(".rp_write_test_%d", Sys.getpid()))
    writable <- isTRUE(tryCatch({ cat("", file = probe); file.remove(probe); TRUE },
                                error = function(e) FALSE))
    if (!writable) {
      stop("combat.progress.dir (", pdir, ") exists but is not writable. Every worker's ",
           "progress write would silently do nothing for the whole run.", call. = FALSE)
    }
  }
  if (!timing && !quiet && !progress && !file_progress) return(NULL)
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

#' Peak resident bytes this process has ever held
#'
#' VmHWM, not VmRSS: the high-water mark is what the fit needed, and it is still readable
#' after the memory has been released. `ps::ps_memory_info()$peak_wset` (Windows/macOS peak
#' working set) as a cross-platform fallback, same reasoning as `rp_mem_available()` and
#' `rp_mem_rss()` in `helper_seq_parallel.R`: without it, the peak-RSS column in the
#' `combat.timing` line never appears off Linux at all. NA when neither is available.
#' @noRd
rp_mem_peak <- function() {
  if (file.exists("/proc/self/status")) {
    l <- tryCatch(readLines("/proc/self/status"), error = function(e) character())
    m <- grep("^VmHWM:", l, value = TRUE)
    if (length(m)) return(as.numeric(gsub("\\D", "", m[1L])) * 1024)
    return(NA_real_)
  }
  if (requireNamespace("ps", quietly = TRUE)) {
    v <- tryCatch(ps::ps_memory_info()[["peak_wset"]], error = function(e) NA_real_)
    if (!is.null(v) && !is.na(v)) return(as.numeric(v))
  }
  NA_real_
}


# ---- setting R_MAX_VSIZE for people who want a second safety net -------------

# rp_mem_cap() degrades the WORKER COUNT before a fork based on live /proc readings; it is a
# preventive guess, tuned by combat.mem.divergence, and the PR that added it disclosed the
# 0.25 default missing a real case (a 4-worker fit that dirtied more per worker than assumed).
# R_MAX_VSIZE is a different, independent net: R's own vector-heap ceiling, checked on every
# allocation, that turns an overshoot into a catchable "cannot allocate vector of size X"
# error instead of leaving the kernel to SIGKILL the process with no condition raised at all.
# Setting it is standard R practice, not something this package invents, but nobody sets it
# without being told to, and the number to pick depends on a machine's own RAM, which most
# people do not have memorized in GB let alone bytes.

#' Set `R_MAX_VSIZE` to half the machine's RAM, rounded to the nearest whole tier
#'
#' Reads total RAM (`/proc/meminfo` on Linux, PowerShell's `Get-CimInstance` on Windows),
#' halves it, and rounds to the nearest of 8/16/32/64/128/256/512/1024 GB, R's own
#' vector-heap ceiling. This does NOT replace `rp_mem_cap()`: that guard degrades the worker
#' count before a fork based on a live reading of what is available right now; this sets a
#' fixed ceiling R itself enforces on every allocation, in every session, whether or not
#' this package's dispatch code is what allocated the memory. Two independent nets against
#' the same failure mode (a silent kernel SIGKILL with no R condition to catch), not one
#' superseding the other.
#'
#' Writes (or updates) the `R_MAX_VSIZE` line in the target `.Renviron`, which only takes
#' effect on the NEXT R session; `.Renviron` is read once at startup, so this cannot change
#' the limit for the session that calls it. Existing lines for other variables are left
#' untouched; only a pre-existing `R_MAX_VSIZE=` line, if any, is replaced.
#'
#' `path` defaults to `Sys.getenv("R_ENVIRON_USER", path.expand("~/.Renviron"))`, R's own
#' resolution order for the per-user file, and on Windows `path.expand("~")` resolves via
#' `USERPROFILE`, not the `HOME` environment variable: a test that tried to redirect this by
#' setting `HOME` still wrote to the real file, because `path.expand()` never consulted it.
#' Pass `path` explicitly to target anything else, which is also how to test this function
#' without touching a real `.Renviron`.
#'
#' @param fraction Fraction of total RAM to use before rounding to a tier. Default 0.5 (half),
#'   matching the standard "leave the other half for the OS, other processes, and anything not
#'   going through this package" rule of thumb.
#' @param dry_run If `TRUE` (default `FALSE`), print what would be written without touching
#'   `path`. Use this first: the computed value is worth checking before it becomes the
#'   ceiling every future R session on this machine runs under.
#' @param path `.Renviron` path to write. Defaults to R's own per-user file; override for
#'   testing or to target `Renviron.site` / a project-local `.Renviron` instead.
#' @return Invisibly, the chosen value in bytes, or `NA` if total RAM could not be read.
#' @examples
#' \dontrun{
#' rnaparallel_set_mem_limit(dry_run = TRUE)   # see what it would write, change nothing
#' rnaparallel_set_mem_limit()                 # write it; restart R for it to take effect
#' }
#' @export
rnaparallel_set_mem_limit <- function(fraction = 0.5, dry_run = FALSE,
                                      path = Sys.getenv("R_ENVIRON_USER",
                                                        path.expand("~/.Renviron"))) {
  if (!(is.numeric(fraction) && length(fraction) == 1L && !is.na(fraction) &&
        fraction > 0 && fraction <= 1)) {
    stop("`fraction` must be a single number in (0, 1]; got ", deparse(fraction), call. = FALSE)
  }
  if (!(is.character(path) && length(path) == 1L && !is.na(path) && nzchar(path))) {
    stop("`path` must be a single non-empty string; got ", deparse(path), call. = FALSE)
  }
  total <- rp_mem_total()
  if (is.na(total)) {
    message("could not read total RAM on this platform (checked /proc/meminfo and ",
            "PowerShell). Set R_MAX_VSIZE yourself in ", path, ", e.g. R_MAX_VSIZE=64Gb")
    return(invisible(NA_real_))
  }
  target <- total * fraction
  tiers_gb <- c(8, 16, 32, 64, 128, 256, 512, 1024)
  tiers_b <- tiers_gb * 2^30
  # Nearest by ratio in log space, not absolute difference: 100 GB is meant to round to 128,
  # not sit exactly between 64 and 128 by raw GB and get pulled to whichever is closer in a
  # way that ignores how differently 64->128 and 512->1024 both double.
  chosen_gb <- tiers_gb[which.min(abs(log(tiers_b / target)))]
  chosen_b <- chosen_gb * 2^30

  line <- sprintf("R_MAX_VSIZE=%dGb", chosen_gb)
  message(sprintf("total RAM ~%.0f GB, %.0f%% -> nearest tier: %s (target: %s)",
                  total / 2^30, fraction * 100, line, path))

  if (isTRUE(dry_run)) {
    message("dry_run = TRUE: nothing written. Re-run with dry_run = FALSE to set it.")
    return(invisible(chosen_b))
  }

  existing <- if (file.exists(path)) readLines(path, warn = FALSE) else character()
  # Exact-prefix match only: a line that merely mentions R_MAX_VSIZE in a comment is left
  # alone, and only a real `R_MAX_VSIZE=...` assignment is replaced.
  is_vsize <- grepl("^\\s*R_MAX_VSIZE\\s*=", existing)
  if (any(is_vsize)) {
    existing[is_vsize] <- line
  } else {
    existing <- c(existing, line)
  }
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(existing, path)
  message("wrote ", path, ". Restart R for this to take effect; .Renviron is read once ",
          "at session start.")
  invisible(chosen_b)
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

  # VmHWM is the high-water mark, so it survives the gc() that follows a big fit and reports
  # what the run actually needed rather than what it happens to hold when it finishes. Without
  # it a memory problem is invisible in the only log the run produces.
  peak <- rp_mem_peak()
  pk <- if (is.na(peak)) "" else sprintf("  peak %.0f GB", peak / 2^30)
  message(sprintf("  %-*s %-16s %8s%s%s",
                  RP_LABEL_WIDTH, substr(h$label, 1L, RP_LABEL_WIDTH), engine, rp_secs(secs), pk, note))
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

# One width, shared by every surface that prints a companion's label: the console tick
# (below), the timing line (rp_step_end), and the watch-mode bar's stage field. Longest
# realistic default label is duplicateCorrelation's own: "duplicateCorrelation 18,270 x
# 1,500" at 36 chars. Sized to 40 so the widest real default and its dimensions both always
# fit; previously the timing line used 34 and the watch bar used 28, both narrower than
# dupcor's default, so those two surfaces silently truncated it mid-number
# ("...18,270 x 1,50") while the console tick (unbounded) did not -- three different label
# widths for what is meant to be one consistent format across every parallel companion.
RP_LABEL_WIDTH <- 40L

