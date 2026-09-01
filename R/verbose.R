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


# ---- single-line progress (opt-in) --------------------------------------------

# combat.timing prints one line at the END of a call. ComBat-seq alone dispatches its hot
# paths up to 2*n_batch + 3 times per call, so on a large cohort (hundreds of batches) there
# is nothing on screen between "computing" and the final line, and a stuck run looks exactly
# like a slow one. This is the opt-in fix: one line, overwritten in place with a carriage
# return, so it never scrolls and never floods a log. Off by default, same as timing and quiet.

#' @noRd
rp_progress_tick <- function() {
  if (!rp_opt_flag("combat.progress")) return(invisible(NULL))
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
  if (!rp_opt_flag("combat.progress")) return(invisible(NULL))
  # Clear the line rather than leaving a stale count sitting there once the step's own
  # combat.timing line (if any) prints below it.
  cr <- "\r"
  cat(cr, strrep(" ", 60L), cr, sep = "", file = stderr())
  utils::flush.console()
  invisible(NULL)
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

#' Begin a timed, optionally quiet companion step
#'
#' Returns a handle for `rp_step_end()`, or NULL when neither option is on, which is the
#' default and costs one `getOption` per call. Registering the teardown with `on.exit()` in the
#' caller is what makes this exception-safe: the sink unwinds and the elapsed line still prints
#' when the original throws, so a failed run reports where it failed rather than vanishing.
#' @noRd
rp_step_begin <- function(label, what, x, backend, workers) {
  timing <- rp_opt_flag("combat.timing")
  quiet  <- rp_opt_flag("combat.quiet")
  progress <- rp_opt_flag("combat.progress")
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
