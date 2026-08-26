## zzz.R -- build identity, recorded at load

.rp_build <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  # Stamped once, cheaply. Deliberately NOT a force() sweep over the namespace: that costs
  # every consumer on every library() call to detect a condition only a reinstall creates,
  # and it cannot catch the case that actually bites, which is a reinstall landing AFTER
  # library() has already run.
  # suppressWarnings because pkgload::load_all() has no installed copy to describe, and a
  # warning on every source load would be noise for a stamp that is simply NA there
  .rp_build$stamp <- tryCatch(
    suppressWarnings(utils::packageDescription(pkgname, lib.loc = libname)$Built),
    error = function(e) NA_character_)
  invisible(NULL)
}

#' Has this package been reinstalled underneath the running session?
#'
#' R lazy-loads function bodies by byte offset into `rnaparallel.rdb`. Reinstalling the
#' package rewrites that file, so a session that was already running gets
#' `lazy-load database ... is corrupt` the first time it touches a symbol it had not yet
#' forced. Inside a forked worker that surfaces as a failed dispatch, and a caller with a
#' fallback path can quietly proceed on uncorrected data: three notebook runs were lost to
#' exactly this before anyone realised the package, not the data, was the problem.
#'
#' Nothing here prevents it. The fix is to restart R after installing. This only lets a long
#' run notice, by comparing the build stamp this session loaded against the one on disk now.
#'
#' Call it defensively. If the reinstall crossed the release that introduced this function, the
#' loaded namespace does not contain it and the call raises `could not find function` instead of
#' returning `TRUE` -- the checker defeated by the very condition it detects. That error IS the
#' stale answer, so treat it as one:
#'
#' ```r
#' if (isTRUE(tryCatch(rnaparallel_stale(), error = function(e) TRUE))) {
#'   stop("rnaparallel was reinstalled under this session; restart R.")
#' }
#' ```
#'
#' @return `TRUE` when the installed build differs from the one this session loaded, `FALSE`
#'   when they agree, and `NA` when the installed build cannot be read at all.
#' @examples
#' if (isTRUE(tryCatch(rnaparallel_stale(), error = function(e) TRUE))) {
#'   message("rnaparallel was reinstalled under this session; restart R.")
#' }
#' @export
rnaparallel_stale <- function() {
  # Version first, because it needs NOTHING recorded at load time. `packageVersion()` reads
  # DESCRIPTION from disk and so reports the NEW version while the running namespace is the old
  # one, which is why it is useless on its own and decisive when compared against the namespace.
  loaded <- tryCatch(as.character(getNamespaceVersion("rnaparallel")),
                     error = function(e) NA_character_)
  disk <- tryCatch(as.character(utils::packageVersion("rnaparallel")),
                   error = function(e) NA_character_)
  if (!is.na(loaded) && !is.na(disk) && !identical(loaded, disk)) return(TRUE)

  # Same version can still be a different build: reinstalling a version over itself rewrites the
  # lazy-load database just as surely, and the version comparison cannot see that.
  now <- tryCatch(utils::packageDescription("rnaparallel")$Built,
                  error = function(e) NA_character_)
  if (is.null(now) || is.na(now) || is.na(.rp_build$stamp %||% NA_character_)) return(NA)
  !identical(.rp_build$stamp, now)
}
