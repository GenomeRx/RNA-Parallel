#!/usr/bin/env Rscript
## Checks the claims named in the numbered sections below, in README.md, REFERENCE.md, DESCRIPTION and NEWS.md, against R/, NAMESPACE and the rendered reports in docs/; run it from the repository root before a release, like tools/exactproof.R.
root <- "."
if (!dir.exists(file.path(root, "R"))) stop("run this from the repository root", call. = FALSE)

## The source tree is loaded, not the installed copy, so every check that calls into the package or reads its exports tests the release candidate.
suppressMessages(pkgload::load_all(root, quiet = TRUE))
rd  <- readLines(file.path(root, "README.md"), warn = FALSE)
rdt <- paste(rd, collapse = "\n")
rf  <- readLines(file.path(root, "REFERENCE.md"), warn = FALSE)
docs <- paste(c(rd, rf), collapse = "\n")
## read every rendered report, not only the macOS one. A claim measured on Linux is verifiable
## only against the report that produced it, and reading one file failed four true claims.
htmlf <- list.files(file.path(root, "docs"), pattern = "\\.html$", full.names = TRUE)
if (!length(htmlf)) stop("docs/ holds no rendered report; render one first", call. = FALSE)
html <- lapply(setNames(htmlf, basename(htmlf)), readLines, warn = FALSE)

pass <- 0L; fail <- 0L
chk <- function(ok, what, detail = "") {
  if (isTRUE(ok)) pass <<- pass + 1L
  else { fail <<- fail + 1L; cat(sprintf("  *** %s%s\n", what, if (nzchar(detail)) paste0(": ", detail) else "")) }
}

## 1a. Each README speedup table cell must match a row of its platform's report with the same companion, the same worker count and the same speedup rounded to two decimals.
reports <- c(macOS = "index.html", Linux = "linux.html", Windows = "windows.html")
cells <- function(l) trimws(strsplit(sub("^\\|", "", sub("\\|\\s*$", "", l)), "|", fixed = TRUE)[[1]])
report_rows <- function(l) {
  td <- grepl("^\\s*<td", l)
  split(sub("^\\s*<td[^>]*>(.*)</td>\\s*$", "\\1", l[td]), cumsum(grepl("^\\s*<tr", l))[td])
}
rows_by_report <- lapply(html, report_rows)
pipe_lines <- which(startsWith(rd, "|"))
hdr <- pipe_lines[vapply(rd[pipe_lines], function(l) any(cells(l) %in% names(reports)), logical(1))][1]
speedup_rows <- list()
if (is.na(hdr)) {
  chk(FALSE, "README has no speedup table with a macOS, Linux or Windows column")
} else {
  hcells <- cells(rd[hdr])
  body <- rd[seq(hdr + 2L, length(rd))]
  body <- body[seq_len(match(FALSE, startsWith(body, "|"), nomatch = length(body) + 1L) - 1L)]
  speedup_rows <- lapply(body, cells)
  chk(length(speedup_rows) > 0L, "the README speedup table has no rows")
  for (r in speedup_rows) {
    stem <- sub("_parallel\\(\\)$", "", gsub("`", "", r[1L]))
    for (j in which(hcells %in% names(reports))) {
      plat <- hcells[j]
      cell <- gsub("\\*", "", r[j])
      what <- sprintf("README %s %s cell '%s'", stem, plat, cell)
      if (!grepl("^[0-9]+\\.[0-9]+x @ [0-9]+w$", cell)) {
        chk(FALSE, paste(what, "does not read N.NNx @ Kw"))
        next
      }
      val <- round(as.numeric(sub("x @.*$", "", cell)), 2L)
      k <- sub("^.*@ ([0-9]+)w$", "\\1", cell)
      rows <- rows_by_report[[reports[[plat]]]]
      if (is.null(rows)) {
        chk(FALSE, paste(what, "has no report"), file.path("docs", reports[[plat]]))
        next
      }
      hit <- vapply(rows, function(x) {
        grepl(stem, x[1L], fixed = TRUE) &&
          (k %in% x[-1L] || grepl(paste0("\\b", k, " workers\\b"), x[1L])) &&
          any(round(suppressWarnings(as.numeric(sub("x$", "", x[-1L]))), 2L) == val, na.rm = TRUE)
      }, logical(1))
      chk(any(hit), paste(what, "matches no row of", reports[[plat]]))
    }
  }
}
documented <- sub("_parallel\\(\\)$", "", gsub("`", "", vapply(speedup_rows, `[`, "", 1L)))
for (f in grep("_parallel$", getNamespaceExports("rnaparallel"), value = TRUE))
  chk(sub("_parallel$", "", f) %in% documented, sprintf("the README speedup table has no row for %s()", f))

## 1b. Any other N.NNx in README.md or REFERENCE.md must appear in a report or an R/ comment, which catches a number found nowhere but not one tied to the wrong platform or worker count.
hl <- unlist(html, use.names = FALSE)
nums <- suppressWarnings(as.numeric(unlist(regmatches(hl, gregexpr("[0-9]+\\.[0-9]+", hl)))))
seen <- unique(round(nums[is.finite(nums)], 2L))
src0 <- unlist(lapply(list.files(file.path(root, "R"), full.names = TRUE), readLines, warn = FALSE))
elsewhere <- round(as.numeric(sub("x$", "", unlist(regmatches(src0, gregexpr("[0-9]+\\.[0-9]+x", src0))))), 2L)
in_table <- if (is.na(hdr)) integer() else hdr + seq_along(speedup_rows) + 1L
prose <- c(README = paste(rd[setdiff(seq_along(rd), in_table)], collapse = "\n"),
           REFERENCE.md = paste(rf, collapse = "\n"))
for (nm in names(prose)) {
  for (c in unique(regmatches(prose[[nm]], gregexpr("[0-9]+\\.[0-9]+x", prose[[nm]]))[[1]])) {
    n <- round(as.numeric(sub("x$", "", c)), 2L)
    chk(n %in% seen || n %in% elsewhere,
        sprintf("%s claims %s, which is in neither a rendered report nor R/", nm, c))
  }
}

## 2. worker default -------------------------------------------------------------------
## every entry point must defer to the same resolver, or the README's one sentence about the
## default is true of some of them and quietly false of the rest
fns <- grep("_parallel$", getNamespaceExports("rnaparallel"), value = TRUE)
defs <- vapply(fns, function(f) is.null(formals(get(f))$workers), logical(1))
chk(all(defs), "workers no longer defaults to NULL everywhere",
    paste(names(defs)[!defs], collapse = "/"))
want <- max(1L, min(8L, max(1L, parallel::detectCores()) - 2L))
## Without fork() the pick is capped again at the performance-core count, because an efficiency
## core that has to be handed a serialized copy stops paying for itself. The check follows the
## code rather than asserting one formula on every platform, which would fail on Windows for a
## default that is deliberately lower there.
if (identical(.Platform$OS.type, "windows")) {
  want <- max(1L, min(want, rnaparallel:::rp_perf_cores()))
}
chk(identical(rnaparallel:::combat_default_workers(NULL), want),
    "the resolved worker default is not min(8, detectCores() - 2), capped at performance cores")
chk(identical(rnaparallel:::combat_default_workers(3L), 3L),
    "an explicit workers value is not passed through untouched")
chk(grepl("min(8, detectCores() - 2)", rdt, fixed = TRUE),
    "README no longer states the workers default")

## 3. backend list ---------------------------------------------------------------------
# require the quoted form the backend list uses. Accepting the name anywhere in the file let a
# backend be dropped from the list and still pass, because several are also named in prose.
for (b in combat_backends())
  chk(grepl(paste0("`\"", b, "\"`"), rdt, fixed = FALSE),
      sprintf("backend %s is offered by the code but not listed in README", b))
named <- regmatches(rdt, gregexpr('`"[a-zA-Z]+"`', rdt))[[1]]
named <- gsub('[`"]', "", named)
for (b in setdiff(named, c(combat_backends(), "rob", "ls", "robust", "topbottom")))
  chk(FALSE, sprintf("README names backend %s which combat_backends() does not offer", b))

## 4. options: every combat.* the code reads is named in README.md or REFERENCE.md, and the docs name nothing the code does not read
src <- unlist(lapply(list.files(file.path(root, "R"), full.names = TRUE), readLines, warn = FALSE))
# require a real trailing segment: a bare "combat.min." prefix appears inside message text
opt_re <- 'combat\\.[a-z]+(\\.[a-z]+)+'
in_code <- sort(unique(regmatches(paste(src, collapse = "\n"),
                gregexpr(opt_re, paste(src, collapse = "\n")))[[1]]))
in_doc <- sort(unique(regmatches(docs, gregexpr(opt_re, docs))[[1]]))
for (o in in_doc) {
  # a retired name is allowed where the README is explicitly documenting the rename
  retired <- grepl(paste0(o, "` became|", o, "` is no longer|", o, "`, renamed"), docs)
  chk(o %in% in_code || retired,
      sprintf("README or REFERENCE.md names option %s which no longer exists in R/", o))
}
for (o in in_code)
  chk(o %in% in_doc, sprintf("R/ reads option %s, which neither README nor REFERENCE.md names", o))
gates <- sum(grepl("^combat\\.min\\.", in_code))
counts <- c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve")
said <- regmatches(docs, gregexpr(paste0("\\b(", paste(counts, collapse = "|"), ") (options set the|`combat[.]min)"),
                                  docs, ignore.case = TRUE))[[1]]
for (s in said)
  chk(match(tolower(sub(" .*$", "", s)), counts) == gates,
      sprintf("the docs say '%s' but the code reads %d combat.min.* gates", gsub("\\s+", " ", s), gates))

## 5. every package function README.md names must be exported, and every one REFERENCE.md names must exist, since REFERENCE.md may describe internals
ex <- getNamespaceExports("rnaparallel")
fn_re <- "`([a-zA-Z_]+_parallel|combat_[a-z_]+|rnaparallel_[a-z_]+|rp_[a-z_]+)\\(\\)`"
for (f in unique(regmatches(rdt, gregexpr(fn_re, rdt))[[1]])) {
  f <- gsub("[`()]", "", f)
  chk(f %in% ex, sprintf("README names %s() which is not exported", f))
}
rft <- paste(rf, collapse = "\n")
for (f in unique(regmatches(rft, gregexpr(fn_re, rft))[[1]])) {
  f <- gsub("[`()]", "", f)
  chk(exists(f, envir = asNamespace("rnaparallel"), inherits = FALSE),
      sprintf("REFERENCE.md names %s() which does not exist in R/", f))
}

## 6. no brittle assertion counts ------------------------------------------------------
## An exact count drifts the moment a test is added and no gate short of running the suite can
## see it. It went stale twice. State the shape of the suite, not its arithmetic.
## "Over 350 assertions" is fine, it does not drift. A bare count preceding the word, or a
## count paired with a second one for R CMD check, is the brittle form.
brittle <- grepl("(^|[^a-zA-Z]) ?[0-9]{3} assertions under|assertions under [^,]+, [0-9]{3}", rdt)
chk(!brittle, "README states an exact assertion count, which drifts on every added test")

## 7. version consistency ---------------------------------------------------------------
## DESCRIPTION, not packageVersion(): this gate runs from the source tree before a release,
## and the installed copy is incidental to it. Reading the installed version made the check
## fail whenever source was ahead of the library, which is the normal state while preparing
## a release, and was unfixable at all while an install freeze was in force.
ver <- as.character(read.dcf(file.path(root, "DESCRIPTION"), fields = "Version")[[1L]])
news <- readLines(file.path(root, "NEWS.md"), warn = FALSE)[1]
chk(grepl(ver, news, fixed = TRUE),
    "NEWS.md does not open with the DESCRIPTION version", sprintf("DESCRIPTION %s, NEWS '%s'", ver, news))

## 8. the README's version badge --------------------------------------------------------
## A badge is the first thing a reader sees and the first thing to go stale, because nothing
## about it breaks when it is wrong. It is the version people will quote back, so it is checked
## against DESCRIPTION exactly as NEWS.md is.
badge <- regmatches(rdt, regexpr("badge/version-[0-9][0-9.]*-", rdt))
chk(length(badge) == 1L && identical(gsub("^badge/version-|-$", "", badge), ver),
    "the README version badge does not match DESCRIPTION",
    sprintf("badge %s, DESCRIPTION %s",
            if (length(badge)) sQuote(gsub("^badge/version-|-$", "", badge)) else "<none found>",
            ver))

## 9. gated symbols ---------------------------------------------------------------------
## A reachability gate is the one guard against a rebind that silently stops reaching, so a doc
## that names the wrong set is worse than one that names none. This drifted already: the gate
## was cut from five symbols to three, and NEWS went on claiming five for two releases, because
## nothing compared the sentence to the list.
newst <- paste(readLines(file.path(root, "NEWS.md"), warn = FALSE), collapse = "\n")
current <- sub("\n# rnaparallel.*$", "", newst)          # this release's section only
src <- paste(unlist(lapply(list.files(file.path(root, "R"), full.names = TRUE),
                           readLines, warn = FALSE)), collapse = "\n")
gated <- unique(unlist(regmatches(src, gregexpr('rebound (<-|=) c\\([^)]*\\)', src))))
gated <- unique(unlist(regmatches(gated, gregexpr('"[A-Za-z_.]+"', gated))))
gated <- gsub('"', "", gated)
chk(length(gated) > 0L, "no rebound list found in R/")
## The direction that matters is docs -> code: NEWS naming a symbol as gated when it is not is
## a false assurance about the one guard that catches a silently-serial rebind. The reverse
## (a gated symbol NEWS does not mention) is normal, since NEWS describes what changed.
claimed <- unlist(regmatches(current, gregexpr("gate covers[^.]*", current)))
claimed <- gsub("`", "", unlist(regmatches(claimed, gregexpr("`[A-Za-z_.]+`", claimed))))
for (g in unique(claimed)) {
  chk(g %in% gated,
      sprintf("NEWS says the reachability gate covers `%s`, which is in no rebound list in R/", g))
}
if (!length(claimed)) cat("  note: this release's NEWS section has no 'gate covers' sentence, so section 9 compared nothing\n")

## 10. DESCRIPTION's Description names the original behind every speedup table row and every backend
desc <- read.dcf(file.path(root, "DESCRIPTION"), fields = "Description")[[1L]]
flat <- function(s) gsub("[-_[:space:]]", "", tolower(s))
runs_col <- if (is.na(hdr)) NA_integer_ else match("runs", hcells)
chk(!is.na(runs_col), "the README speedup table has no 'runs' column naming each original")
for (r in if (is.na(runs_col)) list() else speedup_rows) {
  orig <- sub("^.*::", "", gsub("`", "", r[runs_col]))
  chk(grepl(flat(orig), flat(desc), fixed = TRUE),
      sprintf("DESCRIPTION's Description does not name %s, the original behind %s", orig, r[1L]))
}
for (b in combat_backends())
  chk(grepl(b, desc, fixed = TRUE), sprintf("DESCRIPTION's Description does not name backend %s", b))

cat(sprintf("\n==== doccheck: %d checks, %d passed, %d failed ====\n", pass + fail, pass, fail))
quit(status = if (fail == 0L) 0L else 1L)
