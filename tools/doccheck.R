#!/usr/bin/env Rscript
## Cross-checks every factual claim in README.md and DESCRIPTION against the code and the
## rendered report. Written because stale numbers and renamed options kept surviving hand
## review: prose that restates a fact living in code drifts the moment the code moves, and
## nothing here noticed until a reader did.
##
## Run before any release, like tools/exactproof.R.
suppressMessages({library(rnaparallel)})
## run from the repository root, as tools/exactproof.R is
root <- "."
if (!dir.exists(file.path(root, "R"))) stop("run this from the repository root", call. = FALSE)
rd  <- readLines(file.path(root, "README.md"), warn = FALSE)
rdt <- paste(rd, collapse = "\n")
## read every rendered report, not only the macOS one. A claim measured on Linux is verifiable
## only against the report that produced it, and reading one file failed four true claims.
htmlf <- list.files(file.path(root, "docs"), pattern = "\\.html$", full.names = TRUE)
if (!length(htmlf)) stop("docs/ holds no rendered report; render one first", call. = FALSE)
html <- paste(unlist(lapply(htmlf, readLines, warn = FALSE)), collapse = "\n")

pass <- 0L; fail <- 0L
chk <- function(ok, what, detail = "") {
  if (isTRUE(ok)) pass <<- pass + 1L
  else { fail <<- fail + 1L; cat(sprintf("  *** %s%s\n", what, if (nzchar(detail)) paste0(" -- ", detail) else "")) }
}

## 1. every speedup the README asserts must appear in the rendered report ----------------
if (nzchar(html)) {
  # The report prints full precision (4.408462) where the README rounds (4.41x), so compare
  # rounded values rather than substrings. A substring test silently passed everything that
  # happened to share a prefix and failed everything that did not, which is worse than nothing.
  nums <- as.numeric(regmatches(html, gregexpr("[0-9]+\\.[0-9]+", html))[[1]])
  nums <- nums[is.finite(nums)]
  seen <- unique(round(nums, 2L))
  # numbers measured elsewhere in the repo and documented there, not produced by the report
  elsewhere <- as.numeric(sub("x$", "", unique(regmatches(
    paste(src0 <- unlist(lapply(list.files(file.path(root, "R"), full.names = TRUE),
                                readLines, warn = FALSE)), collapse = "\n"),
    gregexpr("[0-9]+\\.[0-9]+x", paste(src0, collapse = "\n")))[[1]])))
  claims <- unique(regmatches(rdt, gregexpr("[0-9]+\\.[0-9]+x", rdt))[[1]])
  for (c in claims) {
    n <- round(as.numeric(sub("x$", "", c)), 2L)
    chk(n %in% seen || n %in% round(elsewhere, 2L),
        sprintf("README claims %s, which is in neither the rendered report nor R/", c))
  }
} else stop("the rendered reports are empty; render them first", call. = FALSE)

## 2. worker default -------------------------------------------------------------------
## every entry point must defer to the same resolver, or the README's one sentence about the
## default is true of some of them and quietly false of the rest
fns <- c("ComBat_seq_parallel", "lmFit_parallel", "duplicateCorrelation_parallel",
         "calcNormFactors_parallel")
defs <- vapply(fns, function(f) is.null(formals(get(f))$workers), logical(1))
chk(all(defs), "workers no longer defaults to NULL everywhere",
    paste(names(defs)[!defs], collapse = "/"))
want <- max(1L, min(8L, max(1L, parallel::detectCores()) - 2L))
## Without fork() the pick is capped again at the performance-core count, because an efficiency
## core that has to be handed a serialised copy stops paying for itself. The check follows the
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

## 4. options: every combat.* the code reads, and nothing the code does not -------------
src <- unlist(lapply(list.files(file.path(root, "R"), full.names = TRUE), readLines, warn = FALSE))
# require a real trailing segment: a bare "combat.min." prefix appears inside message text
opt_re <- 'combat\\.[a-z]+(\\.[a-z]+)+'
in_code <- sort(unique(regmatches(paste(src, collapse = "\n"),
                gregexpr(opt_re, paste(src, collapse = "\n")))[[1]]))
in_doc <- sort(unique(regmatches(rdt, gregexpr(opt_re, rdt))[[1]]))
for (o in in_doc) {
  # a retired name is allowed where the README is explicitly documenting the rename
  retired <- grepl(paste0(o, "` became|", o, "` is no longer|", o, "`, renamed"), rdt)
  chk(o %in% in_code || retired,
      sprintf("README names option %s which no longer exists in R/", o))
}
gates <- sum(grepl("^combat\\.min\\.", in_code))
if (grepl("Seven `combat[.]min", rdt))
  chk(gates == 7L, "README says seven size gates but the code disagrees",
      sprintf("code has %d", gates))

## 5. exported functions named in README must exist -------------------------------------
ex <- getNamespaceExports("rnaparallel")
for (f in unique(regmatches(rdt, gregexpr("`[a-zA-Z_]+_parallel\\(\\)`|`combat_[a-z_]+\\(\\)`", rdt))[[1]])) {
  nm <- gsub("[`()]", "", f)
  chk(nm %in% ex, sprintf("README documents %s which is not exported", nm))
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
## fail whenever source was ahead of the library -- which is the normal state while preparing
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
gated <- unique(unlist(regmatches(src, gregexpr('rebound = c\\([^)]*\\)', src))))
gated <- unique(unlist(regmatches(gated, gregexpr('"[A-Za-z_.]+"', gated))))
gated <- gsub('"', "", gated)
chk(length(gated) > 0L, "no rebound = c(...) gate list found in R/")
## The direction that matters is docs -> code: NEWS naming a symbol as gated when it is not is
## a false assurance about the one guard that catches a silently-serial rebind. The reverse
## (a gated symbol NEWS does not mention) is normal, since NEWS describes what changed.
claimed <- unlist(regmatches(current, gregexpr("gate covers[^.]*", current)))
claimed <- gsub("`", "", unlist(regmatches(claimed, gregexpr("`[A-Za-z_.]+`", claimed))))
for (g in unique(claimed)) {
  chk(g %in% gated,
      sprintf("NEWS says the reachability gate covers `%s`, which is in no rebound list in R/", g))
}

cat(sprintf("\n==== doccheck: %d checks, %d passed, %d failed ====\n", pass + fail, pass, fail))
quit(status = if (fail == 0L) 0L else 1L)
