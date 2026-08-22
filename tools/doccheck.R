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
htmlf <- file.path(root, "inst/examples/RNA_Parallel.html")
html <- if (file.exists(htmlf)) paste(readLines(htmlf, warn = FALSE), collapse = "\n") else ""

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
} else cat("  (no rendered report found, skipping number cross-check)\n")

## 2. worker default -------------------------------------------------------------------
defs <- vapply(c("ComBat_seq_parallel", "lmFit_parallel", "duplicateCorrelation_parallel",
                 "calcNormFactors_parallel"),
               function(f) as.integer(eval(formals(get(f))$workers)), integer(1))
chk(all(defs == 4L), "workers default is not four everywhere", paste(defs, collapse = "/"))
chk(grepl("defaults to four", rdt), "README no longer states the workers default")

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
ver <- as.character(packageVersion("rnaparallel"))
news <- readLines(file.path(root, "NEWS.md"), warn = FALSE)[1]
chk(grepl(ver, news, fixed = TRUE),
    "NEWS.md does not open with the installed version", sprintf("installed %s, NEWS '%s'", ver, news))

cat(sprintf("\n==== doccheck: %d checks, %d passed, %d failed ====\n", pass + fail, pass, fail))
quit(status = if (fail == 0L) 0L else 1L)
