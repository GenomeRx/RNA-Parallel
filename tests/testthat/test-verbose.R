# This file defines its own skip_if_no_limma() and sim(), so it does not depend on another file's load order.
skip_if_no_limma <- function() {
  testthat::skip_if_not_installed("limma")
  testthat::skip_if_not_installed("edgeR")
}

sim <- function(G = 200L, S = 12L, seed = 1L) {
  set.seed(seed)
  lib <- round(exp(stats::rnorm(S, log(2e6), 0.3)))
  p   <- 2^pmax(stats::rnorm(G, 3, 2), -1); p <- p / sum(p)
  cts <- matrix(stats::rnbinom(G * S, mu = outer(p, lib), size = 10), G, S)
  rownames(cts) <- sprintf("g%04d", seq_len(G))
  colnames(cts) <- sprintf("s%02d", seq_len(S))
  list(counts = cts)
}

# Resolved here, not taken from helper-counts.R, so the original this file compares against does not depend on load order.
vendor_norm <- function(object, ...) {
  ns <- asNamespace("edgeR")
  nm <- if (exists("normLibSizes.default", envir = ns, inherits = FALSE)) "normLibSizes"
        else "calcNormFactors"
  get(nm, envir = ns, inherits = FALSE)(object, ...)
}

hush <- function(expr) {
  out <- NULL
  invisible(utils::capture.output(out <- suppressMessages(expr)))
  out
}

test_that("timing and quieting are off by default and change nothing", {
  skip_if_no_limma()
  d <- sim()
  expect_null(getOption("combat.timing"))
  expect_null(getOption("combat.quiet"))
  ref <- vendor_norm(d$counts)
  expect_silent(got <- calcNormFactors_parallel(d$counts, workers = 2L))
  expect_identical(got, ref)
})

test_that("a timed run returns exactly what an untimed one does", {
  skip_if_no_limma()
  d <- sim()
  quiet_run <- calcNormFactors_parallel(d$counts, workers = 2L)
  withr::local_options(list(combat.timing = TRUE))
  expect_message(loud_run <- calcNormFactors_parallel(d$counts, workers = 2L))
  # the whole product is identical() output; a progress line must not touch the value
  expect_identical(loud_run, quiet_run)
})

test_that("one user-facing call prints one line, not one per nested call", {
  skip_if_no_limma()
  d <- sim()
  withr::local_options(list(combat.timing = TRUE))
# calcNormFactors_parallel on a DGEList reaches edgeR's DGEList method, which calls the companion again on the counts matrix.
  msgs <- testthat::capture_messages(
    calcNormFactors_parallel(edgeR::DGEList(d$counts), workers = 2L))
  expect_length(grep("TMM", msgs), 1L)
})

test_that("the timing line names what actually ran, not what was asked for", {
  skip_if_no_limma()
  d <- sim()
  withr::local_options(list(combat.timing = TRUE, combat.fork = FALSE))
# Every stage is gated to serial here, so the line must say serial and never "mclapply x2".
  msgs <- testthat::capture_messages(calcNormFactors_parallel(d$counts, workers = 2L))
  expect_match(paste(msgs, collapse = ""), "serial")
  expect_match(paste(msgs, collapse = ""), "gated")
})

test_that("combat.timing.min hides steps faster than the threshold", {
  skip_if_no_limma()
  d <- sim()
  withr::local_options(list(combat.timing = TRUE, combat.timing.min = 3600))
  expect_silent(calcNormFactors_parallel(d$counts, workers = 2L))
})

test_that("quieting swallows the original's cat() but never a warning or an error", {
  skip_if_no_limma()
  skip_if_not_installed("sva")
  set.seed(4)
  cts <- matrix(rnbinom(1600, mu = 50, size = 5), nrow = 200)
  bt  <- rep(1:2, each = 4)
  withr::local_options(list(combat.quiet = TRUE))
# sva announces itself with cat(), which only a stdout sink silences, and stderr stays open so a warning about the data still lands.
  out <- utils::capture.output(
    got <- ComBat_seq_parallel(cts, batch = bt, group = NULL, workers = 2L))
  expect_length(out, 0L)
  expect_identical(got, hush(sva::ComBat_seq(cts, batch = bt, group = NULL)))
})

test_that("a failed call unwinds the sink instead of silencing the session", {
  skip_if_no_limma()
  withr::local_options(list(combat.timing = TRUE, combat.quiet = TRUE))
  before <- sink.number()
# A dangling sink would swallow every later print in the caller's session.
  expect_error(suppressMessages(calcNormFactors_parallel("not a matrix", workers = 2L)))
  expect_identical(sink.number(), before)
  expect_identical(rnaparallel:::rp_or0(rnaparallel:::.rp_dispatch$depth), 0L)
})

test_that("a garbage timing option is refused rather than ignored", {
  withr::local_options(list(combat.timing = "yes"))
  expect_error(rnaparallel:::rp_opt_flag("combat.timing"), "TRUE or FALSE")
  withr::local_options(list(combat.timing = TRUE, combat.timing.min = "soon"))
  expect_error(rnaparallel:::rp_opt_secs("combat.timing.min"), "non-negative")
})

test_that("label overrides the derived name", {
  skip_if_no_limma()
  d <- sim()
  withr::local_options(list(combat.timing = TRUE))
  msgs <- testthat::capture_messages(
    calcNormFactors_parallel(d$counts, workers = 2L, label = "Breast cohort"))
  expect_match(paste(msgs, collapse = ""), "Breast cohort", fixed = TRUE)
})

test_that("rnaparallel_stale reports on the build, not on nothing", {
# A package loaded by load_all() has no Built field, so NA is the answer there, and against a real install it must be FALSE.
  v <- rnaparallel_stale()
  expect_true(is.logical(v) && length(v) == 1L)
  expect_false(isTRUE(v))
})

test_that("a backend that falls back to serial is not reported as parallel", {
  skip_if_no_limma()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  d <- sim()
  withr::local_options(list(combat.timing = TRUE, combat.min.norm.cols = 0,
                            combat.min.norm.cells = 0))
# A future plan that resolves in one process must be counted serial, never as "future x2" beside the warning that it ran serially.
  future::plan(future::sequential)
  withr::defer(future::plan(future::sequential))
  msgs <- testthat::capture_messages(suppressWarnings(
    calcNormFactors_parallel(d$counts, workers = 2L, parallel_backend = "future")))
  line <- paste(msgs, collapse = "")
  expect_match(line, "serial")
  expect_false(grepl("future x", line, fixed = TRUE))
})

test_that("every backend produces a timing line and the same answer", {
  skip_if_no_limma()
  d <- sim()
  ref <- vendor_norm(d$counts)
  withr::local_options(list(combat.timing = TRUE, combat.min.norm.cols = 0,
                            combat.min.norm.cells = 0))
  for (b in c("serial", "mclapply")) {
    msgs <- testthat::capture_messages(
      got <- calcNormFactors_parallel(d$counts, workers = 2L, parallel_backend = b))
    expect_identical(got, ref, info = b)
    expect_match(paste(msgs, collapse = ""), "TMM", info = b)
  }
  # a function backend is labeled "custom", not by a name it does not have
  msgs <- testthat::capture_messages(
    got <- calcNormFactors_parallel(d$counts, workers = 2L,
                                    parallel_backend = function(idx, f, w) lapply(idx, f)))
  expect_identical(got, ref)
  expect_match(paste(msgs, collapse = ""), "custom")
})

test_that("a pinned original excerpt standing down is reported, not silent", {
  skip_if_not_installed("sva")
  withr::local_options(list(combat.timing = TRUE, combat.quiet = TRUE,
                            combat.min.cells = 0, combat.min.glm.cells = 0,
                            combat.min.disp.cells = 0))
  set.seed(3)
  cts <- matrix(rnbinom(3200, mu = 50, size = 5), 400, 8)
  bt <- rep(1:2, each = 4)

  clean <- testthat::capture_messages(
    a <- ComBat_seq_parallel(cts, batch = bt, group = NULL, workers = 2L))
  expect_false(grepl("stood down", paste(clean, collapse = "")))

# An sva whose match_quantiles body only reformats shuts the byte-exact gate, and that slower run with identical output must be reported.
  ns <- asNamespace("sva")
  orig <- get("match_quantiles", envir = ns)
  drift <- orig
  body(drift) <- parse(text = paste0("{ ", paste(deparse(body(orig)), collapse = "\n"), " }"))[[1L]]
  withr::defer({ assign("match_quantiles", orig, envir = ns); lockBinding("match_quantiles", ns) })
  unlockBinding("match_quantiles", ns); assign("match_quantiles", drift, envir = ns)

  loud <- testthat::capture_messages(
    b <- ComBat_seq_parallel(cts, batch = bt, group = NULL, workers = 2L))
  expect_match(paste(loud, collapse = ""), "match_quantiles stood down", fixed = TRUE)
  expect_identical(a, b)          # the whole point: slower, never different
})

test_that("rnaparallel_stale catches a namespace/disk version mismatch", {
  skip_on_cran()
# The comparison needs nothing recorded at load, because an old namespace still has a version and packageVersion() reads the new one off disk.
  expect_false(isTRUE(rnaparallel_stale()))

  ns <- asNamespace("rnaparallel")
  real <- get("getNamespaceVersion", envir = baseenv())
  local_ns <- new.env(parent = ns)
  local_ns$getNamespaceVersion <- function(...) "0.0.1"     # pretend the session is ancient
  f <- rnaparallel_stale
  environment(f) <- local_ns
  expect_true(f())
  expect_identical(real(ns), getNamespaceVersion("rnaparallel"))   # nothing was mutated
})

test_that("the help page shows the defensive call, and it reads FALSE on a current session", {
  documented <- quote(isTRUE(tryCatch(rnaparallel_stale(), error = function(e) TRUE)))
  root <- rp_dev_root()
  db <- if (nzchar(root)) tools::Rd_db(dir = root) else tools::Rd_db("rnaparallel")
  rd <- paste(as.character(db[["rnaparallel_stale.Rd"]]), collapse = "")
  expect_true(grepl(paste(deparse(documented), collapse = ""), rd, fixed = TRUE))
  expect_false(eval(documented, new.env(parent = asNamespace("rnaparallel"))))
})
