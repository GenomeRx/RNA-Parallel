## calcNormFactors_parallel against edgeR on upstream drift, edgeR-accepted inputs and S3 methods.

skip_if_no_edger <- function() {
  testthat::skip_if_not_installed("edgeR")
}

sim_g3 <- function(G = 400L, S = 12L, seed = 1L) {
  set.seed(seed)
  cts <- matrix(stats::rnbinom(G * S, mu = 50, size = 5), G, S)
  dimnames(cts) <- list(sprintf("g%04d", seq_len(G)), sprintf("s%02d", seq_len(S)))
  cts
}

edger_default <- function() {
  get(paste0(rnaparallel:::rp_edger_generic(), ".default"), envir = asNamespace("edgeR"),
      inherits = FALSE)
}

## An edgeR internal with one piece of body text replaced, wired into the default method as `backend =`.
drifted_backend <- function(internal, from, to) {
  ns <- asNamespace("edgeR")
  f <- get(internal, envir = ns, inherits = FALSE)
  txt <- paste(deparse(body(f), width.cutoff = 500L), collapse = "\n")
  testthat::skip_if_not(grepl(from, txt, fixed = TRUE),
                        paste("edgeR's", internal, "no longer contains", from))
  body(f) <- parse(text = sub(from, to, txt, fixed = TRUE))[[1L]]
  e <- new.env(parent = ns)
  assign(internal, f, envir = e)
  fn <- edger_default()
  environment(fn) <- e
  fn
}

spy <- function() {
  n <- 0L
  list(backend = function(idx, f, workers) { n <<- n + 1L; lapply(idx, f) },
       count = function() n)
}

## The value or error message plus every warning and message, so two calls compare on all a user sees.
conditions <- function(expr) {
  w <- character()
  m <- character()
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "g3_error")),
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") },
    message = function(c) { m <<- c(m, conditionMessage(c)); invokeRestart("muffleMessage") })
  list(value = value, warnings = w, messages = m)
}

timing_line <- function(expr) {
  withr::local_options(list(combat.timing = TRUE, combat.timing.min = 0))
  paste(testthat::capture_messages(suppressWarnings(expr)), collapse = "")
}


# ---- defaults the caller did not pass are the backend's -----------------------

test_that("a scalar default the caller did not pass comes from the backend, not the companion", {
  skip_if_no_edger()
  cts <- sim_g3()
  base <- edger_default()
  drifts <- list(list("logratioTrim", 0.25, "TMM"), list("sumTrim", 0.1, "TMM"),
                 list("doWeighting", FALSE, "TMM"), list("Acutoff", -9, "TMM"),
                 list("p", 0.9, "upperquartile"))
  for (d in drifts) {
    b <- base
    formals(b)[[d[[1L]]]] <- d[[2L]]
    ref <- b(cts, method = d[[3L]])
    expect_false(identical(ref, base(cts, method = d[[3L]])), info = d[[1L]])
    expect_identical(calcNormFactors_parallel(cts, method = d[[3L]], workers = 2L, backend = b),
                     ref, info = d[[1L]])
  }
  b <- base
  formals(b)$logratioTrim <- 0.25
  expect_identical(calcNormFactors_parallel(cts, workers = 2L, backend = b), b(cts))
})

test_that("a reordered method default follows the backend instead of erroring", {
  skip_if_no_edger()
  cts <- sim_g3()
  b <- edger_default()
  formals(b)$method <- c("TMMwsp", "TMM", "RLE", "upperquartile", "none")
  expect_identical(calcNormFactors_parallel(cts, workers = 2L, backend = b), b(cts))
})

test_that("a DGEList takes its defaults from edgeR's own DGEList method", {
  skip_if_no_edger()
  ns <- asNamespace("edgeR")
  nm <- paste0(rnaparallel:::rp_edger_generic(), ".DGEList")
  skip_if(!exists(nm, envir = ns, inherits = FALSE), "no DGEList method in this edgeR")
  dge <- edgeR::DGEList(sim_g3())
  orig <- get(nm, envir = ns, inherits = FALSE)
  drifted <- orig
  formals(drifted)$logratioTrim <- 0.25
# UseMethod reads the S3 table, which still holds the original, so the reference calls the drift directly.
  ref <- drifted(dge)
  expect_false(identical(ref, orig(dge)))
  withr::defer({ assign(nm, orig, envir = ns); lockBinding(nm, ns) })
  unlockBinding(nm, ns)
  assign(nm, drifted, envir = ns)
  expect_identical(calcNormFactors_parallel(dge, workers = 2L), ref)
})

test_that("every normLibSizes.default argument exists on the companion with the same default", {
  skip_if_no_edger()
  vf <- formals(edger_default())
  pf <- formals(calcNormFactors_parallel)
  expect_true(all(names(vf) %in% names(pf)))
  for (a in names(vf)) {
    expect_identical(deparse(pf[[a]]), deparse(vf[[a]]), info = paste("default differs for", a))
  }
})


# ---- the TMMwzp alias --------------------------------------------------------

test_that("method = \"TMMwzp\" behaves exactly as edgeR does, rename message included", {
  skip_if_no_edger()
  cts <- sim_g3()
  for (obj in list(cts, edgeR::DGEList(cts))) {
    ref <- conditions(edgeR_norm(obj, method = "TMMwzp"))
    expect_false(inherits(ref$value, "g3_error"))
    expect_identical(conditions(calcNormFactors_parallel(obj, method = "TMMwzp", workers = 2L)),
                     ref, info = class(obj)[1L])
  }
})

test_that("an invalid method fails with edgeR's own error, not the companion's", {
  skip_if_no_edger()
  cts <- sim_g3()
  expect_identical(conditions(calcNormFactors_parallel(cts, method = "nope", workers = 2L)),
                   conditions(edgeR_norm(cts, method = "nope")))
})


# ---- the RLE apply shim ------------------------------------------------------

test_that("stock RLE, upperquartile and TMMwsp each reach the column split once", {
  skip_if_no_edger()
  cts <- sim_g3()
  for (m in c("RLE", "upperquartile", "TMMwsp")) {
    s <- spy()
    got <- calcNormFactors_parallel(cts, method = m, workers = 2L, chunks = 4L,
                                    parallel_backend = s$backend)
    expect_identical(got, edgeR_norm(cts, method = m), info = m)
    expect_identical(s$count(), 1L, info = m)
  }
})

test_that("the RLE column split survives an integer MARGIN and extra apply arguments", {
  skip_if_no_edger()
  cts <- sim_g3()
  drifts <- list(
    c("apply(data, 2,", "apply(data, 2L,"),
    c("apply(data, 2, function(u) median((u/gm)[gm > 0]))",
      "apply(data, 2, function(u, g) median((u/g)[g > 0]), g = gm)"))
  for (d in drifts) {
    b <- drifted_backend(".calcFactorRLE", d[1L], d[2L])
    s <- spy()
    got <- calcNormFactors_parallel(cts, method = "RLE", workers = 2L, chunks = 4L,
                                    parallel_backend = s$backend, backend = b)
    expect_identical(got, b(cts, method = "RLE"), info = d[2L])
    expect_identical(s$count(), 1L, info = d[2L])
  }
})

test_that("an RLE apply the shim cannot split is reported as stood down", {
  skip_if_no_edger()
  cts <- sim_g3()
  stock <- timing_line(calcNormFactors_parallel(cts, method = "RLE", workers = 2L))
  expect_false(grepl("stood down", stock, fixed = TRUE))
  b <- drifted_backend(".calcFactorRLE", "apply(data, 2,", "apply(t(data), 1,")
  got <- NULL
  line <- timing_line(got <- calcNormFactors_parallel(cts, method = "RLE", workers = 2L,
                                                      backend = b))
  expect_identical(got, b(cts, method = "RLE"))
  expect_match(line, "RLE columns stood down", fixed = TRUE)
})


# ---- the quantile shim -------------------------------------------------------

test_that("the quantile column split survives extra quantile arguments", {
  skip_if_no_edger()
  cts <- sim_g3()
  for (to in c("probs = p, names = FALSE)", "probs = p, type = 7)")) {
    b <- drifted_backend(".calcFactorQuantile", "probs = p)", to)
    for (m in c("upperquartile", "TMM")) {
      s <- spy()
      got <- calcNormFactors_parallel(cts, method = m, workers = 2L, chunks = 4L,
                                      parallel_backend = s$backend, backend = b)
      expect_identical(got, b(cts, method = m), info = paste(to, m))
      expect_identical(s$count(), if (m == "TMM") 2L else 1L, info = paste(to, m))
    }
  }
})

test_that("a quantile call the shim cannot batch runs as edgeR does and is reported", {
  skip_if_no_edger()
  cts <- sim_g3()
  two_p <- function() calcNormFactors_parallel(cts, method = "upperquartile", p = c(0.5, 0.75),
                                               workers = 2L)
  expect_identical(conditions(two_p()),
                   conditions(edgeR_norm(cts, method = "upperquartile", p = c(0.5, 0.75))))
  expect_match(timing_line(two_p()), "quantile columns stood down", fixed = TRUE)

  b <- drifted_backend(".calcFactorQuantile", "quantile(data[, j], probs = p)",
                       "quantile(data[, j])")
  expect_identical(
    conditions(calcNormFactors_parallel(cts, method = "upperquartile", workers = 2L, backend = b)),
    conditions(b(cts, method = "upperquartile")))
})


# ---- refColumn of length other than one --------------------------------------

test_that("a refColumn of length other than one is left to edgeR, warnings and errors included", {
  skip_if_no_edger()
  cts <- sim_g3(S = 8L)
  cases <- list(list(refColumn = c(1, 2)), list(refColumn = integer(0)),
                list(refColumn = rep_len(1:2, nrow(cts))),
                list(method = "TMMwsp", refColumn = c(1, 2)))
  for (a in cases) {
    ref <- conditions(do.call(edgeR_norm, c(list(cts), a)))
    got <- conditions(do.call(calcNormFactors_parallel, c(list(cts, workers = 2L), a)))
    expect_identical(got, ref, info = substr(deparse1(a), 1L, 80L))
  }
  dge <- edgeR::DGEList(cts)
  expect_identical(conditions(calcNormFactors_parallel(dge, refColumn = c(1, 2), workers = 2L)),
                   conditions(edgeR_norm(dge, refColumn = c(1, 2))))
})


# ---- S3 methods registered by another package --------------------------------

test_that("an S3 method registered on edgeR's generic by another package is refused", {
  skip_if_no_edger()
  ns <- asNamespace("edgeR")
  generic <- rnaparallel:::rp_edger_generic()
  cl <- "rnaparallelG3Counts"
  tbl <- get(".__S3MethodsTable__.", envir = ns, inherits = FALSE)
# Written straight into the table, because registerS3method() also adds a row to edgeR's S3methods info that the defer cannot undo and that breaks every later assignInNamespace() on edgeR.
  assign(paste0(generic, ".", cl), function(object, ...) stop("the registered method ran"),
         envir = tbl)
  withr::defer(rm(list = paste0(generic, ".", cl), envir = tbl))
  cts <- sim_g3()
  obj <- structure(cts, class = c(cl, "matrix", "array"))
  expect_error(calcNormFactors_parallel(obj, workers = 2L), "does not wrap")
  expect_identical(calcNormFactors_parallel(cts, workers = 2L), edgeR_norm(cts))
})
