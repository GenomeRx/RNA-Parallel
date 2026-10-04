test_that("a new_mu narrower than the counts errors as the original does, not with NA cells", {
  skip_if_not_installed("sva")
  mq <- get("match_quantiles", envir = asNamespace("sva"))
  set.seed(3)
  cs <- matrix(rnbinom(300, mu = 20, size = 4), 50, 6)
  om <- matrix(20, 50, 6)
  op <- rep(0.25, 50)
  nm <- matrix(30, 50, 5)
  np <- rep(0.2, 50)
  expect_error(mq(cs, om, op, nm, np), "subscript out of bounds")
  expect_error(rnaparallel:::match_quantiles_parallel(mq, cs, om, op, nm, np, workers = 2L,
                                                      parallel_backend = "serial"),
               "subscript out of bounds")
  expect_null(rnaparallel:::combat_mq_dispatch(mq, cs, om, op, nm, np))
  expect_null(rnaparallel:::combat_mq_dispatch(mq, cs, om, op, matrix(30, 50, 6), np[-1]))
  expect_null(rnaparallel:::combat_mq_dispatch(mq, cs, om, op, as.data.frame(matrix(30, 50, 6)), np))
  expect_identical(rnaparallel:::combat_mq_dispatch(mq, cs, om, op, matrix(30, 50, 6), np),
                   rnaparallel:::match_quantiles_rows)
})

test_that("only a started BiocParallel param counts as copying for the performance-core note", {
  skip_on_os("windows")
  skip_if_not_installed("BiocParallel")
  cc <- rnaparallel:::.combat_clusters
  old_p <- cc$perfcores
  old_w <- cc$warned_ecores
  withr::defer({
    cc$perfcores <- old_p
    cc$warned_ecores <- old_w
  })
  cc$perfcores <- 1L
  run <- function(b) {
    rnaparallel:::combat_parallel_lapply(as.list(1:4), function(i) i * 2, 2L, b,
                                         cells = Inf, min_cells = 0)
  }
  cc$warned_ecores <- NULL
  expect_no_message(got <- run("mclapply"), message = "performance core")
  expect_identical(got, list(2, 4, 6, 8))
  cc$warned_ecores <- NULL
  expect_no_message(got <- run("BiocParallel"), message = "performance core")
  expect_identical(got, list(2, 4, 6, 8))
  p <- parallel::mcparallel(Sys.sleep(30))
  withr::defer({
    tools::pskill(p$pid, tools::SIGKILL)
    suppressWarnings(parallel::mccollect(p))
  })
  cc$warned_ecores <- NULL
  expect_message(got <- run("BiocParallel"), "performance core")
  expect_identical(got, list(2, 4, 6, 8))
})


# ---- conditions raised in workers reach the caller -------------------------------

r2_backends <- function() {
  b <- "serial"
  if (.Platform$OS.type != "windows") b <- c(b, "mclapply")
  if (requireNamespace("doParallel", quietly = TRUE)) b <- c(b, "foreach")
  if (requireNamespace("future.apply", quietly = TRUE)) b <- c(b, "future")
  if (.Platform$OS.type != "windows" && requireNamespace("BiocParallel", quietly = TRUE)) {
    b <- c(b, "BiocParallel")
  }
  b
}

r2_plan <- function(env = parent.frame()) {
  if (!requireNamespace("future.apply", quietly = TRUE)) return(invisible(NULL))
  local_socket_plan(2L, .env = env)
}

r2_conds <- function(f) {
  w <- character()
  m <- character()
  out <- NULL
  invisible(utils::capture.output(out <- withCallingHandlers(f(),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      m <<- c(m, conditionMessage(x))
      invokeRestart("muffleMessage")
    })))
  list(value = out, warnings = w, messages = m)
}

r2_drift <- function(f, edits) {
  txt <- paste(deparse(body(f), width.cutoff = 500L), collapse = "\n")
  for (from in names(edits)) {
    stopifnot(grepl(from, txt, fixed = TRUE))
    txt <- sub(from, edits[[from]], txt, fixed = TRUE)
  }
  g <- f
  body(g) <- parse(text = txt)[[1L]]
  environment(g) <- environment(f)
  g
}

r2_rename_counts <- function(f) {
  txt <- paste(deparse(body(f), width.cutoff = 500L), collapse = "\n")
  txt <- gsub("(?<![$\\w.])counts(?![\\w.])(?!\\s*=[^=])", "cts0", txt, perl = TRUE)
  g <- f
  body(g) <- as.call(c(as.name("{"), quote(cts0 <- counts), as.list(parse(text = txt)[[1L]])[-1L]))
  environment(g) <- environment(f)
  g
}

r2_noisy <- function() {
  r2_drift(sva::ComBat_seq, c(
    "sapply(1:n_batch, function(i) {" =
      "sapply(1:n_batch, function(i) {\nwarning(\"common batch \", i)",
    "lapply(1:n_batch, function(j) {" =
      "lapply(1:n_batch, function(j) {\nwarning(\"tagwise batch \", j)\nmessage(\"tagwise note \", j)"))
}

test_that("per-column edgeR warnings reach the caller once each on every backend", {
  set.seed(1)
  x <- matrix(rnbinom(400 * 12, mu = 50, size = 5), 400, 12)
  ref <- r2_conds(function() edgeR_norm(x, refColumn = 0))
  expect_length(ref$warnings, 12L)
  r2_plan()
  for (b in r2_backends()) {
    got <- r2_conds(function() calcNormFactors_parallel(x, refColumn = 0, workers = 2L,
                                                       parallel_backend = b))
    expect_identical(got$value, ref$value, info = b)
    expect_identical(got$warnings, ref$warnings, info = b)
    expect_identical(got$messages, ref$messages, info = b)
  }
})

test_that("warnings and messages raised inside ComBat-seq's batch jobs arrive in batch order", {
  skip_if_not_installed("sva")
  d <- r2_noisy()
  set.seed(4)
  cts <- matrix(rnbinom(300 * 30, mu = 60, size = 4), nrow = 300)
  bat <- rep(1:3, each = 10L)
  ref <- r2_conds(function() d(cts, batch = bat, group = NULL))
  expect_identical(ref$warnings, c(paste("common batch", 1:3), paste("tagwise batch", 1:3)))
  expect_identical(ref$messages, paste0("tagwise note ", 1:3, "\n"))
  r2_plan()
  for (b in r2_backends()) {
    got <- r2_conds(function() ComBat_seq_parallel(cts, batch = bat, group = NULL, workers = 2L,
                                                  parallel_backend = b, backend = d))
    expect_identical(got$value, ref$value, info = b)
    expect_identical(got$warnings, ref$warnings, info = b)
    expect_identical(got$messages, ref$messages, info = b)
  }
})

test_that("a batch dispatch that stands down and reruns serially warns once per batch", {
  skip_if_not_installed("sva")
  d <- r2_rename_counts(r2_noisy())
  set.seed(4)
  cts <- matrix(rnbinom(300 * 30, mu = 60, size = 4), nrow = 300)
  bat <- rep(1:3, each = 10L)
# Genes all zero in one batch are filtered into cts0, so the shape check reads the unfiltered counts, finds too many rows and stands down.
  cts[1:5, 11:20] <- 0L
  ref <- r2_conds(function() d(cts, batch = bat, group = NULL))
  expect_identical(ref$warnings, c(paste("common batch", 1:3), paste("tagwise batch", 1:3)))
  r2_plan()
  for (b in r2_backends()) {
    rnaparallel:::rp_count_reset()
    got <- r2_conds(function() ComBat_seq_parallel(cts, batch = bat, group = NULL, workers = 2L,
                                                  parallel_backend = b, backend = d))
    expect_identical(got$value, ref$value, info = b)
    expect_identical(got$warnings, ref$warnings, info = b)
    expect_identical(got$messages, ref$messages, info = b)
    expect_true("tagwise across batches" %in% rnaparallel:::.rp_dispatch$fallback, info = b)
  }
  got <- r2_conds(function() ComBat_seq_parallel(cts, batch = bat, group = NULL, workers = 1L,
                                                parallel_backend = "mclapply", backend = d))
  expect_identical(got$value, ref$value, info = "workers = 1")
  expect_identical(got$warnings, ref$warnings, info = "workers = 1")
})

test_that("a duplicateCorrelation warning raised in the blocks matches the original on every backend", {
  skip_if_not_installed("limma")
  skip_if_not_installed("statmod")
  set.seed(2)
  M <- matrix(stats::rnorm(60), 10, 6)
  ref <- r2_conds(function() limma::duplicateCorrelation(M, block = 1:6))
  expect_length(ref$warnings, 1L)
  r2_plan()
  for (b in r2_backends()) {
    got <- r2_conds(function() duplicateCorrelation_parallel(M, block = 1:6, workers = 2L,
                                                            parallel_backend = b))
    expect_identical(got$value, ref$value, info = b)
    expect_identical(got$warnings, ref$warnings, info = b)
  }
})

test_that("a fallback note raised in a worker reaches the master's timing line", {
  skip_if_not_installed("sva")
  d <- r2_drift(sva::ComBat_seq, c("prior.df = 0))" = "prior.df = 0, min.row.sum = 1000))"))
  set.seed(4)
  cts <- matrix(rnbinom(300 * 30, mu = 60, size = 4), nrow = 300)
  bat <- rep(1:3, each = 10L)
  ref <- quietly(d(cts, batch = bat, group = NULL))
  r2_plan()
  for (b in setdiff(r2_backends(), "serial")) {
    for (round in 1:2) {
      rnaparallel:::rp_count_reset()
      got <- quietly(ComBat_seq_parallel(cts, batch = bat, group = NULL, workers = 2L,
                                         parallel_backend = b, backend = d))
      expect_identical(got, ref, info = b)
      expect_true("tagwise rows" %in% rnaparallel:::.rp_dispatch$fallback,
                  info = paste(b, "round", round))
    }
  }
})
