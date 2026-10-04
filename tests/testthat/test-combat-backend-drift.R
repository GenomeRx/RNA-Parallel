# Each drift test passes a body-edited sva::ComBat_seq through `backend =` and expects the companion to match that copy.

drift_body <- function(f, from, to) {
  txt <- paste(deparse(body(f), width.cutoff = 500L), collapse = "\n")
  stopifnot(grepl(from, txt, fixed = TRUE))
  g <- f
  body(g) <- parse(text = sub(from, to, txt, fixed = TRUE))[[1L]]
  environment(g) <- environment(f)
  g
}

drift_counts_name <- function(f) {
  txt <- paste(deparse(body(f), width.cutoff = 500L), collapse = "\n")
  txt <- gsub("(?<![$\\w.])counts(?![\\w.])(?!\\s*=[^=])", "cts0", txt, perl = TRUE)
  g <- f
  body(g) <- as.call(c(as.name("{"), quote(cts0 <- counts), as.list(parse(text = txt)[[1L]])[-1L]))
  environment(g) <- environment(f)
  g
}

in_process <- function(idx, f, workers) lapply(idx, f)

test_that("extra arguments on the tagwise call reach edgeR instead of being dropped", {
  skip_if_not_installed("sva")
  d <- drift_body(sva::ComBat_seq, "prior.df = 0))", "prior.df = 0, min.row.sum = 1000))")
  set.seed(4)
  cts <- matrix(rnbinom(300 * 30, mu = 60, size = 4), nrow = 300)
  bat <- rep(1:3, each = 10L)
  ref <- quietly(d(cts, batch = bat, group = NULL))
  expect_false(identical(ref, quietly(sva::ComBat_seq(cts, batch = bat, group = NULL))))

  rnaparallel:::rp_count_reset()
  got <- quietly(ComBat_seq_parallel(cts, batch = bat, group = NULL, workers = 2L,
                                     parallel_backend = in_process, backend = d))
  expect_identical(got, ref)
  expect_true("tagwise rows" %in% rnaparallel:::.rp_dispatch$fallback)

  expect_identical(quietly(ComBat_seq_parallel(cts, batch = bat, group = NULL, workers = 2L,
                                               backend = d)), ref)
})

test_that("extra tagwise arguments are honored inside a socket worker too", {
  skip_on_cran()
  skip_if_not_installed("sva")
  skip_if_not_installed("future.apply")
  d <- drift_body(sva::ComBat_seq, "prior.df = 0))", "prior.df = 0, min.row.sum = 1000))")
  set.seed(4)
  cts <- matrix(rnbinom(300 * 30, mu = 60, size = 4), nrow = 300)
  bat <- rep(1:3, each = 10L)
  ref <- quietly(d(cts, batch = bat, group = NULL))
  local_socket_plan(2L)
  expect_identical(quietly(ComBat_seq_parallel(cts, batch = bat, group = NULL, workers = 2L,
                                               parallel_backend = "future", backend = d)), ref)
})

test_that("an argument the caller omits takes the backend's own default", {
  skip_if_not_installed("sva")
  b2 <- sva::ComBat_seq
  formals(b2)$full_mod <- FALSE
  set.seed(4)
  cts <- matrix(rnbinom(300 * 30, mu = 60, size = 4), nrow = 300)
  bat <- rep(1:3, each = 10L)
  grp <- rep(1:2, 15)
  ref <- quietly(b2(cts, batch = bat, group = grp))
  expect_false(identical(ref, quietly(sva::ComBat_seq(cts, batch = bat, group = grp))))

  expect_identical(quietly(ComBat_seq_parallel(cts, batch = bat, group = grp, workers = 2L,
                                               backend = b2)), ref)
  expect_identical(quietly(ComBat_seq_parallel(cts, batch = bat, group = grp, full_mod = TRUE,
                                               workers = 2L, backend = b2)),
                   quietly(b2(cts, batch = bat, group = grp, full_mod = TRUE)))
})

test_that("a tagwise shape check that cannot find the filtered matrix says it stood down", {
  skip_if_not_installed("sva")
  d <- drift_counts_name(sva::ComBat_seq)
  set.seed(7)
  y <- matrix(rnbinom(300 * 12, mu = 40, size = 4), 300, 12)
  batch <- rep(1:2, each = 6)
  y[5, batch == 1] <- 0L
  ref <- quietly(sva::ComBat_seq(y, batch = batch, group = NULL))
  expect_identical(quietly(d(y, batch = batch, group = NULL)), ref)

  rnaparallel:::rp_count_reset()
  stock <- quietly(ComBat_seq_parallel(y, batch = batch, group = NULL, workers = 2L,
                                       parallel_backend = in_process))
  expect_identical(stock, ref)
  expect_false("tagwise across batches" %in% rnaparallel:::.rp_dispatch$fallback)

  rnaparallel:::rp_count_reset()
  drifted <- quietly(ComBat_seq_parallel(y, batch = batch, group = NULL, workers = 2L,
                                         parallel_backend = in_process, backend = d))
  expect_identical(drifted, ref)
  expect_true("tagwise across batches" %in% rnaparallel:::.rp_dispatch$fallback)
})

test_that("the per-batch jobs ship only what the original's closure reads", {
  skip_if_not_installed("sva")
  set.seed(3)
  cts <- matrix(rnbinom(2000 * 30, mu = 60, size = 4), 2000, 30)
  storage.mode(cts) <- "integer"
  bat <- rep(1:3, each = 10L)
  box <- new.env(parent = globalenv())
  box$mat_bytes <- length(serialize(cts, NULL))
  box$ratios <- numeric()
# load_all keeps srcrefs whose parse data is several matrices' worth here and absent once installed, so srcfile environments are written as a name.
  box$nosrc <- function(e) if (inherits(e, "srcfile")) "srcfile" else NULL
  rec <- function(idx, f, workers) {
    if (all(lengths(idx) == 1L)) {
      ratios <<- c(ratios, length(serialize(f, NULL, refhook = nosrc)) / mat_bytes)
    }
    lapply(idx, f)
  }
  environment(rec) <- box
  environment(box$nosrc) <- box
  out <- quietly(ComBat_seq_parallel(cts, batch = bat, group = NULL, workers = 2L,
                                     parallel_backend = rec))
  expect_identical(out, quietly(sva::ComBat_seq(cts, batch = bat, group = NULL)))
  expect_length(box$ratios, 2L)
  expect_true(all(box$ratios < 2),
              info = paste("payload / matrix:", paste(round(box$ratios, 2), collapse = " ")))
})

test_that("a per-batch closure that looks its matrix up by string ships unchanged", {
  skip_if_not_installed("sva")
  txt <- paste(deparse(body(sva::ComBat_seq), width.cutoff = 500L), collapse = "\n")
  txt <- gsub("estimateGLMCommonDisp(counts[", "estimateGLMCommonDisp(get(\"counts\")[", txt, fixed = TRUE)
  txt <- gsub("subset = nrow(counts)", "subset = nrow(get(\"counts\"))", txt, fixed = TRUE)
  d <- sva::ComBat_seq
  body(d) <- parse(text = txt)[[1L]]
  environment(d) <- environment(sva::ComBat_seq)
  set.seed(4)
  cts <- matrix(rnbinom(300 * 30, mu = 60, size = 4), nrow = 300)
  bat <- rep(1:3, each = 10L)
  ref <- quietly(d(cts, batch = bat, group = NULL))
  expect_identical(ref, quietly(sva::ComBat_seq(cts, batch = bat, group = NULL)))
  expect_identical(quietly(ComBat_seq_parallel(cts, batch = bat, group = NULL, workers = 2L,
                                               parallel_backend = in_process, backend = d)), ref)
  expect_identical(quietly(ComBat_seq_parallel(cts, batch = bat, group = NULL, workers = 2L,
                                               backend = d)), ref)
})
