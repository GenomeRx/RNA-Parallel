# Every equivalence test in this suite asserts identical(par, ref). A wrapper that did nothing
# but forward to sva::ComBat_seq would satisfy all of them, so none of them can tell a working
# parallel layer from a dead one. These tests watch the dispatch itself.

test_that("the entry point actually reaches the parallel layer", {
  skip_if_not_installed("sva")
  set.seed(4)
  cts <- matrix(rnbinom(300 * 30, mu = 60, size = 4), nrow = 300)
  bat <- rep(1:3, each = 10L)

  calls <- 0L
  spy <- function(idx, f, workers) { calls <<- calls + 1L; lapply(idx, f) }

  out <- ComBat_seq_parallel(cts, batch = bat, group = NULL, workers = 2L,
                             parallel_backend = spy)

  # These counts fell when the one-group gate landed. ComBat-seq builds a batch-only
  # design, which is exactly the layout edgeR sends to its mglmOneGroup kernel, and that
  # kernel is not a pure function of the gene it fits when the fit does not converge. So
  # glmFit and the tagwise dispersion now run whole and only the quantile match and the
  # common dispersion still dispatch: one per batch, plus one across batches.
  # Exact, not a bound: one match_quantiles per batch, one common dispersion across batches,
  # and one tagwise dispersion across batches. A lower bound is also satisfied with a rebind
  # deleted, so it would not prove the rebinds are reached.
  expect_identical(calls, 2L + length(unique(bat)))
  expect_identical(out, sva::ComBat_seq(cts, batch = bat, group = NULL))
})

test_that("the common dispersion is dispatched across batches, not once per run", {
  skip_if_not_installed("sva")
  set.seed(5)
  count_for <- function(nb) {
    cts <- matrix(rnbinom(300 * (nb * 10L), mu = 60, size = 4), nrow = 300)
    bat <- rep(seq_len(nb), each = 10L)
    calls <- 0L
    spy <- function(idx, f, workers) { calls <<- calls + 1L; lapply(idx, f) }
    ComBat_seq_parallel(cts, batch = bat, group = NULL, workers = 2L,
                        parallel_backend = spy)
    calls
  }
  # exact counts. "more batches means more dispatches" is true without the shim too, since
  # the quantile match already scales with batch count.
  expect_identical(count_for(2L), 4L)
  expect_identical(count_for(4L), 6L)
})

test_that("a backend that no longer calls the rebound symbols is refused", {
  # the failure this guards is silent: an unreachable rebind still returns identical output
  fake <- function(counts, batch, group = NULL, covar_mod = NULL, full_mod = TRUE,
                   shrink = FALSE, shrink.disp = FALSE, gene.subset.n = NULL) counts
  environment(fake) <- new.env(parent = globalenv())
  assign("match_quantiles",
         function(counts_sub, old_mu, old_phi, new_mu, new_phi) counts_sub,
         envir = environment(fake))

  expect_error(rnaparallel:::combat_backend(fake), "no longer calls")
})

test_that("cluster bookkeeping survives a non-cluster entry in the cache", {
  # The named flags are excluded from the sweep by name, so planting one of those would
  # never reach the is.list() guard this test exists for. Plant an unrecognised entry.
  assign("bogus_entry", TRUE, envir = rnaparallel:::.combat_clusters)
  on.exit(suppressWarnings(rm(list = "bogus_entry",
                              envir = rnaparallel:::.combat_clusters)), add = TRUE)
  expect_silent(combat_cluster_stop())
})

test_that("the foreach backend stays parallel after its first dispatch", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("doParallel")
  skip_if_not_installed("foreach")
  # registerDoSEQ() makes getDoParRegistered() TRUE, so a guard that tests only that
  # registered the pool on the first dispatch and ran every later one sequentially in the
  # parent with the cached cluster idle. Measured before the fix: 4 worker PIDs, then 1.
  on.exit(combat_cluster_stop(), add = TRUE)
  parent <- Sys.getpid()
  idx <- rnaparallel:::combat_row_chunks(200L, workers = 2L)
  pids_for <- function() unlist(rnaparallel:::combat_parallel_lapply(
    idx, function(i) Sys.getpid(), workers = 2L, parallel_backend = "foreach",
    cells = Inf, min_cells = 0))

  first  <- pids_for()
  second <- pids_for()
  third  <- pids_for()

  expect_false(all(first == parent))
  expect_false(all(second == parent))
  expect_false(all(third == parent))
  expect_gt(length(unique(second)), 1L)
})

test_that("a backend that namespace-qualifies a rebound call is refused", {
  # all.names() flattens edgeR::glmFit to include "glmFit", so a gate built on it would
  # accept exactly the degradation it exists to catch.
  mk <- function(qualified) {
    f <- function(counts, batch, group = NULL, covar_mod = NULL, full_mod = TRUE,
                  shrink = FALSE, shrink.disp = FALSE, gene.subset.n = NULL) NULL
    body(f) <- parse(text = sprintf(
      "{ match_quantiles(a, b, c, d, e); %s; glmFit.default(y); estimateGLMTagwiseDisp(y); sapply(z, g); lapply(z, g) }",
      if (qualified) "edgeR::glmFit(y)" else "glmFit(y)"))[[1]]
    e <- new.env(parent = globalenv())
    assign("match_quantiles",
           function(counts_sub, old_mu, old_phi, new_mu, new_phi) counts_sub, envir = e)
    environment(f) <- e
    f
  }
  expect_silent(rnaparallel:::combat_backend(mk(FALSE)))
  expect_error(rnaparallel:::combat_backend(mk(TRUE)), "no longer calls")
  # the weaker check really would have let it through
  expect_true("glmFit" %in% all.names(body(mk(TRUE))))
})

test_that("the glmFit row split is reached, and returns what the original returns", {
  # Both tests above pass group = NULL, which builds a batch-only design that
  # combat_design_oneway refuses to split, so glmFit contributes nothing to their counts.
  # A design with a condition column is the only way into this path.
  set.seed(11)
  ng <- 400L; ns <- 18L
  y <- matrix(rnbinom(ng * ns, mu = 50, size = 2), ng, ns,
              dimnames = list(paste0("g", seq_len(ng)), paste0("s", seq_len(ns))))
  bat <- rep(1:3, each = ns / 3)
  grp <- rep(0:1, times = ns / 2)
  des <- model.matrix(~ as.factor(bat) + as.factor(grp))
  expect_false(rnaparallel:::combat_design_oneway(des))

  # The shapes ComBat-seq actually passes, not the convenient ones. A scalar dispersion slips
  # through rp_per_gene() untouched and a row-constant offset makes any row permutation
  # arithmetically invisible, so both misalignment bugs this test exists to catch would pass
  # against them. Per-gene dispersion and a row-varying offset make the alignment observable.
  common <- edgeR::estimateGLMCommonDisp(y, design = des)
  disp <- matrix(rep(edgeR::estimateGLMTagwiseDisp(y, design = des, dispersion = common,
                                                   prior.df = 0), ns), ng, ns)
  set.seed(23)
  off <- matrix(log(colSums(y)), ng, ns, byrow = TRUE) + matrix(rnorm(ng, 0, 0.4), ng, ns)
  original <- edgeR::glmFit.default(y, design = des, dispersion = disp, offset = off,
                                  lib.size = NULL, weights = NULL, prior.count = 1e-04,
                                  start = NULL)

  # the five fields the parent binds and ComBat-seq reads. `failed` is deliberately absent
  # from the assembled result: the parent consumes it to decide whether to discard the split,
  # and sva::ComBat_seq never reads it.
  keep <- c("coefficients", "fitted.values", "df.residual", "unshrunk.coefficients", "method")
  for (w in c(2L, 4L)) for (k in c(2L, 3L, 7L)) {
    p <- rnaparallel:::glmFit_rows_parallel(y, design = des, dispersion = disp, offset = off,
                                            prior.count = 1e-04, workers = w, chunks = k)
    # the split has two escapes back to the whole-matrix original fit, the one-way gate and the
    # failed-fit fallback, and either would make every assertion below pass trivially. The
    # split assembles a plain list; both escapes return edgeR's DGEGLM.
    expect_false(inherits(p, "DGEGLM"), info = sprintf("split not taken, w=%d k=%d", w, k))
    for (nm in keep) expect_identical(p[[nm]], original[[nm]], info = sprintf("%s w=%d k=%d", nm, w, k))
  }
})

test_that("a condition column adds glmFit dispatches a batch-only design never makes", {
  # guards the parallel layer going dead: identical() alone cannot tell a working split from
  # one that quietly calls the original whole, because ComBat-seq returns whole-number counts
  # and rounds small differences away. A spy backend rather than trace(): the counted function
  # forks, and a traced function inherited by a forked child wedges the worker.
  skip_if_not_installed("sva")
  set.seed(7)
  ng <- 300L; ns <- 24L
  cts <- matrix(rnbinom(ng * ns, mu = 40, size = 1.5), ng, ns,
                dimnames = list(paste0("g", seq_len(ng)), paste0("s", seq_len(ns))))
  bat <- rep(1:3, each = ns / 3)
  grp <- rep(0:1, times = ns / 2)

  count_for <- function(group) {
    calls <- 0L
    spy <- function(idx, f, workers) { calls <<- calls + 1L; lapply(idx, f) }
    invisible(quietly(ComBat_seq_parallel(cts, batch = bat, group = group, workers = 2L,
                                         parallel_backend = spy)))
    calls
  }
  # the batch-only design is refused by combat_design_oneway, so glmFit contributes nothing.
  # adding a condition column makes the design two-way and the row split runs, twice.
  expect_identical(count_for(grp) - count_for(NULL), 2L)
})
