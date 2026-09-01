## The identical() proof.
##
## That ComBat_seq_parallel returns exactly what the original ComBat-seq returns, on every
## argument path and every chunk layout, and that it delegates rather than reimplements.
## If this file passes, the package's central claim holds.

## identical(), not all.equal(). The whole promise of the rebind is that there is
## no second copy of the algorithm, so any difference at all is a defect.

test_that("default path is identical to a serial run", {
  original <- backend_fn()
  d <- make_counts(1, G = 400, n_per_batch = c(8, 7, 6))
  ref <- quietly(original(d$counts, d$batch, group = NULL))
  par <- quietly(ComBat_seq_parallel(d$counts, d$batch, group = NULL, workers = 2L))
  expect_identical(par, ref)
})

test_that("identical with a group, which changes the design and the dispersion path", {
  original <- backend_fn()
  d <- make_counts(2, G = 400, n_per_batch = c(9, 8), with_group = TRUE)
  ref <- quietly(original(d$counts, d$batch, group = d$group))
  par <- quietly(ComBat_seq_parallel(d$counts, d$batch, group = d$group, workers = 2L))
  expect_identical(par, ref)
})

test_that("identical on the argument paths", {
  original <- backend_fn()
  d <- make_counts(3, G = 300, n_per_batch = c(8, 8), with_group = TRUE)

  # full_mod = FALSE
  expect_identical(
    quietly(ComBat_seq_parallel(d$counts, d$batch, group = d$group, full_mod = FALSE, workers = 2L)),
    quietly(original(d$counts, d$batch, group = d$group, full_mod = FALSE)))

  # shrink. seed before EACH arm: the original's monte_carlo_int_NB calls sample(),
  # so unseeded back-to-back runs draw different gene subsets and differ for a
  # reason that is not a bug. An earlier comparison reported a fake max|diff| of
  # 1096 from exactly this.
  set.seed(99); ref <- quietly(original(d$counts, d$batch, group = d$group,
                                      shrink = TRUE, shrink.disp = TRUE))
  set.seed(99); par <- quietly(ComBat_seq_parallel(d$counts, d$batch, group = d$group,
                                             shrink = TRUE, shrink.disp = TRUE, workers = 2L))
  expect_identical(par, ref)

  set.seed(7); ref2 <- quietly(original(d$counts, d$batch, group = d$group,
                                      shrink = TRUE, gene.subset.n = 50))
  set.seed(7); par2 <- quietly(ComBat_seq_parallel(d$counts, d$batch, group = d$group,
                                              shrink = TRUE, gene.subset.n = 50, workers = 2L))
  expect_identical(par2, ref2)
})

test_that("covar_mod is identical, upstream's own example", {
  # The example from the ComBat-seq README, which is how this actually gets called.
  original <- backend_fn()
  set.seed(42)
  cov1 <- rep(c(0, 1), 4)
  cov2 <- c(0, 0, 1, 1, 0, 0, 1, 1)
  covar_mat <- cbind(cov1, cov2)
  n <- nrow(covar_mat)
  batch <- rep(1:2, each = 4)
  cts <- matrix(rnbinom(300 * n, mu = 60, size = 4), nrow = 300,
                dimnames = list(paste0("g", 1:300), paste0("s", 1:n)))

  # two columns, and one column: neither collapses on the way through
  for (cm in list(covar_mat, cbind(cov1))) {
    expect_identical(
      quietly(ComBat_seq_parallel(cts, batch = batch, group = NULL,
                                  covar_mod = cm, workers = 2L)),
      quietly(original(cts, batch = batch, group = NULL, covar_mod = cm)))
  }
})

test_that("a confounded covar_mod is rejected the same way by both", {
  # cov1 alternates exactly like a 2-level group, so group + cov1 is collinear and
  # ComBat-seq refuses. Correct behaviour, and the companion must refuse identically.
  original <- backend_fn()
  set.seed(43)
  cov1 <- rep(c(0, 1), 4)
  covar_mat <- cbind(cov1, cov2 = c(0, 0, 1, 1, 0, 0, 1, 1))
  batch <- rep(1:2, each = 4)
  grp <- factor(rep(c("a", "b"), 4))
  cts <- matrix(rnbinom(300 * 8, mu = 60, size = 4), nrow = 300)

  ref <- tryCatch(quietly(original(cts, batch, group = grp, covar_mod = covar_mat)),
                  error = function(e) conditionMessage(e))
  par <- tryCatch(quietly(ComBat_seq_parallel(cts, batch, group = grp,
                                              covar_mod = covar_mat, workers = 2L)),
                  error = function(e) conditionMessage(e))
  expect_match(ref, "confounded")
  expect_identical(par, ref)
})

test_that("the serial escape hatch gives the same answer", {
  original <- backend_fn()
  d <- make_counts(4, G = 250, n_per_batch = c(6, 6))
  ref <- quietly(original(d$counts, d$batch, group = NULL))
  withr_fork <- getOption("combat.fork")
  options(combat.fork = FALSE)
  on.exit(options(combat.fork = withr_fork), add = TRUE)
  expect_identical(quietly(ComBat_seq_parallel(d$counts, d$batch, group = NULL, workers = 4L)), ref)
})

test_that("the input matrix is not modified by either arm", {
  original <- backend_fn()
  d <- make_counts(5, G = 200, n_per_batch = c(6, 6))
  # + 0L forces a real copy. Aliasing the same SEXP would make a C-level in-place write
  # change both sides, so the comparison would pass while the property it checks is false.
  before <- d$counts + 0L
  invisible(quietly(original(d$counts, d$batch, group = NULL)))
  invisible(quietly(ComBat_seq_parallel(d$counts, d$batch, group = NULL, workers = 2L)))
  expect_identical(d$counts, before)
})

test_that("batch given as character and as integer both work", {
  d <- make_counts(6, G = 200, n_per_batch = c(6, 6))
  as_chr <- as.character(d$batch)
  as_int <- as.integer(d$batch)
  a <- quietly(ComBat_seq_parallel(d$counts, as_chr, group = NULL, workers = 2L))
  b <- quietly(ComBat_seq_parallel(d$counts, as_int, group = NULL, workers = 2L))
  expect_identical(a, b)
})

test_that("an unused batch level breaks the original, and breaks us the same way", {
  # Worth knowing before it bites a real cohort: ComBat_seq counts samples per
  # FACTOR LEVEL, so a level left over from subsetting looks like a batch with no
  # samples and the original stops with "doesn't support 1 sample per batch yet".
  # Nothing in that message points at the unused level. droplevels() first.
  original <- backend_fn()
  d <- make_counts(7, G = 200, n_per_batch = c(6, 6))
  padded <- factor(as.character(d$batch), levels = c(levels(d$batch), "b_unused"))

  ref <- tryCatch(quietly(original(d$counts, padded, group = NULL)),
                  error = function(e) conditionMessage(e))
  par <- tryCatch(quietly(ComBat_seq_parallel(d$counts, padded, group = NULL, workers = 2L)),
                  error = function(e) conditionMessage(e))
  expect_identical(par, ref)

  # and dropping the level makes both work and agree
  dropped <- droplevels(padded)
  expect_identical(quietly(ComBat_seq_parallel(d$counts, dropped, group = NULL, workers = 2L)),
                   quietly(original(d$counts, dropped, group = NULL)))
})
## The anti-regression suite. An earlier version of this code transcribed
## ComBat-seq by hand and drifted from it twice in ways synthetic data could not
## expose. These tests fail if anyone turns the companion back into a copy.

test_that("every ComBat_seq argument exists on the companion with the same default", {
  original <- backend_fn()
  vf <- formals(original)
  pf <- formals(ComBat_seq_parallel)

  expect_true(all(names(vf) %in% names(pf)))
  for (a in names(vf)) {
    expect_identical(deparse(pf[[a]]), deparse(vf[[a]]),
                     info = paste("default differs for argument", a))
  }
})

test_that("the companion adds only parallel controls, nothing else", {
  expect_identical(setdiff(names(formals(ComBat_seq_parallel)), names(formals(backend_fn()))),
                   c("workers", "chunks", "parallel_backend", "backend", "label"))
})

test_that("argument order matches the original for the shared arguments", {
  original <- backend_fn()
  shared <- intersect(names(formals(ComBat_seq_parallel)), names(formals(original)))
  expect_identical(shared, names(formals(original)))
})

test_that("the companion does not reimplement the algorithm", {
  src <- paste(deparse(ComBat_seq_parallel), collapse = "\n")
  # Match CALLS, not bare mentions. `inherits(y, "DGEList")` is a class check and
  # entirely legitimate here; only `DGEList(` would mean the companion is building
  # the object itself. A looser grep flagged the class check as a reimplementation.
  vendor_calls <- c("estimateGLMCommonDisp", "estimateGLMTagwiseDisp", "vec2mat",
                    "monte_carlo_int_NB", "DGEList", "mapDisp")
  for (sym in vendor_calls) {
    expect_false(grepl(paste0(sym, "("), src, fixed = TRUE),
                 info = paste("companion body CALLS original internal", sym))
  }
  # and the rebinding mechanism must still be the mechanism
  expect_true(grepl("environment(f) <- env", src, fixed = TRUE))
  expect_true(grepl("new.env(parent = be$env)", src, fixed = TRUE))
})

test_that("no ::: appears anywhere in the package sources", {
  # An installed package ships R/rnaparallel.rdb, not R/*.R, so reading every file
  # in that directory feeds binary to readLines and only produces locale warnings.
  # Deparse the closures instead, which is what actually got installed.
  ns <- asNamespace("rnaparallel")
  src <- unlist(lapply(ls(ns, all.names = TRUE), function(n) {
    obj <- get(n, envir = ns)
    if (is.function(obj)) deparse(obj) else character(0)
  }))
  expect_gt(length(src), 0)
  expect_false(any(grepl(":::", src, fixed = TRUE)))
})

test_that("the function actually run is the original's, byte for byte", {
  original <- backend_fn()
  be <- rnaparallel:::combat_backend()
  expect_identical(body(be$fn), body(original))
  expect_identical(formals(be$fn), formals(original))
})

test_that("the export surface is exactly the entry points and the shared controls", {
  # Pinned deliberately. A companion that is added without a measured speedup, or a helper
  # that leaks out of the namespace, changes this line first.
  expect_setequal(getNamespaceExports("rnaparallel"),
                  c("ComBat_seq_parallel",
                    "lmFit_parallel", "duplicateCorrelation_parallel",
                    "calcNormFactors_parallel",
                    "combat_backends", "combat_cluster_stop", "rnaparallel_stale",
                    "rnaparallel_progress",
                    "removeBatchEffect_parallel"))
})
## The backend resolver is what makes this work against sva OR a sourced copy of
## upstream zhangyuqing/ComBat-seq, and what removes the need for a ::: call.

test_that("the default backend resolves with the helper from its own environment", {
  skip_if_not_installed("sva")
  be <- rnaparallel:::combat_backend()
  expect_true(is.function(be$fn))
  expect_true(is.function(be$match_quantiles))
  expect_identical(names(formals(be$match_quantiles)),
                   c("counts_sub", "old_mu", "old_phi", "new_mu", "new_phi"))
  # the helper must come from the same place as the function, never a second copy
  expect_identical(be$env, environment(be$fn))
})

test_that("a sourced upstream copy works, helpers in a plain environment", {
  skip_if_not_installed("sva")
  up <- new.env(parent = globalenv())
  # borrow the real original pieces but present them the way source() would
  assign("match_quantiles", get("match_quantiles", envir = asNamespace("sva")), envir = up)
  fake <- sva::ComBat_seq
  environment(fake) <- up

  be <- rnaparallel:::combat_backend(fake)
  expect_identical(be$env, up)
  expect_false(isNamespace(be$env))
  expect_true(is.function(be$match_quantiles))
})

test_that("a backend missing an argument is refused, not called", {
  bad <- function(counts, batch, group = NULL) NULL
  expect_error(rnaparallel:::combat_backend(bad), "missing argument")
})

test_that("a renamed or reshaped match_quantiles is refused", {
  env <- new.env(parent = globalenv())
  assign("match_quantiles", function(a, b) NULL, envir = env)
  f <- function(counts, batch, group = NULL, covar_mod = NULL, full_mod = TRUE,
                shrink = FALSE, shrink.disp = FALSE, gene.subset.n = NULL) NULL
  environment(f) <- env
  expect_error(rnaparallel:::combat_backend(f), "splits it by row")
})

test_that("a backend with no match_quantiles in scope is refused", {
  env <- new.env(parent = baseenv())
  f <- function(counts, batch, group = NULL, covar_mod = NULL, full_mod = TRUE,
                shrink = FALSE, shrink.disp = FALSE, gene.subset.n = NULL) NULL
  environment(f) <- env
  expect_error(rnaparallel:::combat_backend(f), "could not find")
})

test_that("a non-function backend is refused", {
  expect_error(rnaparallel:::combat_backend("sva::ComBat_seq"), "must be a function")
})
## Tagwise dispersion is the third parallelised path, and the only one whose exactness
## depends on an argument value rather than on the shape of the computation. It is exact
## because ComBat-seq passes prior.df = 0, which zeroes the weight on the across-gene
## moderation term. These tests fail if that reasoning stops holding.

test_that("splitting tagwise dispersion by row is exact at prior.df = 0", {
  skip_if_not_installed("edgeR")
  set.seed(31)
  G <- 1200L; n <- 20L
  y <- matrix(rnbinom(G * n, mu = rep(2^runif(G, 0, 11), n), size = 4), G, n)
  storage.mode(y) <- "integer"
  # A continuous column, so the design is not the one-group layout the split refuses. An
  # intercept-only design here would make this test assert 0 dispatches twice and prove nothing
  # about the gate it exists to check.
  design <- cbind(Intercept = 1, x = as.numeric(scale(seq_len(n))))
  dc <- edgeR::estimateGLMCommonDisp(y, design = design, subset = nrow(y))

  ref <- edgeR::estimateGLMTagwiseDisp(y, design = design, dispersion = dc, prior.df = 0)
  got <- rnaparallel:::estimateGLMTagwiseDisp_rows_parallel(
    y, design = design, dispersion = dc, prior.df = 0,
    workers = 2L, chunks = 4L, parallel_backend = "serial")

  expect_identical(got, ref)
})

test_that("it stays exact when dead, constant and near-empty genes are present", {
  # the failure mode that synthetic uniform counts cannot expose: a gene whose
  # adjusted profile likelihood is degenerate would poison an across-gene mean, and
  # the chunk mean and the whole-matrix mean are not the same set of genes
  skip_if_not_installed("edgeR")
  set.seed(32)
  G <- 1500L; n <- 24L
  y <- matrix(rnbinom(G * n, mu = rep(2^runif(G, 0, 11), n), size = 4), G, n)
  k <- G %/% 30L
  y[seq_len(k), ] <- 0L
  y[k + seq_len(k), ] <- 5L
  y[2 * k + seq_len(k), ] <- rbinom(k * n, 3, 0.02)
  y[3 * k + seq_len(k), seq_len(n %/% 2)] <- 0L
  storage.mode(y) <- "integer"
  # A continuous column, so the design is not the one-group layout the split refuses. An
  # intercept-only design here would make this test assert 0 dispatches twice and prove nothing
  # about the gate it exists to check.
  design <- cbind(Intercept = 1, x = as.numeric(scale(seq_len(n))))
  dc <- edgeR::estimateGLMCommonDisp(y, design = design, subset = nrow(y))

  ref <- edgeR::estimateGLMTagwiseDisp(y, design = design, dispersion = dc, prior.df = 0)
  for (nch in c(2L, 4L, 7L)) {
    got <- rnaparallel:::estimateGLMTagwiseDisp_rows_parallel(
      y, design = design, dispersion = dc, prior.df = 0,
      workers = 2L, chunks = nch, parallel_backend = "serial")
    expect_identical(got, ref, info = paste("chunks:", nch))
  }
})

test_that("any prior.df other than 0 is handed to edgeR whole, not split", {
  # moderation towards the across-gene curve is a real cross-row term, so splitting
  # would silently change every gene. The guard must refuse rather than approximate.
  skip_if_not_installed("edgeR")
  set.seed(33)
  G <- 600L; n <- 16L
  y <- matrix(rnbinom(G * n, mu = rep(2^runif(G, 0, 11), n), size = 4), G, n)
  storage.mode(y) <- "integer"
  # A continuous column, so the design is not the one-group layout the split refuses. An
  # intercept-only design here would make this test assert 0 dispatches twice and prove nothing
  # about the gate it exists to check.
  design <- cbind(Intercept = 1, x = as.numeric(scale(seq_len(n))))
  dc <- edgeR::estimateGLMCommonDisp(y, design = design, subset = nrow(y))

  for (pdf in c(10, 5)) {
    ref <- edgeR::estimateGLMTagwiseDisp(y, design = design, dispersion = dc, prior.df = pdf)
    got <- rnaparallel:::estimateGLMTagwiseDisp_rows_parallel(
      y, design = design, dispersion = dc, prior.df = pdf,
      workers = 2L, chunks = 4L, parallel_backend = "serial")
    expect_identical(got, ref, info = paste("prior.df:", pdf))
  }
})

test_that("the end to end result is unchanged now dispersion is parallel too", {
  original <- backend_fn()
  d <- make_counts(34, G = 900L, n_per_batch = c(14, 12, 13, 11))
  ref <- quietly(original(d$counts, d$batch, group = NULL))

  for (w in c(1L, 2L, 4L)) {
    got <- quietly(ComBat_seq_parallel(d$counts, d$batch, group = NULL, workers = w))
    expect_identical(got, ref, info = paste("workers:", w))
  }
})

test_that("the dispersion path has its own, higher size threshold", {
  # One shared threshold made this split a net loss: dispersion costs less per cell than
  # match_quantiles, so 20,000 cells pays for a fork on one path and not the other.
  # Measured on 10-column slices: 0.65x at 10,000 cells, 0.79x at 20,000, 1.08x at 30,000.
  calls <- 0L
  spy <- function(idx, f, workers) { calls <<- calls + 1L; lapply(idx, f) }
  set.seed(35)
  G <- 400L; n <- 12L
  y <- matrix(rnbinom(G * n, mu = rep(2^runif(G, 0, 11), n), size = 4), G, n)
  storage.mode(y) <- "integer"
  # A continuous column, so the design is not the one-group layout the split refuses. An
  # intercept-only design here would make this test assert 0 dispatches twice and prove nothing
  # about the gate it exists to check.
  design <- cbind(Intercept = 1, x = as.numeric(scale(seq_len(n))))
  dc <- edgeR::estimateGLMCommonDisp(y, design = design, subset = nrow(y))

  old <- options(combat.min.cells = 0, combat.min.disp.cells = 5e4)
  on.exit(options(old), add = TRUE)

  # 4,800 cells: under the dispersion threshold even though combat.min.cells allows it
  rnaparallel:::estimateGLMTagwiseDisp_rows_parallel(
    y, design = design, dispersion = dc, prior.df = 0,
    workers = 2L, chunks = 2L, parallel_backend = spy)
  expect_identical(calls, 0L)

  options(combat.min.disp.cells = 0)
  rnaparallel:::estimateGLMTagwiseDisp_rows_parallel(
    y, design = design, dispersion = dc, prior.df = 0,
    workers = 2L, chunks = 2L, parallel_backend = spy)
  expect_identical(calls, 1L)
})

test_that("a one-group design is refused the tagwise split even under the gate", {
  # The gate is what keeps ComBat-seq exact, so it needs its own assertion rather than being
  # inferred from the dispatch counts elsewhere.
  skip_if_not_installed("edgeR")
  set.seed(77)
  y <- matrix(rnbinom(60 * 8, mu = 40, size = 3), 60, 8)
  calls <- 0L
  spy <- function(idx, f, w) { calls <<- calls + 1L; lapply(idx, f) }
  withr::local_options(combat.min.disp.cells = 0)
  rnaparallel:::estimateGLMTagwiseDisp_rows_parallel(
    y, design = cbind(rep(1, 8)), dispersion = 0.1, prior.df = 0,
    workers = 2L, chunks = 2L, parallel_backend = spy)
  expect_identical(calls, 0L)
})

## Every other fixture in this suite uses benign counts with near-equal library sizes, and that
## is exactly why they all passed while ComBat_seq_parallel was returning chunk-dependent output.
## Low counts make edgeR's fits fail to converge; one over-sequenced library makes the offsets
## extreme. Together they reach the one-group kernel, whose result depends on which other genes
## share the matrix. Measured before the gate: 1, 2 and 3 differing cells at chunks 2, 4 and 8.

oversequenced <- function(seed = 2024L, G = 240L, S = 12L, mult = 1000L) {
  set.seed(seed)
  y <- matrix(stats::rnbinom(G * S, mu = 2, size = 0.3), G, S,
              dimnames = list(paste0("g", seq_len(G)), paste0("s", seq_len(S))))
  y[, 5] <- as.integer(y[, 5] * mult)
  y
}

test_that("an over-sequenced library with non-convergent fits stays identical", {
  skip_if_not_installed("sva")
  batch <- rep(1:2, each = 6)
  for (sd in c(2024L, 1L, 42L)) {
    y <- oversequenced(sd)
    ref <- suppressMessages(sva::ComBat_seq(y, batch = batch, group = NULL))
    for (k in c(2L, 4L, 8L)) {
      expect_identical(
        suppressMessages(ComBat_seq_parallel(y, batch = batch, group = NULL,
                                             workers = 4L, chunks = k)),
        ref, info = paste("seed", sd, "chunks", k))
    }
  }
})

test_that("the same fixture with a group supplied stays identical", {
  skip_if_not_installed("sva")
  y <- oversequenced()
  batch <- rep(1:2, each = 6); group <- rep(rep(1:2, each = 3), 2)
  ref <- suppressMessages(sva::ComBat_seq(y, batch = batch, group = group))
  for (k in c(2L, 4L, 8L)) {
    expect_identical(
      suppressMessages(ComBat_seq_parallel(y, batch = batch, group = group,
                                           workers = 4L, chunks = k)),
      ref, info = paste("chunks", k))
  }
})

test_that("a one-group design is refused the row split rather than trusted to it", {
  # The gate is design-level because the bad values are finite and carry no per-gene flag,
  # so no result-level filter can see them.
  expect_true(combat_design_oneway(NULL))
  expect_true(combat_design_oneway(matrix(1, 6, 1)))
  expect_true(combat_design_oneway(stats::model.matrix(~ 0 + factor(rep(1:2, each = 3)))))
  expect_false(combat_design_oneway(
    cbind(1, rep(0:1, each = 3), stats::rnorm(6))))
})

## The one copy of original code this package holds is the row-vectorised match_quantiles.
## The gate that guards it had no coverage, so a permanently shut gate would leave every
## suite green while silently running the slow path forever.

test_that("the match_quantiles gate opens on the real backend, even under hostile options", {
  skip_if_not_installed("sva")
  set.seed(9)
  cs <- matrix(rnbinom(200 * 8, mu = 50, size = 2), 200, 8)
  om <- matrix(runif(200 * 8, 20, 80), 200, 8)
  op <- runif(200, 0.05, 0.5)
  expect_identical(
    rnaparallel:::combat_mq_dispatch(sva:::match_quantiles, cs, om, op),
    rnaparallel:::match_quantiles_rows)
  # deparse honours scipen; an analyst's .Rprofile must not silently shut the gate
  withr::local_options(scipen = 999)
  expect_identical(
    rnaparallel:::combat_mq_dispatch(sva:::match_quantiles, cs, om, op),
    rnaparallel:::match_quantiles_rows)
})

test_that("the match_quantiles gate closes on a genuinely edited original body", {
  skip_if_not_installed("sva")
  set.seed(9)
  cs <- matrix(rnbinom(50 * 4, mu = 50, size = 2), 50, 4)
  om <- matrix(runif(50 * 4, 20, 80), 50, 4)
  op <- runif(50, 0.05, 0.5)
  f <- sva:::match_quantiles
  body(f)[[2]] <- quote(new_counts_sub <- matrix(NA_integer_, nrow = 0, ncol = 0))
  expect_identical(rnaparallel:::combat_mq_dispatch(f, cs, om, op), f)
})

test_that("the row-vectorised transcription matches the original cell loop bit for bit", {
  skip_if_not_installed("sva")
  # adversarial: zeros, ones, huge counts, and cells landing in the outlier branch
  set.seed(41)
  cs <- matrix(rnbinom(300 * 10, mu = 3, size = 0.3), 300, 10)
  cs[1:40, ] <- 0L; cs[41:60, 1] <- 1L; cs[61:70, ] <- 50000L
  om <- matrix(runif(300 * 10, 0.5, 200), 300, 10)
  op <- runif(300, 0.01, 4)
  nm <- matrix(runif(300 * 10, 0.5, 200), 300, 10)
  np <- runif(300, 0.01, 4)
  expect_identical(rnaparallel:::match_quantiles_rows(cs, om, op, nm, np),
                   sva:::match_quantiles(cs, om, op, nm, np))
})

test_that("an NA count reproduces the original's own error, not a wrapped one", {
  skip_if_not_installed("sva")
  set.seed(9)
  cs <- matrix(rnbinom(40 * 4, mu = 50, size = 2), 40, 4); cs[3, 2] <- NA
  om <- matrix(runif(40 * 4, 20, 80), 40, 4)
  op <- runif(40, 0.05, 0.5)
  vendor_err <- tryCatch(sva:::match_quantiles(cs, om, op, om, op),
                         error = conditionMessage)
  par_err <- tryCatch(
    rnaparallel:::match_quantiles_parallel(sva:::match_quantiles, cs, om, op, om, op,
                                           workers = 2L),
    error = conditionMessage)
  expect_identical(par_err, vendor_err)
})

test_that("data.frame counts go to the original whole instead of crashing the gate", {
  skip_if_not_installed("sva")
  # is.finite has no data.frame method; the original accepts data.frames end to end
  set.seed(6)
  y <- matrix(rnbinom(200 * 8, mu = 60, size = 5), 200, 8,
              dimnames = list(paste0("g", 1:200), paste0("s", 1:8)))
  df <- as.data.frame(y); batch <- rep(1:2, each = 4)
  ref <- suppressMessages(sva::ComBat_seq(df, batch = batch, group = NULL))
  expect_identical(
    suppressMessages(ComBat_seq_parallel(df, batch = batch, group = NULL, workers = 2L)),
    ref)
})

test_that("the tagwise batch dispatch survives the original's own gene filter", {
  skip_if_not_installed("sva")
  # ComBat-seq drops genes all-zero within a batch before the tagwise lapply, so a shape
  # check against the unfiltered row count failed every batch and the dominant stage was
  # silently recomputed serially.
  set.seed(7)
  y <- matrix(rnbinom(300 * 12, mu = 40, size = 4), 300, 12,
              dimnames = list(paste0("g", 1:300), paste0("s", 1:12)))
  batch <- rep(1:2, each = 6)
  y[5, batch == 1] <- 0L                      # filtered by the original
  calls <- 0L
  spy <- function(idx, f, w) { calls <<- calls + 1L; lapply(idx, f) }
  out <- suppressMessages(ComBat_seq_parallel(y, batch = batch, group = NULL,
                                              workers = 2L, parallel_backend = spy))
  expect_identical(out, suppressMessages(sva::ComBat_seq(y, batch = batch, group = NULL)))
  # one match_quantiles per batch, one common dispersion, one tagwise across batches: the
  # tagwise dispatch must still be among them rather than falling back to base::lapply
  expect_identical(calls, 4L)
})

test_that("group= is identical across many batches without confounding", {
  skip_if_not_installed("sva")
  original <- backend_fn()

  # The pooled pan-cancer arm passes several hundred batches WITH group= set to cancer type,
  # and that combination had no test. Building the fixture by recycling batch and group
  # independently is what fails: sva refuses with "At least one covariate is confounded with
  # batch" whenever some batch does not see more than one group level. Laying the group cycle
  # INSIDE each batch guarantees every batch carries every group, so the design keeps full
  # rank and the original accepts it.
  set.seed(11)
  B <- 12L; G <- 3L; per <- 6L                       # per batch: 6 samples, 2 of each group
  N <- B * per
  batch <- rep(seq_len(B), each = per)
  group <- rep(rep(seq_len(G), each = per / G), times = B)
  expect_true(all(table(batch, group) > 0L))         # the property sva actually requires

  cts <- matrix(rnbinom(400L * N, mu = 50, size = 5), 400L, N,
                dimnames = list(paste0("g", 1:400), paste0("s", seq_len(N))))
  expect_identical(
    quietly(ComBat_seq_parallel(cts, batch = batch, group = group, workers = 2L)),
    quietly(original(cts, batch = batch, group = group)))
})

test_that("group= with a covariate matrix is identical too", {
  skip_if_not_installed("sva")
  original <- backend_fn()
  set.seed(12)
  B <- 8L; G <- 2L; per <- 6L
  N <- B * per
  batch <- rep(seq_len(B), each = per)
  group <- rep(rep(seq_len(G), each = per / G), times = B)
  covar <- cbind(stage = rep(rep(1:3, each = per / 3), times = B))
  cts <- matrix(rnbinom(300L * N, mu = 50, size = 5), 300L, N,
                dimnames = list(paste0("g", 1:300), paste0("s", seq_len(N))))
  expect_identical(
    quietly(ComBat_seq_parallel(cts, batch = batch, group = group,
                                covar_mod = covar, workers = 2L)),
    quietly(original(cts, batch = batch, group = group, covar_mod = covar)))
})
