## Synthetic count matrices for the test suite.
##
## `adversarial` exists because a generator built to STRESS the algorithm and one
## built to DEMONSTRATE it working are different jobs. Genes zeroed in one batch or
## pinned at a constant carry batch signal that is structurally uncorrectable, so
## leaving them in a before/after summary makes a correct correction look inert:
## batch R-squared on PC1 barely moved, 0.997 to 0.994, while clean data on the same
## code gave 0.9935 to 0.0002. Equivalence tests want them ON, since that is where
## a hand-transcribed copy diverges. Visual checks want them OFF.
make_counts <- function(seed, G, n_per_batch, with_group = FALSE,
                        batch_effect = 0.7, group_effect = 0.5, adversarial = TRUE) {
  set.seed(seed)
  batch <- factor(rep(paste0("b", seq_along(n_per_batch)), times = n_per_batch))
  n <- length(batch)
  grp <- if (with_group) factor(rep_len(c("ctl", "trt"), n)) else NULL

  mu_g <- 2^runif(G, 0, 11)
  mu_mat <- outer(mu_g, rep(1, n))
  ## the signal ComBat-seq exists to remove: each gene shifted by a different amount in each batch
  if (batch_effect > 0)
    mu_mat <- mu_mat * matrix(exp(rnorm(G * nlevels(batch), 0, batch_effect)),
                              G, nlevels(batch))[, as.integer(batch)]
  ## biology that must survive the correction, on a tenth of the genes
  if (!is.null(grp) && group_effect > 0) {
    de <- seq_len(max(1L, G %/% 10L))
    mu_mat[de, grp == "trt"] <- mu_mat[de, grp == "trt"] * exp(rnorm(length(de), 0, group_effect))
  }
  cts <- matrix(rnbinom(G * n, mu = as.vector(mu_mat), size = rep(1 / runif(G, 0.05, 2), n)), G, n)

  if (adversarial) {
    b2 <- batch == levels(batch)[2]
    b3 <- batch == levels(batch)[nlevels(batch)]
    k <- max(1L, G %/% 30L)
    cts[seq_len(k), b2] <- 0L                                    # all zero in exactly one batch
    cts[k + seq_len(k), ] <- 0L                                  # dead everywhere
    cts[2 * k + seq_len(k), b3] <- sample(0:1, k * sum(b3), TRUE) # all 0 or 1 in one batch
    cts[3 * k + seq_len(k), b2] <- 5L                             # zero variance, nonzero
    cts[4 * k + seq_len(2 * k), ] <- rbinom(2 * k * n, 3, 0.05)   # sparse tail
  }

  storage.mode(cts) <- "integer"
  dimnames(cts) <- list(paste0("g", seq_len(G)), paste0("s", seq_len(n)))
  list(counts = cts, batch = batch, group = grp)
}

## The backend, resolved the same way the package resolves it, so a test never
## compares against a different copy of ComBat-seq than the companion used.
backend_fn <- function() {
  testthat::skip_if_not_installed("sva")
  sva::ComBat_seq
}

edgeR_norm <- function(object, ...) {
  ns <- asNamespace("edgeR")
  generic <- if (exists("normLibSizes.default", envir = ns, inherits = FALSE)) {
    "normLibSizes"
  } else {
    "calcNormFactors"
  }
  get(generic, envir = ns, inherits = FALSE)(object, ...)
}

quietly <- function(expr) {
  out <- NULL
  invisible(utils::capture.output(out <- suppressMessages(expr)))
  out
}

## R CMD check sets _R_CHECK_LIMIT_CORES_ and then errors above two cores. The package
## honours that itself, but a test that builds its OWN executor bypasses the package's
## clamp, so the custom executors below have to honour it too. This only started
## mattering once setup-parallel.R made the suite actually fork.
test_max_cores <- function() {
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && !identical(tolower(chk), "false")) 2L else max(1L, parallel::detectCores(), na.rm = TRUE)
}
