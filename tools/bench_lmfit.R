#!/usr/bin/env Rscript
## Times limma::lmFit against lmFit_parallel() on a simulated voom EList of 60,000 genes by 48 samples, unblocked and with block + correlation, at the package's default worker count; the source of the 3.70x and 4.35x in REFERENCE.md and the lmFit_parallel() help page, run as `Rscript tools/bench_lmfit.R` with rnaparallel installed.
suppressMessages({library(rnaparallel); library(edgeR); library(limma)})
options(combat.progress = FALSE)
set.seed(20261004)
G <- 60000L; N <- 48L
mu <- exp(rnorm(G, mean = 4, sd = 1.8))
counts <- matrix(rnbinom(G * N, mu = mu, size = 1 / 0.15), G, N)
group <- factor(rep(c("ctl", "trt"), length.out = N))
subj <- factor(rep(seq_len(N / 2L), each = 2L))
design <- model.matrix(~ group)
v <- voom(DGEList(counts), design)
cor <- duplicateCorrelation(v, design, block = subj)$consensus.correlation
best <- function(f) min(vapply(1:3, function(i) system.time(f())[["elapsed"]], numeric(1)))
w <- rnaparallel:::combat_default_workers(NULL)
res <- list()
for (arm in c("unblocked", "blocked")) {
  orig <- if (arm == "unblocked") function() lmFit(v, design) else
    function() lmFit(v, design, block = subj, correlation = cor)
  comp <- if (arm == "unblocked") function() lmFit_parallel(v, design, workers = w) else
    function() lmFit_parallel(v, design, block = subj, correlation = cor, workers = w)
  stopifnot(identical(orig(), comp()))
  t0 <- best(orig); t1 <- best(comp)
  res[[arm]] <- data.frame(arm = arm, workers = w, original_s = t0, companion_s = t1,
                           speedup = round(t0 / t1, 2))
}
print(do.call(rbind, res), row.names = FALSE)
cat("machine:", Sys.info()[["machine"]], " cores:", parallel::detectCores(), " R:", R.version.string,
    " limma:", as.character(packageVersion("limma")), "\n")
