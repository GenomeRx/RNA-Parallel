#!/usr/bin/env Rscript
## Serial time of every step in three RNA-seq pipelines and two batch corrections on one simulated cohort, the table behind README "Why these five", run as `Rscript tools/profile_pipelines.R [genes] [samples]` (default 20000 300, originals only, each at its own defaults except that topTable and topTags return the full table).
suppressMessages({library(edgeR); library(limma); library(DESeq2); library(sva)})
args <- commandArgs(trailingOnly = TRUE)
G <- if (length(args) >= 1L) suppressWarnings(as.integer(args[1L])) else 20000L
N <- if (length(args) >= 2L) suppressWarnings(as.integer(args[2L])) else 300L
if (is.na(G) || G < 1L || is.na(N) || N < 8L || N %% 4L != 0L)
  stop("samples must be a multiple of 4 and at least 8, because each subject gives two samples ",
       "and subjects alternate between the two groups; genes must be a positive integer", call. = FALSE)
set.seed(20261003)
mu0 <- exp(rnorm(G, mean = 4, sd = 1.8))
batch <- factor(rep(c("b1", "b2", "b3"), length.out = N))
group <- factor(rep(c("ctl", "trt"), each = 2L, length.out = N))
subj <- factor(rep(seq_len(N / 2L), each = 2L))
bfx <- matrix(exp(rnorm(G * 3L, 0, 0.3)), G, 3L)
gfx <- ifelse(runif(G) < 0.1, exp(rnorm(G, 0, 0.7)), 1)
mu <- mu0 * bfx[, as.integer(batch)] * outer(gfx, as.integer(group) - 1L, `^`)
counts <- matrix(rnbinom(G * N, mu = mu, size = 1 / 0.15), G, N,
                 dimnames = list(paste0("g", seq_len(G)), paste0("s", seq_len(N))))
storage.mode(counts) <- "integer"
design <- model.matrix(~ group + batch)

rows <- list()
tm <- function(pipeline, step, expr, reps = 3L) {
  e <- substitute(expr)
  env <- parent.frame()
  val <- NULL
  t1 <- system.time(val <- eval(e, env))[["elapsed"]]
  ts <- if (t1 < 5) c(t1, vapply(seq_len(reps - 1L), function(i) system.time(eval(e, env))[["elapsed"]], numeric(1))) else t1
  rows[[length(rows) + 1L]] <<- data.frame(pipeline = pipeline, step = step, seconds = min(ts), runs = length(ts))
  cat(sprintf("%-8s %-38s %9.2f s\n", pipeline, step, min(ts)))
  invisible(val)
}

cat(sprintf("%d genes x %d samples, 3 batches, 2 groups, %d subjects\n", G, N, nlevels(subj)))

tm("batch", "sva::ComBat_seq", suppressMessages(ComBat_seq(counts, batch = batch, group = group)))
lcpm <- cpm(counts, log = TRUE)
tm("batch", "sva::ComBat on log-CPM", suppressMessages(ComBat(lcpm, batch = batch)))

dge <- DGEList(counts)
keep <- tm("limma", "edgeR::filterByExpr", filterByExpr(dge, design))
dge <- dge[keep, , keep.lib.sizes = FALSE]
cat(sprintf("%d genes pass filterByExpr\n", nrow(dge)))
dge <- tm("limma", "edgeR::normLibSizes (TMM)", normLibSizes(dge))
v <- tm("limma", "limma::voom", voom(dge, design))
tm("limma", "limma::voomWithQualityWeights", voomWithQualityWeights(dge, design))
fit <- tm("limma", "limma::lmFit, voom weights", lmFit(v, design))
tm("limma", "limma::contrasts.fit", contrasts.fit(fit, coefficients = 2))
eb <- tm("limma", "limma::eBayes", eBayes(fit))
tm("limma", "limma::topTable", topTable(eb, coef = 2, number = Inf))
tm("limma", "limma::decideTests", decideTests(eb))
dc <- tm("limma", "limma::duplicateCorrelation", duplicateCorrelation(v, design, block = subj))
tm("limma", "limma::lmFit, block + correlation", lmFit(v, design, block = subj, correlation = dc$consensus.correlation))
tm("limma", "limma::removeBatchEffect", removeBatchEffect(v$E, batch = batch, design = model.matrix(~ group)))

dd <- tm("edgeR", "edgeR::estimateDisp", estimateDisp(dge, design))
qf <- tm("edgeR", "edgeR::glmQLFit", glmQLFit(dd, design))
qt <- tm("edgeR", "edgeR::glmQLFTest", glmQLFTest(qf, coef = 2))
tm("edgeR", "edgeR::topTags", topTags(qt, n = Inf))

dds <- DESeqDataSetFromMatrix(counts[keep, ], data.frame(group, batch), ~ batch + group)
dds <- tm("DESeq2", "DESeq2::DESeq, serial", suppressMessages(DESeq(dds)))
tm("DESeq2", "DESeq2::results", results(dds, name = "group_trt_vs_ctl"))
tm("DESeq2", "DESeq2::vst", suppressMessages(vst(dds)))

cat("\n")
print(do.call(rbind, rows), row.names = FALSE)
blas <- sessionInfo()$BLAS
cat("\n", R.version.string, " (", R.version$platform, ")  BLAS ",
    if (is.null(blas) || !nzchar(blas)) "R internal" else basename(blas),
    "  limma ", as.character(packageVersion("limma")),
    "  edgeR ", as.character(packageVersion("edgeR")), "  DESeq2 ", as.character(packageVersion("DESeq2")),
    "  sva ", as.character(packageVersion("sva")), "\n", sep = "")
