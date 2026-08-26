## Compare sva::ComBat_seq with ComBat_seq_parallel: identical output, measured speedup.
##
## Simulated counts, so it runs anywhere in a few minutes with no download and no
## configuration. Counts follow the generative model from the ComBat-seq paper's simulations.
## The rendered TCGA report covers the same comparison at cohort scale:
## https://genomerx.github.io/RNA-Parallel/
##
##   Rscript run_example.R [genes] [samples per batch-condition cell] [workers]
##
## Nothing is written to disk.

args    <- commandArgs(trailingOnly = TRUE)
G       <- as.integer(if (length(args) >= 1) args[1] else 20000L)
n_cell  <- as.integer(if (length(args) >= 2) args[2] else 120L)   # samples per batch-condition cell
workers <- as.integer(if (length(args) >= 3) args[3] else 4L)

need <- c(rnaparallel = 'remotes::install_github("GenomeRx/RNA-Parallel")',
          sva = 'BiocManager::install("sva")',
          edgeR = 'BiocManager::install("edgeR")',
          limma = 'BiocManager::install("limma")')
missing <- names(need)[!vapply(names(need), requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("install first:\n  ", paste(need[missing], collapse = "\n  "), call. = FALSE)
}
suppressMessages({library(sva); library(edgeR); library(limma); library(rnaparallel)})

set.seed(123)
bio_fold <- 1.5; batch_fold <- 3; disp_fold <- 4; disp_1 <- 0.15

batch <- factor(rep(c("1", "2"), each = 2 * n_cell))
group <- factor(rep(rep(c("control", "treated"), each = n_cell), 2))
lv    <- c("1.control", "2.control", "1.treated", "2.treated")
cell  <- factor(paste(batch, group, sep = "."), levels = lv)

base     <- round(runif(G, 20, 600))
de_idx   <- sample.int(G, max(2L, G %/% 9L))
ups      <- de_idx[seq_len(length(de_idx) %/% 2)]
downs    <- setdiff(de_idx, ups)
batch_up <- sample.int(G, G %/% 2)

fc <- matrix(1, G, 4, dimnames = list(NULL, lv))
fc[batch_up,  c(2, 4)] <- batch_fold        # batch 2 raises one half of the genes
fc[-batch_up, c(2, 4)] <- 1 / batch_fold    # and lowers the other
fc[ups,   c(3, 4)] <- fc[ups,   c(3, 4)] * bio_fold
fc[downs, c(1, 2)] <- fc[downs, c(1, 2)] * bio_fold

size <- matrix(1 / disp_1, G, 4); size[, c(2, 4)] <- 1 / (disp_1 * disp_fold)

k <- as.integer(cell)
counts <- vapply(seq_along(k), function(j) {
  x <- rnbinom(G, mu = base * fc[, k[j]], size = size[, k[j]])
  x[x == 0] <- 1L                            # polyester replaces zero draws with 1
  x
}, numeric(G))
storage.mode(counts) <- "integer"
rownames(counts) <- sprintf("gene%06d", seq_len(G))

cat(sprintf("\n%s genes x %d samples, 2 batches, %s truly differential\n",
            format(G, big.mark = ","), ncol(counts), format(length(de_idx), big.mark = ",")))
cat(sprintf("mean batch effect %gx, dispersion ratio %gx, biological signal %gx\n\n",
            batch_fold, disp_fold, bio_fold))

t0  <- Sys.time(); ref <- sva::ComBat_seq(counts, batch = batch, group = group); t1 <- Sys.time()
par <- ComBat_seq_parallel(counts, batch = batch, group = group, workers = workers); t2 <- Sys.time()
s_ref <- as.numeric(difftime(t1, t0, units = "secs"))
s_par <- as.numeric(difftime(t2, t1, units = "secs"))

cat(sprintf("sva::ComBat_seq                  %7.1f s\n", s_ref))
cat(sprintf("ComBat_seq_parallel (%d workers)  %7.1f s   %.2fx\n\n", workers, s_par, s_ref / s_par))

pca <- function(m) {
  p <- prcomp(t(cpm(m, log = TRUE, prior.count = 1)), scale. = FALSE)
  list(x = p$x, ve = p$sdev^2 / sum(p$sdev^2))
}
p_ref <- pca(ref); p_par <- pca(par)

des <- model.matrix(~ group)
norm_edgeR <- if (exists("normLibSizes.default", envir = asNamespace("edgeR"), inherits = FALSE)) {
  edgeR::normLibSizes
} else {
  edgeR::calcNormFactors
}
de <- function(m) topTable(eBayes(lmFit(voom(norm_edgeR(DGEList(m)), des), des)),
                            coef = 2, number = Inf, sort.by = "none")
t_ref <- de(ref); t_par <- de(par)

checks <- c(
  "corrected counts"        = identical(ref, par),
  "principal components"    = identical(p_ref$x, p_par$x),
  "variance explained"      = identical(p_ref$ve, p_par$ve),
  "log fold changes"        = identical(t_ref$logFC, t_par$logFC),
  "adjusted p-values"       = identical(t_ref$adj.P.Val, t_par$adj.P.Val),
  "significant gene set"    = identical(rownames(t_ref)[t_ref$adj.P.Val < 0.05],
                                        rownames(t_par)[t_par$adj.P.Val < 0.05]))
for (nm in names(checks)) cat(sprintf("  %-24s identical(): %s\n", nm, checks[[nm]]))

sig   <- t_par$adj.P.Val < 0.05
is_de <- seq_len(G) %in% de_idx
cat(sprintf("\nrecovery against known truth: %s of %s planted genes, %s false positives\n",
            format(sum(sig & is_de), big.mark = ","), format(length(de_idx), big.mark = ","),
            format(sum(sig & !is_de), big.mark = ",")))
cat(sprintf("batch R2 on PC1: %.4f uncorrected, %.4f corrected (chance floor %.4f)\n",
            summary(lm(pca(counts)$x[, 1] ~ batch))$r.squared,
            summary(lm(p_ref$x[, 1] ~ batch))$r.squared, 1 / (ncol(counts) - 1)))

if (!all(unlist(checks))) stop("implementations diverged", call. = FALSE)
cat("\nall stages identical\n")
