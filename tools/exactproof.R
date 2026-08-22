suppressMessages({library(sva); library(edgeR); library(limma); library(statmod)})
suppressMessages(pkgload::load_all(".", quiet = TRUE))
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1)

# Without this the harness runs at production defaults, every fixture here falls under a size
# gate, and the checks pass through the serial fallback while reporting the parallel path green.
# Same list as tests/testthat/setup-parallel.R, and it has to move with it.
options(combat.min.cells = 0, combat.min.disp.cells = 0, combat.min.ls.cells = 0,
        combat.min.norm.cells = 0, combat.min.order.cells = 0,
        combat.min.dupcor.cells = 0, combat.min.glm.cells = 0)
pass <- 0L; fail <- 0L
chk <- function(lbl, ok) {
  if (isTRUE(ok)) { pass <<- pass + 1L } else { fail <<- fail + 1L; cat(sprintf("  *** FAIL: %s\n", lbl)) }
}
q <- function(e) suppressMessages(e)

cat("=== ComBat_seq_parallel against sva::ComBat_seq ===\n")
mk <- function(G, nb, per, seed, mu = 200, size = 5, outlier = 0) {
  set.seed(seed); S <- nb * per
  b <- factor(rep(seq_len(nb), each = per))
  m <- matrix(rnbinom(G*S, mu = mu, size = size), G, S,
              dimnames = list(paste0("g", 1:G), paste0("s", 1:S)))
  if (outlier > 0) m[, 1] <- as.integer(m[, 1] * outlier)
  list(counts = m, batch = b)
}
cfgs <- list(
  list(lbl="ordinary 800x24 2 batches",      d=mk(800L,2L,12L,1), grp=NULL,  cov=NULL),
  list(lbl="ordinary 600x30 5 batches",      d=mk(600L,5L,6L,2),  grp=NULL,  cov=NULL),
  list(lbl="low counts, non-convergent",     d=mk(400L,2L,10L,3, mu=2, size=0.3), grp=NULL, cov=NULL),
  list(lbl="over-sequenced library x1000",   d=mk(240L,2L,6L,2024, mu=2, size=0.3, outlier=1000), grp=NULL, cov=NULL),
  list(lbl="over-sequenced x60, 4 batches",  d=mk(500L,4L,8L,7, mu=3, size=0.4, outlier=60), grp=NULL, cov=NULL)
)
for (cf in cfgs) {
  d <- cf$d
  ref <- q(sva::ComBat_seq(d$counts, batch = d$batch, group = cf$grp))
  for (w in c(2L,4L,8L)) for (k in list(NULL, 2L, 3L, 7L)) {
    p <- q(ComBat_seq_parallel(d$counts, batch = d$batch, group = cf$grp, workers = w, chunks = k))
    chk(sprintf("%s w=%d k=%s", cf$lbl, w, if (is.null(k)) "default" else k), identical(ref, p))
  }
}
# with a group, and with covar_mod
d <- mk(600L, 3L, 12L, 11)
g <- factor(rep(rep(c("a","b"), each = 6), 3))
ref <- q(sva::ComBat_seq(d$counts, batch = d$batch, group = g))
for (w in c(2L,4L,8L)) chk(sprintf("group supplied w=%d", w),
  identical(ref, q(ComBat_seq_parallel(d$counts, batch = d$batch, group = g, workers = w))))
cv <- model.matrix(~ factor(rep(rep(c("x","y","z"), length.out = 12), 3)))
ref <- q(sva::ComBat_seq(d$counts, batch = d$batch, group = NULL, covar_mod = cv))
for (w in c(2L,4L,8L)) chk(sprintf("covar_mod w=%d", w),
  identical(ref, q(ComBat_seq_parallel(d$counts, batch = d$batch, group = NULL, covar_mod = cv, workers = w))))
# every backend
for (bk in c("mclapply","serial")) chk(sprintf("backend %s", bk),
  identical(ref, q(ComBat_seq_parallel(d$counts, batch = d$batch, group = NULL, covar_mod = cv,
                                       workers = 4L, parallel_backend = bk))))
cat(sprintf("  ComBat: %d checks\n", pass + fail))

cat("\n=== limma and edgeR companions ===\n")
set.seed(21); G <- 900L; S <- 24L
grp <- factor(rep(c("A","B"), each = S/2)); blk <- factor(rep(1:12, each = 2))
y <- matrix(rnbinom(G*S, mu = 200, size = 5), G, S,
            dimnames = list(paste0("g",1:G), paste0("s",1:S)))
y[, 1] <- as.integer(y[, 1] * 1000)
des <- model.matrix(~ grp)
for (m in c("TMM","TMMwsp","RLE","upperquartile","none"))
  chk(paste("calcNormFactors", m),
      identical(calcNormFactors(y, method = m), calcNormFactors_parallel(y, method = m, workers = 4L)))
chk("calcNormFactors DGEList", identical(calcNormFactors(DGEList(y)), calcNormFactors_parallel(DGEList(y), workers = 4L)))
v <- voom(calcNormFactors(DGEList(y)), des)
for (k in c(1L,2L,4L,7L)) chk(paste("lmFit chunks", k),
  identical(lmFit(v, des), lmFit_parallel(v, des, workers = 4L, chunks = k)))
aw <- arrayWeights(v$E, des)
chk("lmFit array weights", identical(lmFit(v$E, des, weights = aw), lmFit_parallel(v$E, des, weights = aw, workers = 4L)))
dcr <- duplicateCorrelation(v, des, block = blk)
for (k in c(1L,2L,4L,7L)) chk(paste("dupcor chunks", k),
  identical(dcr, duplicateCorrelation_parallel(v, des, block = blk, workers = 4L, chunks = k)))
chk("lmFit blocked", identical(lmFit(v, des, block = blk, correlation = dcr$consensus.correlation),
     lmFit_parallel(v, des, block = blk, correlation = dcr$consensus.correlation, workers = 4L)))
cm <- suppressWarnings(makeContrasts(grpB, levels = des))
ref_tt <- topTable(eBayes(contrasts.fit(lmFit(voom(calcNormFactors(DGEList(y)), des), des), cm)), number = Inf, sort.by = "none")
par_tt <- topTable(eBayes(contrasts.fit(lmFit_parallel(voom(calcNormFactors_parallel(DGEList(y), workers=4L), des), des, workers=4L), cm)), number = Inf, sort.by = "none")
chk("END TO END final gene list", identical(ref_tt, par_tt))

cat(sprintf("\n==== %d checks, %d passed, %d failed ====\n", pass + fail, pass, fail))
if (fail > 0) quit(status = 1)
