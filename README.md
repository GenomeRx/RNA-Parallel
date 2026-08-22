# RNA-Parallel

Companion plugins for RNA-seq tools. Each one calls the original function unmodified and returns
output `identical()` to it, bit for bit. Same arguments, same defaults, same result, faster.

**[Read the rendered analysis](https://genomerx.github.io/RNA-Parallel/)**

| companion | runs | axis | 2 workers | 4 workers | 8 workers |
|---|---|---|---:|---:|---:|
| `ComBat_seq_parallel()` | `sva::ComBat_seq` | gene rows, batches | 1.85x | 3.50x | **5.45x** |
| `calcNormFactors_parallel()` | `edgeR::calcNormFactors` | sample columns | 3.74x | 6.53x | **8.79x** |
| `lmFit_parallel()` | `limma::lmFit` | gene rows | 2.64x | 4.41x | **4.44x** |
| `duplicateCorrelation_parallel()` | `limma::duplicateCorrelation` | gene rows | 1.96x | 3.18x | **4.45x** |

TCGA, 18,270 genes by 1,500 tumours, Apple M3, every arm `identical()` to the vendor. Each
original is timed on both sides of the companion arms and averaged, so drift in machine load
moves both halves of a ratio together; the two readings differed by at most 3.2% here. `lmFit`
gains almost nothing between four workers and eight. `calcNormFactors` clears the worker count
because part of its gain is a serial fix to TMM's double `rank()` rather than parallelism. Full
method, cohort and figures in the [report](https://genomerx.github.io/RNA-Parallel/).

**Nothing here is reimplemented.** The vendor function is the one that executes. It is called
with its hot paths rebound in a child of its own environment, so every other symbol still
resolves to the vendor's code. A companion earns its place only by returning exactly what it
replaces, so the suite asserts `identical()` rather than a tolerance.

One inner loop is the exception, pinned to the vendor body it came from and standing down on any
upstream change: `sva::match_quantiles`, 66.5% of ComBat-seq's serial time. See
[License](#7-license). Everywhere else, including the design-invariant work inside
`statmod::mixedModel2Fit`, the vendor's own object runs and only the primitives it calls are
memoised.

## 1. Install

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("sva", "edgeR", "limma"))

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("GenomeRx/RNA-Parallel")
```

## 2. Use

Each companion takes its vendor's arguments in the same order with the same defaults, and adds
`workers`, `chunks`, `parallel_backend` and `backend`. Swapping one in changes runtime and
nothing else.

```r
library(rnaparallel); library(edgeR); library(limma)

adjusted <- ComBat_seq_parallel(counts, batch = batch, group = NULL, workers = 4L)

# limma-voom differential expression
dge <- calcNormFactors_parallel(DGEList(counts), workers = 8L)
v   <- voom(dge, design)
fit <- lmFit_parallel(v, design, workers = 8L)
tt  <- topTable(eBayes(fit), coef = 2, number = Inf)

# blocked design
cor <- duplicateCorrelation_parallel(v, design, block = subject, workers = 8L)
fit <- lmFit_parallel(v, design, block = subject,
                      correlation = cor$consensus.correlation, workers = 8L)
```

The one interface difference: `duplicateCorrelation_parallel` requires `block`, where the vendor
defaults it to `NULL`. Section 3 says why.

Check the claim yourself in ten lines:

```r
set.seed(1)
counts <- matrix(rnbinom(1600, mu = 50, size = 5), nrow = 200)
batch  <- rep(1:2, each = 4)

identical(ComBat_seq_parallel(counts, batch, group = NULL, workers = 4L),
          sva::ComBat_seq(counts, batch, group = NULL))
#> TRUE
```

## 3. How it works

`ComBat_seq_parallel()` creates a child of the vendor's own environment, binds six names in it,
and calls `sva::ComBat_seq` unchanged. Shares are of serial time on 10,000 genes by 500 samples.

| ComBat-seq calls | share | split by | why that is exact |
|---|---|---|---|
| `match_quantiles` | 66.5% | gene rows | each output cell reads only its own gene |
| `estimateGLMTagwiseDisp` | 14.0% | gene rows | valid at `prior.df = 0`; other values go to `edgeR` unsplit |
| `estimateGLMCommonDisp` | 13.1% | batches | it sums over all genes, so a row split changes accumulation order; batches are independent and RNG-free |
| `glmFit`, `glmFit.default` | remainder | gene rows | offset, dispersion, weights and start arrive explicitly and slice with the rows |
| `monte_carlo_int_NB` | small | serial | draws depend on the previous batch's state |

Assembly is what makes a split safe rather than merely fast. Chunks are interleaved so a sorted
matrix does not leave workers idle, each carries a tag and is reassembled by it, and a dead
worker, duplicate chunk or short result halts the run instead of being patched.

**What is deliberately not parallelised.** This list is the reason the output is identical, not
an admission. `voom`'s `lowess` trend takes its span as a *fraction* of gene count, so a block
fits a different curve (measured 0.99x). `eBayes`, `squeezeVar` and `fitFDist` pool across all
genes. `p.adjust` is length-dependent, so `topTable` and `decideTests` cannot be blocked.
`arrayWeights`, `normalizeBetweenArrays` and `normalizeQuantiles` reduce across genes by
construction. `contrasts.fit` is exactly splittable and takes 0.001 s, so a fork costs more than
the work. `lmFit(method = "robust")` and `ndups >= 2` reshape the row axis through `unwrapdups`,
so a split would pair different genes; both are refused with an error rather than left silent,
which is also why `block` is required for `duplicateCorrelation_parallel`.

Two edgeR companions were built, measured and then deleted. `estimateDisp` reached 1.4x but moved
19,999 of 20,000 tagwise dispersions once one library was heavily over-sequenced. `glmQLFit`
reached 1.6x, passed every small fixture, then failed on real TCGA data at 22 of 18,270 genes:
`mglmLevenberg` records a deviance whose value depends on which genes share the block when a fit
does not converge, with no flag set. A modest speedup does not buy a companion that returns
different numbers. The same kernel is why a batch-only ComBat-seq design is dispatched across
batches rather than gene rows.

**Three ways a limma row split goes wrong,** each reproduced against limma 3.62.2 before being
guarded, and each returning a wrong answer quietly rather than failing. `asMatrixWeights`
dispatches on the *block's* row count, so a per-array weight vector is read as per-gene weights
by any block whose height equals the sample count (coefficients moved by 0.325). `NoProbeWts` is
an AND-reduction over every cell that selects between two numerically different algorithms, so an
all-finite block flips to the fast path while the whole matrix took the slow one (114 of 400
`sigma` differed). `stats::lm.fit` drops a one-column response to a vector, so a one-gene block
swaps `colMeans` for `mean` (4,080 of 16,000 one-gene blocks were not `identical()`). Every chunk
now carries at least two genes, weights are expanded once against the whole matrix, and the
branch is proven unable to flip before anything is split.

## 4. Verify it yourself

- [**Rendered report**](https://genomerx.github.io/RNA-Parallel/). Cohort, argument parity, both
  corrections, PCA before and after, differential expression, worker sweep.
- **[run_example.R](inst/examples/run_example.R)**. The whole claim in a few minutes, no download.

  ```sh
  Rscript inst/examples/run_example.R
  ```
- **[RNA_Parallel.Rmd](inst/examples/RNA_Parallel.Rmd)**. Source of the report. Open it and knit;
  the first run downloads the HNSC, LUAD and LUSC cohorts to the per-user cache, or to
  `RNAPARALLEL_TCGA_DIR` if set.
- **[tests/](tests/testthat)**. Every argument path, chunk layout and backend, plus a dispatch
  count through the public entry point, because `identical()` alone cannot tell a working
  parallel layer from a dead one. Over 350 assertions; `R CMD check` runs fewer, skipping the
  three that need a live cluster.

## 5. Tuning

`workers` is the primary parameter and defaults to four everywhere. Six workers alongside a second
forking R session has caused a kernel panic on a 24 GB machine. Added workers also return less
throughput each, so more is not always faster: `lmFit` above gains almost nothing between four and
eight.

`parallel_backend` selects the framework: `"mclapply"` (default, forks), `"future"`,
`"BiocParallel"`, `"foreach"`, `"serial"`, or any `function(idx, f, workers)`. All return
identical results. Forking is the default because a forked worker reads the count matrix through
copy-on-write instead of being sent a copy.

On Windows nothing forks. `mclapply` runs serially and says so once per session. `MulticoreParam`
does not exist there, so `BiocParallel` substitutes a serial parameter of its own and this package
says nothing. `"foreach"` does run in parallel on Windows, over a PSOCK cluster, but choose it for
correctness rather than speed: `doParallel` serialises each task's closure and the closure
captures the matrix, which measured 0.24x to 0.41x against `limma::lmFit` at four workers, slower
than not parallelising at all.

Seven `combat.min.*` options set the size below which a call runs serially, because under them a
fork costs more than the work it saves. They are the reason a small input shows no speedup.
`combat.min.dupcor.genes` became `combat.min.dupcor.cells` in 0.4.4, counted in gene-by-array
cells with a default of 5000; the old name is no longer read.

`options(combat.fork = FALSE)` forces serial on every backend including a custom executor, and
refuses a non-logical value rather than silently disabling itself. `combat_cluster_stop()`
releases cached clusters.

## 6. Citation

Cite the method you used **and** this companion. This repository adds parallelism and contributes
no statistics. `citation("rnaparallel")` prints all of it.

**ComBat-seq.** Zhang Y, Parmigiani G, Johnson WE (2020). ComBat-seq: batch effect adjustment for
RNA-seq count data. *NAR Genomics and Bioinformatics* 2(3), lqaa078.
doi:[10.1093/nargab/lqaa078](https://doi.org/10.1093/nargab/lqaa078). Distributed in
[sva](https://bioconductor.org/packages/release/bioc/html/sva.html).

**limma.** Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). limma powers
differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids
Research* 43(7), e47. doi:[10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007).
Package: <https://bioconductor.org/packages/release/bioc/html/limma.html>.

**edgeR.** Robinson MD, McCarthy DJ, Smyth GK (2010). edgeR: a Bioconductor package for
differential expression analysis of digital gene expression data. *Bioinformatics* 26(1),
139-140. doi:[10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616).
Package: <https://bioconductor.org/packages/release/bioc/html/edgeR.html>.

**This companion.** Nguyen N (2026). *rnaparallel: Parallel Companions for RNA-Seq Tools.*
<https://github.com/GenomeRx/RNA-Parallel>.

## 7. License

MIT for this companion, copyright GenomeRx. The vendor packages are executed, not forked.

One file reproduces vendor source, so the companion can detect an upstream change and stand
down: `R/helper_seq_parallel.R` carries the deparsed `sva::match_quantiles` body and a
row-vectorised transcription of it, derived from Artistic-2.0 code by Zhang, Parmigiani and
Johnson. The block is marked inline. Nothing else here reproduces vendor code.
