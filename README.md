# RNA-Parallel

Parallel companions for RNA-seq tools. Each calls the original function unmodified and returns
output `identical()` to it, bit for bit. Same arguments, same defaults, same result, faster.

Rendered analysis: [macOS](https://genomerx.github.io/RNA-Parallel/) ·
[Linux](https://genomerx.github.io/RNA-Parallel/linux.html)

## Speedup

TCGA, 18,270 genes by 1,500 tumours. Every arm `identical()` to the vendor on both platforms.

| companion | runs | macOS<br>M3, 8 cores | Linux<br>2x Xeon, 16 cores |
|---|---|---:|---:|
| `ComBat_seq_parallel()` | `sva::ComBat_seq` | 5.37x @ 8w | **9.38x @ 16w** |
| `calcNormFactors_parallel()` | `edgeR::normLibSizes` | **6.86x @ 8w** | 4.37x @ 8w |
| `lmFit_parallel()` | `limma::lmFit` | 2.66x @ 4w | **3.42x @ 16w** |
| `duplicateCorrelation_parallel()` | `limma::duplicateCorrelation` | 3.63x @ 8w | **7.31x @ 16w** |
| `removeBatchEffect_parallel()` | `limma::removeBatchEffect` | **3.37x @ 8w** | 1.40x @ 8w |

Originals are timed on both sides of the companion arms and averaged; short stages are the best
of three. Bold marks the faster platform. Machine specs and what separates them are in
[Cross-platform](#cross-platform).

**Nothing is reimplemented.** The vendor function is the one that runs. It is called with its
hot paths rebound in a child of its own environment, so every other symbol still resolves to
vendor code. A companion is worth keeping only if it returns exactly what it replaces, so the
suite asserts `identical()` rather than a tolerance.

## Install

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("sva", "edgeR", "limma"))

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("GenomeRx/RNA-Parallel")
```

## Use

Each companion takes its vendor's arguments in the same order with the same defaults, and adds
`workers`, `chunks`, `parallel_backend` and `backend`. All five in one pass:

```r
library(rnaparallel); library(edgeR); library(limma)

# raw counts, batch corrected
adjusted <- ComBat_seq_parallel(counts, batch = batch, group = NULL, workers = 8L)

# limma-voom differential expression
dge <- calcNormFactors_parallel(DGEList(adjusted), workers = 8L)
v   <- voom(dge, design)
fit <- lmFit_parallel(v, design, workers = 8L)
tt  <- topTable(eBayes(fit), coef = 2, number = Inf)

# blocked design, repeated measures on one subject
cor <- duplicateCorrelation_parallel(v, design, block = subject, workers = 8L)
fit <- lmFit_parallel(v, design, block = subject,
                      correlation = cor$consensus.correlation, workers = 8L)

# batch out of a log-expression matrix, for PCA and heatmaps
vis <- removeBatchEffect_parallel(v$E, batch = batch, design = design, workers = 8L)
```

One interface difference: `duplicateCorrelation_parallel` requires `block`, whereas the vendor
defaults it to `NULL`. [How it works](#how-it-works) says why.

### Check it yourself

All five companions against their vendors. No download, nothing to configure:

```r
library(rnaparallel); library(sva); library(edgeR); library(limma)

set.seed(1)
counts  <- matrix(rnbinom(4000, mu = 50, size = 5), nrow = 500)
batch   <- rep(1:2, each = 4)
group   <- rep(0:1, 4)
design  <- model.matrix(~ group)
subject <- rep(1:4, each = 2)

dge <- DGEList(counts)
v   <- voom(normLibSizes(dge), design)

c(ComBat_seq        = identical(ComBat_seq_parallel(counts, batch, group = NULL, workers = 4L),
                                ComBat_seq(counts, batch, group = NULL)),
  normLibSizes      = identical(calcNormFactors_parallel(dge, workers = 4L),
                                normLibSizes(dge)),
  lmFit             = identical(lmFit_parallel(v, design, workers = 4L),
                                lmFit(v, design)),
  duplicateCor      = identical(duplicateCorrelation_parallel(v, design, ndups = 1,
                                                              block = subject, workers = 4L),
                                duplicateCorrelation(v, design, ndups = 1, block = subject)),
  removeBatchEffect = identical(removeBatchEffect_parallel(v$E, batch, design = design,
                                                           workers = 4L),
                                removeBatchEffect(v$E, batch, design = design)))
#>        ComBat_seq      normLibSizes             lmFit      duplicateCor
#>              TRUE              TRUE              TRUE              TRUE
#> removeBatchEffect
#>              TRUE
```

Four thousand cells is under every size gate, so those arms run serially and are still
`identical()`. [run_example.R](inst/examples/run_example.R) runs the same comparison at a size
where the speedup shows, in minutes, still with nothing to download.

```sh
Rscript inst/examples/run_example.R
```

Rendered reports carry it at cohort scale, same sections on both platforms: cohort, argument
parity, corrections, PCA, differential expression, worker sweep.
[macOS](https://genomerx.github.io/RNA-Parallel/) ·
[Linux](https://genomerx.github.io/RNA-Parallel/linux.html). Sources are
[RNA_Parallel.Rmd](inst/examples/RNA_Parallel.Rmd) and
[RNA_Parallel_linux.Rmd](inst/examples/RNA_Parallel_linux.Rmd), the second rendered with
[render_linux.sh](inst/examples/render_linux.sh), which pins BLAS before R starts. The first run
downloads HNSC, LUAD and LUSC to the per-user cache, or `RNAPARALLEL_TCGA_DIR` if set.

[tests/](tests/testthat) covers every argument path, chunk layout and backend, plus a dispatch
count through the public entry point, because `identical()` alone cannot tell a working parallel
layer from a dead one. Over 400 assertions.

### Seeing what is running

```r
options(combat.timing = TRUE)      # one elapsed line per call
options(combat.quiet  = TRUE)      # swallow the vendor's progress chatter
options(combat.timing.min = 1)     # hide anything faster than a second
```

```
  Breast ComBat-seq                  mclapply x6         23.2s
  TMM 12,000 x 700                   mclapply x6          0.9s
  TMM 200 x 8                        serial               0.0s  2 gated
```

The engine column reports what **ran**, not what was asked for: a call that fell under a size
gate says `serial`, because a companion returning `identical()` output at serial pace otherwise
looks exactly like one that forked. If a pinned vendor excerpt stands down after an upstream
change, the line says `[match_quantiles stood down]`.

**`rnaparallel_stale()`** returns TRUE if the package was reinstalled under a running session.
The fix is to restart R.

## How it works

`ComBat_seq_parallel()` creates a child of the vendor's environment, binds six names in it, and
calls `sva::ComBat_seq` unchanged. Shares are of serial time on 10,000 genes by 500 samples.

| ComBat-seq calls | share | split by | why that is exact |
|---|---|---|---|
| `match_quantiles` | 66.5% | gene rows | each output cell reads only its own gene |
| `estimateGLMTagwiseDisp` | 14.0% | gene rows | valid at `prior.df = 0`; other values go to `edgeR` unsplit |
| `estimateGLMCommonDisp` | 13.1% | batches | it sums over all genes, so a row split changes accumulation order; batches are independent and RNG-free |
| `glmFit`, `glmFit.default` | remainder | gene rows | offset, dispersion, weights and start arrive explicitly and slice with the rows |
| `monte_carlo_int_NB` | small | serial | draws depend on the previous batch's state |

`calcNormFactors_parallel()` splits by sample column, `lmFit_parallel()` and
`duplicateCorrelation_parallel()` by gene row. `removeBatchEffect_parallel()` splits nothing: it
rebinds the one `lmFit` call inside the vendor to `lmFit_parallel()` and inherits its curve, which
still pays, because the vendor does not finish inside ten minutes at 9,493 samples. Chunks are
interleaved so a sorted matrix does not leave workers idle, each carries a tag, and a dead worker,
duplicate chunk or short result halts the run.

**What stays serial keeps the output identical:** anything reducing across genes (`eBayes`,
`squeezeVar`, `fitFDist`, `arrayWeights`, `normalizeBetweenArrays`, `normalizeQuantiles`) or
scaling with gene count (`voom`'s `lowess` span, `p.adjust` and so `topTable` and `decideTests`).
`contrasts.fit` takes 0.001 s, less than a fork. `lmFit(method = "robust")` and `ndups >= 2`
reshape the row axis through `unwrapdups` and error rather than split, which is why `block` is
required for `duplicateCorrelation_parallel`.

**Built, measured, deleted,** because a modest speedup does not buy different numbers.

| companion | reached | what it moved |
|---|---|---|
| `estimateDisp` | 1.4x | 19,999 of 20,000 tagwise dispersions, once one library was over-sequenced |
| `glmQLFit` | 1.6x | 22 of 18,270 TCGA genes, where `mglmLevenberg` records a non-convergent deviance that depends on which genes share the block, with no flag set |

**Three ways a limma row split goes wrong,** each reproduced against limma 3.62.2, each returning
a wrong answer quietly rather than failing.

| trap | wrong answer | guard |
|---|---|---|
| `asMatrixWeights` dispatches on the *block's* row count | a per-array weight vector reads as per-gene weights in any block whose height equals the sample count | weights expanded once against the whole matrix |
| `NoProbeWts` is an AND-reduction selecting between two numerically different algorithms | an all-finite block takes the fast path while the whole matrix took the slow one | the branch is proven unable to flip before anything is split |
| `stats::lm.fit` drops a one-column response to a vector | a one-gene block swaps `colMeans` for `mean` | at least two genes per chunk |

One inner loop is pinned to the vendor body it came from and stands down on any upstream change:
`sva::match_quantiles`, 66.5% of ComBat-seq's serial time. See [License](#license).

## Cross-platform

The same package on two unlike machines. Every arm on both returns `identical()` output: macOS
and Linux, fork and socket dispatch, threaded and single-threaded BLAS.

| | macOS | Linux |
|---|---|---|
| CPU | Apple M3 | 2 x Intel Xeon @ 2.80 GHz |
| Physical cores | 8 (4 performance + 4 efficiency) | 16 |
| Logical | 8 | 32 |
| NUMA nodes | 1 | 2 |
| BLAS | Accelerate | OpenBLAS 0.3.26 |
| Workers swept | 2, 4, 8 | 4, 8, 16, 32 |

Absolute wall clock at cohort scale, every companion. Vendor is the serial original; best is the
fastest companion arm on that machine, with its worker count. Bold marks the faster platform.

| stage | macOS vendor | macOS best | Linux vendor | Linux best |
|---|---:|---:|---:|---:|
| ComBat-seq | 1,653.9 s | **308.0 s** (8w) | 3,108.8 s | 331.5 s (16w) |
| duplicateCorrelation | 660.5 s | 182.0 s (8w) | 519.3 s | **71.0 s** (16w) |
| calcNormFactors (TMM) | 10.3 s | **1.5 s** (8w) | 21.4 s | 4.9 s (8w) |
| lmFit | 4.4 s | **1.7 s** (4w) | 11.9 s | 3.5 s (16w) |
| removeBatchEffect | 4.5 s | **1.3 s** (8w) | 2.6 s | 1.9 s (8w) |

Speedup is a ratio and rewards a slower core: Linux scales further, 9.38x against 5.37x, and
still finishes behind at 331.5 s against 308.0 s. Read the seconds. On Linux, pin
`OPENBLAS_NUM_THREADS` before R starts or every forked worker opens its own thread pool, though it
changes little here: unpinning moves DGEMM throughput 5x, from 62 to 314 GFLOPS, and `lmFit` by
0.3%, because limma solves a small QR per gene.

## Tuning

| knob | default | change it when |
|---|---|---|
| `workers` | `min(8, detectCores() - 2)` | rarely. Going past your *physical* core count can be slower, as the 32-worker Linux arm shows. See `?workers` |
| `chunks` | `workers` | only to cut peak memory per worker |
| `parallel_backend` | `"mclapply"` | you cannot fork, or a cluster is already running |

`workers` is not a safety ceiling. Six of them alongside a second forking R session have
kernel-panicked a 24 GB machine, and nothing here can see that session.

**Backends** are `"mclapply"` (forks), `"future"`, `"BiocParallel"`, `"foreach"`, `"serial"`, or any
`function(idx, f, workers)`, all returning identical results. Forking is the default because a
forked worker reads the matrix copy-on-write, while socket backends re-serialise it per chunk and
measured slower than not parallelising at all. **Windows** cannot fork: `mclapply` runs serially and
says so once, `BiocParallel` substitutes a serial parameter, `"foreach"` runs over PSOCK and is
slower than serial here. **Nesting** is safe on the default backend, where `mc.allow.recursive` is
`FALSE`; `foreach` builds its own cluster and does multiply, so spend the worker budget *inside* a
loop body rather than across it. Inverting one 15-cohort screen measured 245 s to 81.6 s.

**Size gates** are why a small input shows no speedup. Seven options set the cell count below which
a call runs serially: `combat.min.cells` (20,000), `combat.min.disp.cells` (30,000),
`combat.min.glm.cells` (100,000), `combat.min.ls.cells` (6e6), `combat.min.norm.cells` (2e5),
`combat.min.order.cells` (4e6), `combat.min.dupcor.cells` (5,000).
`options(combat.fork = FALSE)` forces serial everywhere, `combat_cluster_stop()` releases cached
clusters.

**`calcNormFactors_parallel()`** wraps `normLibSizes` on current edgeR and `calcNormFactors` on
older ones. `normLibSizes` **errors** on negative counts where the old name returned NaN-warned
factors.

## License

MIT for this companion, copyright GenomeRx 2026, in [LICENSE](LICENSE). The vendor packages are
called at run time from your own installation and none of them is redistributed here: sva is
Artistic-2.0, limma and edgeR are GPL (>= 2).

One block is the exception and is marked as such in the source. `R/helper_seq_parallel.R` carries
the deparsed `sva::match_quantiles` body and a row-vectorised transcription of it, so the
companion can detect an upstream change and stand down. It derives from Artistic-2.0 code by
Zhang, Parmigiani and Johnson. Nothing else here reproduces vendor code.

## Citation

Cite the method you used and this companion. This repository adds parallelism and contributes no
statistics, so the method citation is the one carrying the result. `citation("rnaparallel")`
prints every entry below.

**ComBat-seq.** Zhang Y, Parmigiani G, Johnson WE (2020). ComBat-seq: batch effect adjustment for
RNA-seq count data. *NAR Genomics and Bioinformatics* 2(3), lqaa078.
doi:[10.1093/nargab/lqaa078](https://doi.org/10.1093/nargab/lqaa078).
Package: <https://bioconductor.org/packages/release/bioc/html/sva.html>.

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
