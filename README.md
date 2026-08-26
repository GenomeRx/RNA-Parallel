# RNA-Parallel

Parallel companions for RNA-seq tools. Each calls the original function unmodified and returns
output `identical()` to it, bit for bit. Same arguments, same defaults, same result, faster.

Rendered analysis: [macOS](https://genomerx.github.io/RNA-Parallel/) ·
[Linux](https://genomerx.github.io/RNA-Parallel/linux.html)

## Speedup

TCGA, 18,270 genes by 1,500 tumours. Every arm `identical()` to the vendor on both platforms.

| companion | runs | macOS<br>M3, 8 cores | Linux<br>2x Xeon, 16 cores |
|---|---|---:|---:|
| `ComBat_seq_parallel()` | `sva::ComBat_seq` | 5.37x @ 8w | **7.31x @ 16w** |
| `calcNormFactors_parallel()` | `edgeR::normLibSizes` | **6.86x @ 8w** | 3.89x @ 8w |
| `lmFit_parallel()` | `limma::lmFit` | 2.66x @ 4w | **2.80x @ 16w** |
| `duplicateCorrelation_parallel()` | `limma::duplicateCorrelation` | 3.63x @ 8w | **6.00x @ 16w** |
| `removeBatchEffect_parallel()` | `limma::removeBatchEffect` | **3.37x @ 8w** | 1.22x @ 4w |

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
`workers`, `chunks`, `parallel_backend` and `backend`.

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

One interface difference: `duplicateCorrelation_parallel` requires `block`, whereas the vendor
defaults it to `NULL`. [How it works](#how-it-works) says why.

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

The limma and edgeR companions each split one axis: `calcNormFactors_parallel()` by sample
column, `lmFit_parallel()` and `duplicateCorrelation_parallel()` by gene row.

Assembly is what makes a split safe rather than only fast. Chunks are interleaved so a sorted
matrix does not leave workers idle, each carries a tag and is reassembled by it, and a dead
worker, duplicate chunk or short result halts the run instead of being patched.

**One companion is second-order.** `removeBatchEffect_parallel()` splits nothing itself.
`limma::removeBatchEffect` is argument shaping around a single `lmFit` call with no cross-gene
reduction, so the companion rebinds that call to `lmFit_parallel()` and runs the vendor
unchanged. It inherits `lmFit`'s curve rather than having its own. It pays off on absolute
cost: the vendor measures 8.5 s at 948 samples, 44.2 s at 3,000, and does not finish
inside ten minutes at 9,493.

**Deliberately not parallelised.** Leaving these serial is what keeps the output identical.
`voom`'s `lowess` trend takes its span as a *fraction* of gene count, so a block fits a different
curve. `eBayes`, `squeezeVar` and `fitFDist` pool across all genes. `p.adjust` is
length-dependent, so `topTable` and `decideTests` cannot be blocked. `arrayWeights`,
`normalizeBetweenArrays` and `normalizeQuantiles` reduce across genes by construction.
`contrasts.fit` takes 0.001 s, so a fork costs more than the work. `lmFit(method = "robust")` and
`ndups >= 2` reshape the row axis through `unwrapdups`, so a split would pair different genes;
both are refused with an error, which is also why `block` is required for
`duplicateCorrelation_parallel`.

**Two edgeR companions were built, measured and deleted.** `estimateDisp` reached 1.4x but moved
19,999 of 20,000 tagwise dispersions once one library was over-sequenced. `glmQLFit` reached
1.6x, then failed on real TCGA data at 22 of 18,270 genes: `mglmLevenberg` records a deviance
whose value depends on which genes share the block when a fit does not converge, with no flag
set. A modest speedup does not buy a companion that returns different numbers.

**Three ways a limma row split goes wrong,** each reproduced against limma 3.62.2 before being
guarded, each returning a wrong answer quietly rather than failing. `asMatrixWeights` dispatches
on the *block's* row count, so a per-array weight vector is read as per-gene weights by any block
whose height equals the sample count. `NoProbeWts` is an AND-reduction selecting between two
numerically different algorithms, so an all-finite block flips to the fast path while the whole
matrix took the slow one. `stats::lm.fit` drops a one-column response to a vector, so a one-gene
block swaps `colMeans` for `mean`. Every chunk now carries at least two genes, weights are
expanded once against the whole matrix, and the branch is proven unable to flip before anything
is split.

One inner loop is the exception. It is pinned to the vendor body it came from and stands down
on any upstream change: `sva::match_quantiles`, 66.5% of ComBat-seq's serial time. See
[License](#license).

## Cross-platform

The same package on two unlike machines. Every arm on both returns `identical()` output: macOS
and Linux, fork and socket dispatch, threaded and single-threaded BLAS.

| | macOS | Linux |
|---|---|---|
| CPU | Apple M3 | 2 x Intel Xeon Silver 4208 @ 2.10 GHz |
| Physical cores | 8 (4 performance + 4 efficiency) | 16 |
| Logical | 8 | 32 |
| NUMA nodes | 1 | 2 |
| BLAS | Accelerate | OpenBLAS 0.3.20 |
| Workers swept | 2, 4, 8 | 4, 8, 16, 32 |

Absolute wall clock at cohort scale, every companion. Vendor is the serial original; best is the
fastest companion arm on that machine, with its worker count. Bold marks the faster platform.

| stage | macOS vendor | macOS best | Linux vendor | Linux best |
|---|---:|---:|---:|---:|
| ComBat-seq | 1,653.9 s | **308.0 s** (8w) | 3,355.6 s | 459.3 s (16w) |
| duplicateCorrelation | 660.5 s | 182.0 s (8w) | 663.8 s | **110.6 s** (16w) |
| calcNormFactors (TMM) | 10.3 s | **1.5 s** (8w) | 21.6 s | 5.5 s (8w) |
| lmFit | 4.4 s | **1.7 s** (4w) | 11.6 s | 4.1 s (16w) |
| removeBatchEffect | 4.5 s | **1.3 s** (8w) | 2.6 s | 2.1 s (4w) |

**Speedup is a ratio, not a speed.** The M3 runs the serial baseline 2.03x faster than the Xeon
and holds about 2x at every matched worker count. Linux scales further, 7.31x against 5.37x,
because it has twice the physical cores, and still finishes behind on the clock: 459.3 s against
308.0 s. A slower core inflates the speedup column; read the seconds.

**Hyperthreads are not cores.** 32 workers is slower than 16 on this box for every companion:
485.4 s against 459.3 s at cohort scale, 33.2 s against 31.4 s simulated.

**Where each wins.** Linux takes `duplicateCorrelation` outright. It is also the only stage
whose two vendor baselines sit within 0.5% of each other, so scaling decides it rather than
per-core speed. Everything else follows the core.

**Two Linux notes.** OpenBLAS is multi-threaded by default and reads `OPENBLAS_NUM_THREADS` when
it loads, so pin it before R starts or every forked worker opens its own thread pool. Use
`mclapply`: fork shares the parent's pages copy-on-write.

The BLAS is the obvious suspect for the platform gap, and it is not the cause. Unpinning
OpenBLAS moves DGEMM throughput 8x, from 28 to 225 GFLOPS, and moves `lmFit` by 1.3%: 8.82 s
against 8.71 s on 18,270 by 1,500. limma solves a small QR per gene and never reaches the size
where threaded BLAS pays. Per-core speed separates the machines.

## Verify it yourself

Check the claim in ten lines:

```r
set.seed(1)
counts <- matrix(rnbinom(1600, mu = 50, size = 5), nrow = 200)
batch  <- rep(1:2, each = 4)

identical(ComBat_seq_parallel(counts, batch, group = NULL, workers = 4L),
          sva::ComBat_seq(counts, batch, group = NULL))
#> TRUE
```

- **Rendered reports**, same sections on both platforms: cohort, argument parity, corrections,
  PCA, differential expression, worker sweep.
  [macOS](https://genomerx.github.io/RNA-Parallel/) ·
  [Linux](https://genomerx.github.io/RNA-Parallel/linux.html)
- **[run_example.R](inst/examples/run_example.R)** runs the whole claim in minutes without
  downloading anything.
  ```sh
  Rscript inst/examples/run_example.R
  ```
- **Report sources.** [RNA_Parallel.Rmd](inst/examples/RNA_Parallel.Rmd) for macOS,
  [RNA_Parallel_linux.Rmd](inst/examples/RNA_Parallel_linux.Rmd) for Linux, rendered with
  [render_linux.sh](inst/examples/render_linux.sh), which sets the BLAS environment before R
  starts. The first run downloads HNSC, LUAD and LUSC to the per-user cache, or
  `RNAPARALLEL_TCGA_DIR` if set.
- **[tests/](tests/testthat)** covers every argument path, chunk layout and backend, plus a
  dispatch count through the public entry point, because `identical()` alone cannot tell a
  working parallel layer from a dead one. Over 400 assertions.

## Tuning

| knob | default | change it when |
|---|---|---|
| `workers` | `min(8, detectCores() - 2)` | rarely. Going past your *physical* core count can be slower, as the 32-worker Linux arm shows. See `?workers` |
| `chunks` | `workers` | only to cut peak memory per worker |
| `parallel_backend` | `"mclapply"` | you cannot fork, or a cluster is already running |

`workers` is not a safety ceiling. Six of them alongside a second forking R session have
kernel-panicked a 24 GB machine, and nothing here can see that session.

**Backends** are `"mclapply"` (forks), `"future"`, `"BiocParallel"`, `"foreach"`, `"serial"`, or
any `function(idx, f, workers)`, all returning identical results. Forking is the default because
a forked worker reads the matrix through copy-on-write, while socket backends re-serialise it
per chunk and measured slower than not parallelising at all. **Windows** cannot fork, so `mclapply`
runs serially and says so once, `BiocParallel` substitutes a serial parameter, and `"foreach"`
runs over PSOCK but is slower than serial here. **Nesting** is safe on the default backend, where
`mc.allow.recursive` is `FALSE`; `foreach` builds its own PSOCK cluster and does multiply, so
spend the worker budget *inside* a loop body rather than across it. Inverting one 15-cohort
screen measured 245 s to 81.6 s.

**Size gates** are why a small input shows no speedup. Seven `combat.min` options set the size
below which a call runs serially, because under them a fork costs more than the work:
`combat.min.cells`, `combat.min.disp.cells`, `combat.min.glm.cells`, `combat.min.ls.cells`,
`combat.min.norm.cells`, `combat.min.order.cells`, `combat.min.dupcor.cells`.
`options(combat.fork = FALSE)` forces serial everywhere, and `combat_cluster_stop()` releases
cached clusters.

**`calcNormFactors_parallel()`** wraps `normLibSizes` on current edgeR and `calcNormFactors` on
older ones. `normLibSizes` **errors** on negative counts where the old name returned NaN-warned
factors.

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
