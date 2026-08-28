# RNA-Parallel

Parallel companions for RNA-seq tools. Each calls the original function unmodified and returns
output `identical()` to it, bit for bit. Same arguments, same defaults, same result, faster.

Rendered analysis: [macOS](https://genomerx.github.io/RNA-Parallel/) ·
[Linux](https://genomerx.github.io/RNA-Parallel/linux.html) ·
[Windows](https://genomerx.github.io/RNA-Parallel/windows.html)

## Speedup

TCGA, 18,270 genes by 1,500 tumours. Every arm `identical()` to the vendor on all three platforms.

| companion | runs | macOS<br>M3, 4P+4E | Linux<br>2x Xeon, 16 cores | Windows<br>Ultra 9 185H, 6P+10E |
|---|---|---:|---:|---:|
| `ComBat_seq_parallel()` | `sva::ComBat_seq` | 5.37x @ 8w | **7.31x @ 16w** | 2.92x @ 6w |
| `calcNormFactors_parallel()` | `edgeR::normLibSizes` | **6.86x @ 8w** | 3.89x @ 8w | 1.63x @ 2w |
| `lmFit_parallel()` | `limma::lmFit` | 2.66x @ 4w | **2.80x @ 16w** | serial |
| `duplicateCorrelation_parallel()` | `limma::duplicateCorrelation` | 3.63x @ 8w | **6.00x @ 16w** | 2.21x @ 8w |
| `removeBatchEffect_parallel()` | `limma::removeBatchEffect` | **3.37x @ 8w** | 1.22x @ 4w | serial |

Originals are timed on both sides of the companion arms and averaged; short stages are the best
of three. Bold marks the fastest platform. Machine specs and what separates them are in
[Cross-platform](#cross-platform).

**Windows is the odd column and the reason is architectural, not incidental.** It has no
`fork()`, so a worker is a whole process that receives a copy rather than sharing the parent's
pages, and only 6 of its 16 cores are performance cores. Both facts pull the useful worker count
down: the ComBat-seq curve peaks at 6 and turns over after. `serial` is not a missing measurement
— it is the companion declining to split, because `lmFit` is cheap enough per cell that no
dispatch repays the transfer there. It measured 0.52x before that gate closed, so the honest
Windows contribution for those two is parity, not a speedup.

Read the seconds as well as the ratios. Windows has the fastest vendor baseline of the three
(1,506 s against Linux's 3,356 s), and speedup rewards a slow starting point — which is why the
dual Xeon posts 7.31x while finishing *behind* the M3 in wall clock. At 515 s the laptop lands
within 12% of the Xeon's best.

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

Rendered reports carry it at cohort scale, same sections on all three platforms: cohort,
argument parity, corrections, PCA, differential expression, worker sweep.
[macOS](https://genomerx.github.io/RNA-Parallel/) ·
[Linux](https://genomerx.github.io/RNA-Parallel/linux.html) ·
[Windows](https://genomerx.github.io/RNA-Parallel/windows.html). Sources are
[RNA_Parallel.Rmd](inst/examples/RNA_Parallel.Rmd),
[RNA_Parallel_linux.Rmd](inst/examples/RNA_Parallel_linux.Rmd) and
[RNA_Parallel_windows.Rmd](inst/examples/RNA_Parallel_windows.Rmd), the second rendered with
[render_linux.sh](inst/examples/render_linux.sh) and the third with
[render_windows.ps1](inst/examples/render_windows.ps1), both of which pin BLAS before R starts.
The Windows report adds a backend sweep the other two do not need, because without `fork()` the
backend decides whether anything runs in parallel at all. The first run
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

The same package on three unlike machines. Every arm on all three returns `identical()` output:
macOS, Linux and Windows; fork and socket dispatch; threaded and single-threaded BLAS.

| | macOS | Linux | Windows |
|---|---|---|---|
| CPU | Apple M3 | 2 x Intel Xeon Silver 4208 @ 2.10 GHz | Intel Core Ultra 9 185H |
| Physical cores | 8 (4 performance + 4 efficiency) | 16 | 16 (6 performance + 8 efficiency + 2 low-power) |
| Logical | 8 | 32 | 22 |
| NUMA nodes | 1 | 2 | 1 |
| BLAS | Accelerate | OpenBLAS 0.3.20 | reference, pinned to 1 thread |
| `fork()` | yes | yes | **no** |
| Dispatch | forked children | forked children | PSOCK processes |
| Backend | `mclapply` | `mclapply` | `future` (`multisession`) |
| Workers swept | 2, 4, 6, 8 | 4, 8, 16, 32 | 2, 4, 6, 8 |

Absolute wall clock at cohort scale, every companion. Vendor is the serial original; best is the
fastest companion arm on that machine, with its worker count. Bold marks the fastest platform.

| stage | macOS vendor | macOS best | Linux vendor | Linux best | Windows vendor | Windows best |
|---|---:|---:|---:|---:|---:|---:|
| ComBat-seq | 1,653.9 s | **308.0 s** (8w) | 3,355.6 s | 459.3 s (16w) | 1,506.4 s | 515.3 s (6w) |
| duplicateCorrelation | 660.5 s | 182.0 s (8w) | 663.8 s | **110.6 s** (16w) | 439.9 s | 198.8 s (8w) |
| calcNormFactors (TMM) | 10.3 s | **1.5 s** (8w) | 21.6 s | 5.5 s (8w) | 16.3 s | 10.0 s (2w) |
| lmFit | 4.4 s | **1.7 s** (4w) | 11.6 s | 4.1 s (16w) | 6.4 s | 6.3 s (serial) |
| removeBatchEffect | 4.5 s | **1.3 s** (8w) | 2.6 s | 2.1 s (4w) | 2.9 s | 2.9 s (serial) |

Speedup is a ratio and rewards a slower core: Linux scales further, 7.31x against 5.37x, and
still finishes behind at 459.3 s against 308.0 s. Read the seconds. On Linux, pin
`OPENBLAS_NUM_THREADS` before R starts or every forked worker opens its own thread pool, though it
changes little here: unpinning moves DGEMM throughput 8x and `lmFit` by 1.3%, because limma solves
a small QR per gene.

Windows reads lowest on ratio and lands mid-pack on the clock, for the same reason. Its vendor
baseline is the fastest of the three -- 1,506.4 s against Linux's 3,355.6 s -- so there is less to
win back, and at 515.3 s it finishes within 12% of a dual Xeon while running on a laptop. What it
cannot do is scale as far, and that is architectural rather than incidental: no `fork()`, so a
worker is a whole process that receives a copy instead of sharing the parent's pages, and only 6
of its 16 cores are performance cores. Both push the useful worker count down, and the ComBat-seq
curve shows it directly -- 1.80x, 2.77x, **2.92x**, 2.87x at 2, 4, 6 and 8 workers, peaking on the
performance-core count and turning over after. The two `serial` entries are the companion
declining to split rather than a missing measurement: `lmFit` is cheap enough per cell that no
socket dispatch repays the transfer, measured at 0.14x on 21.6M cells and never reaching parity,
so it runs whole and matches the vendor instead of losing to it.

## Tuning

| knob | default | change it when |
|---|---|---|
| `workers` | `min(8, detectCores() - 2)`, capped at performance cores without `fork()` | rarely. Going past your *physical* core count can be slower, as the 32-worker Linux arm shows. See `?workers` |
| `chunks` | `workers` | only to cut peak memory per worker |
| `parallel_backend` | `"mclapply"` | you cannot fork, or a cluster is already running |

`workers` is not a safety ceiling. Six of them alongside a second forking R session have
kernel-panicked a 24 GB machine, and nothing here can see that session.

**Backends** are `"mclapply"` (forks), `"future"`, `"BiocParallel"`, `"foreach"`, `"serial"`, or any
`function(idx, f, workers)`, all returning identical results. Forking is the default because a
forked worker reads the matrix copy-on-write, while socket backends re-serialise it per chunk and
measured slower than not parallelising at all.

**Windows** cannot fork, so the backend is the whole decision there. `mclapply` runs serially and
says so once; `BiocParallel` substitutes a serial parameter. That leaves two that really do run
workers, and they are not interchangeable: `"foreach"` measured 1.18x, 0.94x, 0.57x and 0.28x at 2,
4, 8 and 16 workers, getting *worse* with every worker added, because the cached cluster is rebuilt
whenever the requested width changes and ComBat-seq alternates widths on every dispatch. `"future"`
holds one `multisession` pool across all of them and shows the ordinary shape.

**Set a plan and the package picks `"future"` for you:**

```r
library(future); plan(multisession, workers = 6)   # that is all
```

The default backend is resolved per call, so an active plan selects `"future"` and no plan leaves
`mclapply` and its one-time serial notice. This package will not set a plan for you — a caller's
plan is theirs, and starting worker processes inside someone's session unasked is worse than being
slow — which is also why the default is not simply `"future"`: without a plan that would mean a
warning on every dispatch and no speedup at all.

**Nesting** is blocked on every backend: a dispatch that finds itself already inside one of this
package's workers runs serially. On fork that used to rest on `mc.allow.recursive = FALSE`, which
covered only that branch — over PSOCK the same construction spawned workers + workers^2 processes.
A caller's own parallel loop is unaffected, so spend the worker budget *inside* a loop body rather
than across it. Inverting one 15-cohort screen measured 245 s to 81.6 s.

**Size gates** are why a small input shows no speedup. Seven options set the cell count below which
a call runs serially: `combat.min.cells` (20,000), `combat.min.disp.cells` (30,000),
`combat.min.glm.cells` (100,000), `combat.min.ls.cells` (6e6), `combat.min.norm.cells` (2e5),
`combat.min.order.cells` (4e6), `combat.min.dupcor.cells` (5,000).

Those are *fork* break-evens, and two of them move without fork, in opposite directions because
the underlying work is not alike. `lmFit` is cheap enough per cell that a serialised chunk is
never repaid — measured 0.14x at 21.6M cells and 0.24x at 60M, improving with size and never
reaching parity — so `combat.min.ls.cells` **closes** on Windows rather than merely rising, which
also covers `removeBatchEffect_parallel()` since the vendor rebinds one `lmFit` call. The TMM
column loop does pay once it is big enough — 1.05x at 1.8M cells against 1.58x at 21.6M — so
`combat.min.norm.cells` **rises** to 2e6 instead, an order of magnitude above the fork value.
Set either explicitly to override. Output is `identical()` either way: a gate decides who
computes, never what is computed.
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
