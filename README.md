# RNA-Parallel

[![version](https://img.shields.io/badge/version-0.5.0-blue)](NEWS.md)
[![R](https://img.shields.io/badge/R-%E2%89%A5%204.1-blue)](DESCRIPTION)
[![license](https://img.shields.io/badge/license-MIT-blue)](LICENSE.md)

Parallel companions for RNA-seq tools. Each calls the original function unmodified and returns
output `identical()` to it, bit for bit. Same arguments, same defaults, same result, faster.

Rendered analysis: [macOS](https://genomerx.github.io/RNA-Parallel/) ·
[Linux](https://genomerx.github.io/RNA-Parallel/linux.html) ·
[Windows](https://genomerx.github.io/RNA-Parallel/windows.html)

## Speedup

TCGA, 18,270 genes by 1,500 tumours. Every arm `identical()` to the original on all three platforms.

macOS: M3, 4P+4E. Linux: 2x Xeon, 16 cores. Windows: Ultra 9 185H, 6P+10E.

| companion | runs | macOS | Linux | Windows |
|---|---|---:|---:|---:|
| `ComBat_seq_parallel()` | `sva::ComBat_seq` | 5.43x @ 8w | **9.34x @ 16w** | 3.57x @ 6w |
| `calcNormFactors_parallel()` | `edgeR::normLibSizes` | **6.78x @ 8w** | 4.12x @ 8w | 1.74x @ 2w |
| `lmFit_parallel()` | `limma::lmFit` | 3.02x @ 8w | **3.37x @ 16w** | 1.10x @ 4w |
| `duplicateCorrelation_parallel()` | `limma::duplicateCorrelation` | 3.79x @ 6w | **7.22x @ 16w** | 2.07x @ 4w |
| `removeBatchEffect_parallel()` | `limma::removeBatchEffect` | **2.97x @ 6w** | 1.39x @ 8w | 0.90x @ 2w |

Originals are timed on both sides and averaged (best of three for short stages). Bold = fastest
platform per row. Full specs and the "why" in [Cross-platform](#cross-platform).

**Windows has no `fork()`.** A worker is a whole copied process, not a shared-memory fork, and
only 6 of its 16 cores are performance cores. Both cap the useful worker count. ComBat-seq shows
it directly: 2.12x, 3.34x, **3.57x**, 3.54x at 2/4/6/8 workers, peaking at the P-core count.

**`lmFit` and `removeBatchEffect` on Windows are not speedups (1.10x, 0.90x).** Their size gates
close on this platform, so the companion just runs the original plus wrapper overhead. Call that
parity, not a result.

Ratios don't predict wall clock: the dual Xeon wins 4 of 5 rows but still loses ComBat-seq on time
to the M3 (325.9s vs 298.6s) because it starts from a slower original (3,043.8s vs 1,619.8s).
Windows is last on both counts, 720.1s vs 298.6s.

**Nothing is reimplemented.** The original function runs, called with hot paths rebound in a child
of its own environment. Every other symbol still resolves to original code. `identical()` is
asserted, not a tolerance.

## Install

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("sva", "edgeR", "limma"))

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("GenomeRx/RNA-Parallel")
```

## Use

Each companion takes its original's arguments in the same order with the same defaults, and adds
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

One interface difference: `duplicateCorrelation_parallel` requires `block`, whereas the original
defaults it to `NULL`. [How it works](#how-it-works) says why.

### Check it yourself

All five companions against their originals. No download, nothing to configure:

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

Four thousand cells sits under every size gate, so this runs serially and is still `identical()`.
[run_example.R](inst/examples/run_example.R) runs the same check at a size where speedup shows.

```sh
Rscript inst/examples/run_example.R
```

Rendered reports carry it at cohort scale, same sections on all three platforms:
[macOS](https://genomerx.github.io/RNA-Parallel/) ·
[Linux](https://genomerx.github.io/RNA-Parallel/linux.html) ·
[Windows](https://genomerx.github.io/RNA-Parallel/windows.html). Sources in
[inst/examples/](inst/examples/), rendered with each platform's own script, which pins BLAS
before R starts (a BLAS reads its thread count at load time, too late to set in-session). Windows
also sweeps backends, since without `fork()` the backend decides whether anything runs parallel at
all. First run downloads HNSC/LUAD/LUSC to the per-user cache, or `RNAPARALLEL_TCGA_DIR` if set.

[tests/](tests/testthat) covers every argument path, chunk layout, backend, and dispatch count
through the public entry point. `identical()` alone can't tell a working parallel layer from a
dead one. 400+ assertions.

### Seeing what is running

```r
options(combat.timing = TRUE)      # one elapsed line per call
options(combat.quiet  = TRUE)      # swallow the original's progress chatter
options(combat.timing.min = 1)     # hide anything faster than a second
```

```
  Breast ComBat-seq                  mclapply x6         23.2s
  TMM 12,000 x 700                   mclapply x6          0.9s
  TMM 200 x 8                        serial               0.0s  2 gated
```

The engine column reports what **ran**, not what was asked for: a call that fell under a size
gate says `serial`, because a companion returning `identical()` output at serial pace otherwise
looks exactly like one that forked. If a pinned original excerpt stands down after an upstream
change, the line says `[match_quantiles stood down]`.

**`rnaparallel_stale()`** returns TRUE if the package was reinstalled under a running session.
The fix is to restart R.

## How it works

`ComBat_seq_parallel()` creates a child of the original's environment, binds six names, and calls
`sva::ComBat_seq` unchanged. Shares below are of serial time on 10,000 genes by 500 samples.

| ComBat-seq calls | share | split by | why exact |
|---|---|---|---|
| `match_quantiles` | 66.5% | gene rows | each output cell reads only its own gene |
| `estimateGLMTagwiseDisp` | 14.0% | gene rows | valid at `prior.df = 0`; other values go to `edgeR` unsplit |
| `estimateGLMCommonDisp` | 13.1% | batches | sums over all genes, so batches split independently, RNG-free |
| `glmFit`, `glmFit.default` | remainder | gene rows | offset/dispersion/weights/start slice with the rows |
| `monte_carlo_int_NB` | small | serial | draws depend on the previous batch's state |

`calcNormFactors_parallel()` splits by sample column, `lmFit_parallel()`/`duplicateCorrelation_parallel()`
by gene row. `removeBatchEffect_parallel()` splits nothing: it rebinds the one `lmFit` call
inside the original and inherits its curve (still pays: the original doesn't finish in ten minutes
at 9,493 samples). Chunks interleave so a sorted matrix doesn't idle a worker; a dead worker,
duplicate chunk, or short result halts the run.

**Stays serial → output stays identical:** anything reducing across genes (`eBayes`, `squeezeVar`,
`fitFDist`, `arrayWeights`, `normalizeBetweenArrays`, `normalizeQuantiles`) or scaling with gene
count (`voom`'s `lowess` span, `p.adjust`, `topTable`, `decideTests`). `contrasts.fit` takes
0.001s, less than a fork. `lmFit(method = "robust")` and `ndups >= 2` reshape rows and error
rather than split, which is why `block` is required for `duplicateCorrelation_parallel`.

**Built, measured, deleted:** modest speedup doesn't buy different numbers.

| companion | reached | what it moved |
|---|---|---|
| `estimateDisp` | 1.4x | 19,999 of 20,000 tagwise dispersions, once one library was over-sequenced |
| `glmQLFit` | 1.6x | 22 of 18,270 TCGA genes, non-convergent deviance depending on block membership, no flag set |

**Three limma row-split traps,** reproduced against limma 3.62.2, each returning a wrong answer
quietly:

| trap | wrong answer | guard |
|---|---|---|
| `asMatrixWeights` dispatches on *block's* row count | per-array weights read as per-gene when block height = sample count | weights expanded once against the whole matrix |
| `NoProbeWts` is an AND-reduction between two algorithms | an all-finite block takes the fast path the whole matrix didn't | branch proven unable to flip before split |
| `stats::lm.fit` drops a one-column response to a vector | one-gene block swaps `colMeans` for `mean` | at least two genes per chunk |

One inner loop is pinned to the original body and stands down on any upstream change:
`sva::match_quantiles`, 66.5% of ComBat-seq's serial time. See [License](#license).

## Cross-platform

The same package on three unlike machines. Every arm on all three returns `identical()` output:
macOS, Linux and Windows; fork and socket dispatch; threaded and single-threaded BLAS.

| | macOS | Linux | Windows |
|---|---|---|---|
| CPU | Apple M3 | 2 x Intel Xeon @ 2.80 GHz | Intel Core Ultra 9 185H |
| Physical cores | 8 (4 performance + 4 efficiency) | 16 | 16 (6 performance + 8 efficiency + 2 low-power) |
| Logical | 8 | 32 | 22 |
| NUMA nodes | 1 | 2 | 1 |
| BLAS | Accelerate | OpenBLAS 0.3.8 | reference, pinned to 1 thread |
| `fork()` | yes | yes | **no** |
| Dispatch | forked children | forked children | PSOCK processes |
| Backend | `mclapply` | `mclapply` | `future` (`multisession`) |
| Workers swept | 2, 4, 6, 8 | 2, 4, 8, 16 | 2, 4, 6, 8 |

Absolute wall clock at cohort scale. Original is the serial upstream call; best is the fastest
companion arm, with its worker count. Bold = fastest platform per stage.

| stage | macOS | Linux | Windows |
|---|---:|---:|---:|
| ComBat-seq | 1,619.8s → **298.6s** (8w) | 3,043.8s → 325.9s (16w) | 2,573.3s → 720.1s (6w) |
| duplicateCorrelation | 663.7s → 175.0s (6w) | 469.1s → **64.9s** (16w) | 784.8s → 379.3s (4w) |
| calcNormFactors (TMM) | 10.3s → **1.5s** (8w) | 19.5s → 4.7s (8w) | 26.5s → 15.2s (2w) |
| lmFit | 4.4s → **1.5s** (8w) | 11.4s → 3.4s (16w) | 10.1s → 9.2s (4w) |
| removeBatchEffect | 4.5s → **1.5s** (6w) | 2.5s → 1.8s (8w) | 4.4s → 4.9s (2w) |

Speedup is a ratio and rewards a slower core: Linux scales further, 9.34x against 5.43x, and
still finishes behind at 325.9s vs 298.6s. On Linux, pin `OPENBLAS_NUM_THREADS` before R starts,
otherwise every forked worker opens its own thread pool (unpinning moves DGEMM 3.6x, 61→216
GFLOPS, though `lmFit` barely changes since limma solves a small QR per gene).

Windows reads lowest on ratio and last on the clock, but its original baseline (2,573.3s) sits
between the M3 (1,619.8s) and Xeon (3,043.8s), so this isn't a short-baseline story, it's just
slower to finish (720.1s vs 298.6s/325.9s). What caps it is architectural: no `fork()` means a
worker is a full copied process, and only 6 of 16 cores are performance cores. Both cap the useful
worker count. ComBat-seq peaks at 3.57x at 6 workers, exactly the P-core count, then turns over.

`lmFit` and `removeBatchEffect` measure 1.10x/0.90x on Windows. Neither is a speedup. Both size gates
close where the payload gets copied, so the companion just runs the original plus wrapper
overhead, landing either side of 1.00x on noise rather than climbing with workers. (Left to split,
`lmFit` measured 0.14x at 21.6M cells and never reached parity, which is why the gate exists.)

## When a companion is worth reaching for

Not every companion pays at every size, and the honest answer is per companion rather than one
rule. Measured on the M3, at the default worker count, every arm `identical()` to its original.
The ratio is companion against original, so below 1.00x the companion is the slower of the two.

| companion | smallest measured | crossover | large |
|---|---|---:|---:|
| `duplicateCorrelation_parallel()` | 1.68x at 200 x 12 | pays at every size measured | 9.52x at 3,000 x 100 |
| `calcNormFactors_parallel()` | 1.44x at 2,000 x 20 | pays at every size measured | 6.26x at 20,000 x 500 |
| `ComBat_seq_parallel()` | 0.69x at 300 x 20 | about 1,000 genes | 5.37x on the cohort |
| `lmFit_parallel()`, voom or probe weights | 0.59x at 1,000 x 24 | about 4,000 genes | 2.79x at 60,000 x 48 |
| `lmFit_parallel()`, no probe weights | 1.00x | 6M cells | splits only above the gate |
| `removeBatchEffect_parallel()` | 0.92x at 20,000 x 50 | 6M cells | 1.91x at 20,000 x 500 |

Read it as three groups: **`duplicateCorrelation`/`calcNormFactors` pay unconditionally** (one
gene = one REML fit; TMM's `rank` hoist wins with zero workers). **`ComBat-seq` and weighted
`lmFit` have a floor** below which fork overhead beats the savings. **`lmFit` without probe
weights, and `removeBatchEffect`, are parity until the input is large:** limma vectorises every
gene into one `lm.fit` (milliseconds), and `removeBatchEffect` is just that call plus a BLAS
product.

Choosing wrong costs little: under-gate is one original call plus ~0.3-0.6ms overhead, real only
in a loop over thousands of small units, since gates decide per-call and can't see the loop.

**To check what actually ran:** `options(combat.timing = TRUE)` prints the engine per call:

```
  ComBat-seq 18,270 x 1,500          mclapply x6         308.0s
  lmFit 20,000 x 24                  serial                0.0s  3 gated
```

`serial` with a gated count means every dispatch fell under a threshold and the companion was a
pass-through. On that input, call the original.

## Tuning

| knob | default | change it when |
|---|---|---|
| `workers` | `min(8, detectCores() - 2)`, capped at performance cores without `fork()` | rarely. Going past your *performance* core count can be slower, as the Windows curve shows: 2.12x, 3.34x, **3.57x**, 3.54x at 2, 4, 6 and 8 workers on a chip with six. See `?workers` |
| `chunks` | `workers` | only to cut peak memory per worker |
| `parallel_backend` | `"mclapply"` | you cannot fork, or a cluster is already running |

`workers` is not a safety ceiling. Six workers alongside a second forking R session have
kernel-panicked a 24 GB machine; nothing here can see that other session.

**Backends:** `"mclapply"` (forks, default), `"future"`, `"BiocParallel"`, `"foreach"`, `"serial"`,
or any `function(idx, f, workers)`, all return identical results. Forking wins by default because
a forked worker reads the matrix copy-on-write; socket backends re-serialise per chunk and
measured slower than not parallelising.

**Windows can't fork**, so backend choice is the whole decision. `mclapply` runs serially (says so
once); `BiocParallel` substitutes a serial param. Of the two that actually parallelize, they're not
interchangeable: `"foreach"` gets *worse* with more workers (1.18x→0.28x at 2→16w) because its
cached cluster rebuilds on every width change and ComBat-seq alternates widths per dispatch.
`"future"` holds one `multisession` pool across all of them and scales normally.

**Set a plan and the package picks `"future"` for you:**

```r
library(future); plan(multisession, workers = 6)   # that is all
```

Backend resolves per call: an active plan selects `"future"`; no plan leaves `mclapply` (one-time
serial notice). The package never sets a plan itself. A caller's plan is theirs, and starting
workers inside someone's session unasked is worse than being slow.

**Nesting is blocked on every backend.** A dispatch already inside one of this package's workers
runs serially (via `mc.allow.recursive = FALSE` on fork; an equivalent guard over PSOCK, which
used to spawn workers² processes without it). Your own loop is unaffected either way. Spend the
worker budget *inside* a loop body, not across it. Inverting one 15-cohort screen measured 245s →
81.6s.

**Size gates** are why small inputs show no speedup. Nine options set the cell/gene count below
which a call runs serially: `combat.min.cells` (20,000), `combat.min.disp.cells` (30,000),
`combat.min.glm.cells` (100,000), `combat.min.ls.cells` (6e6), `combat.min.norm.cells` (2e5),
`combat.min.order.cells` (4e6), `combat.min.dupcor.cells` (5,000), `combat.min.batch.cells`
(20,000), `combat.min.wt.genes` (2,000, counted in genes not cells: limma's weighted branch is a
per-gene loop that barely moves with array size).

Two gates move without fork, in opposite directions: `lmFit` never reaches parity on Windows
(0.14x at 21.6M cells, 0.24x at 60M) so `combat.min.ls.cells` **closes** there (covers
`removeBatchEffect_parallel()` too, since it rebinds one `lmFit` call). TMM pays once big enough
(1.05x at 1.8M cells → 1.58x at 21.6M) so `combat.min.norm.cells` **rises** to 2e6. Override either
explicitly. Output is `identical()` regardless: a gate decides who computes, never what.
`options(combat.fork = FALSE)` forces serial everywhere; `combat_cluster_stop()` releases cached
clusters.

**`calcNormFactors_parallel()`** wraps `normLibSizes` on current edgeR, `calcNormFactors` on older
ones. `normLibSizes` **errors** on negative counts where the old name returned NaN-warned factors.

## License

MIT for this companion, copyright GenomeRx 2026, in [LICENSE](LICENSE). The original packages are
called at run time from your own installation and none of them is redistributed here: sva is
Artistic-2.0, limma and edgeR are GPL (>= 2).

One block is the exception and is marked as such in the source. `R/helper_seq_parallel.R` carries
the deparsed `sva::match_quantiles` body and a row-vectorised transcription of it, so the
companion can detect an upstream change and stand down. It derives from Artistic-2.0 code by
Zhang, Parmigiani and Johnson. Nothing else here reproduces original code.

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
