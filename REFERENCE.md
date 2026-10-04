# Reference

Deep detail split out of the README: internals, full cross-platform numbers, and tuning. Read
this when the README's summary is not enough.

## Glossary

- **companion**: a function such as `lmFit_parallel()` that runs one original function, here
  `limma::lmFit`, and returns the same output, faster.
- **worker**: one extra R process that does part of the work.
- **chunk**: the block of gene rows (or sample columns) one worker handles at a time.
- **dispatch**: handing chunks out to workers and collecting the results.
- **fork**: on macOS and Linux, the system clones the running R process. The clone shares memory
  with the parent, so no data is copied to start a worker.
- **backend**: the mechanism used to run workers (`mclapply`, `future`, `BiocParallel`, `foreach`,
  or serial).
- **size gate**: an input-size threshold. Below it, the companion runs the original code
  unchanged, since forking a worker costs more than the split saves at that size.
- **cells**: genes times samples. Most size gates count in this unit; `combat.min.wt.genes`
  counts genes.
- **parity**: running at 1.0x, the same speed as the original.

## When a companion is worth reaching for

Measured on an Apple M3, at the default worker count. Ratio is companion speed against the
original function, same input, same output.

| companion | runs | smallest size tested (genes x samples) | size where it turns faster | largest size tested |
|---|---|---|---:|---:|
| `duplicateCorrelation_parallel()` | `limma::duplicateCorrelation` | 1.68x at 200 x 12 | faster at every size tested | 9.52x at 3,000 x 100 |
| `calcNormFactors_parallel()` | `edgeR::normLibSizes` | 1.44x at 2,000 x 20 (TMM) | TMM: faster at every size tested; other methods: only above their column gate | 6.26x at 20,000 x 500 (TMM) |
| `ComBat_seq_parallel()` | `sva::ComBat_seq` | 1.20x at 300 x 20, nothing dispatched | faster at every size tested | 4.73x at the default 6 workers on the TCGA cohort, 18,270 x 1,500 (5.43x at 8) |
| `lmFit_parallel()`, voom or probe weights | `limma::lmFit` | 1.00x below 2,000 genes (one original call) | 2,000 genes (`combat.min.wt.genes`) | 3.70x at 60,000 x 48 ([tools/bench_lmfit.R](tools/bench_lmfit.R)) |
| `lmFit_parallel()`, no probe weights | `limma::lmFit` | 1.00x | 6M cells | splits only above the gate |
| `removeBatchEffect_parallel()` | `limma::removeBatchEffect` | 1.00x at 2,000 x 20 | 6M cells | 1.91x at 20,000 x 500 |

`duplicateCorrelation_parallel()` and `calcNormFactors_parallel()` with `method = "TMM"` are
faster at every size tested, for reasons that need no workers. `duplicateCorrelation_parallel()`
memoises the `La.svd` call that statmod repeats with identical arguments for every gene, which
measured 47.5% of the original's time at 300 genes. The 1.68x at 200 x 12 comes from that memo:
2,400 cells sit under `combat.min.dupcor.cells` (5,000), so nothing is dispatched. Above that
gate it also splits genes across workers. For `method = "TMM"`, `calcNormFactors_parallel()`
computes each column's two rank vectors once instead of twice. `TMMwsp`, `RLE` and
`upperquartile` get no such step, so below `combat.min.norm.cells` or `combat.min.order.cells`
they run the original plus a fixed overhead, measured at about 0.4x to 0.8x at 2,000 x 20.

`ComBat_seq_parallel()` gates every dispatch on a small input and is still faster there, because
its row-vectorized quantile match beats the original's cell loop in one process.
`lmFit_parallel()` with probe weights makes one plain original call below 2,000 genes, so it
runs at parity there. From 2,000 genes (and 20,000 cells) its per-gene loop splits, measured
1.39x to 1.75x at 2,000 genes with the gate forced open. `lmFit_parallel()` without probe
weights, and `removeBatchEffect_parallel()`, run at the original's speed until the input is
large: limma fits every gene in one vectorized `lm.fit()` call that already takes milliseconds,
and `removeBatchEffect` is that call plus one matrix product.

Using a companion on data that is too small costs one original call plus about 0.3 to 0.6 ms.
That only adds up if it runs thousands of times in a loop, because each call checks its own size
gate and cannot see the loop around it.

## Seeing what is running

```r
options(combat.timing = TRUE)      # one elapsed line per call
options(combat.quiet  = TRUE)      # swallow the original's progress chatter and the progress bar
options(combat.timing.min = 1)     # hide anything faster than a second
```

```
  Breast ComBat-seq                        mclapply x6         23.2s
  calcNormFactors TMM 12,000 x 700         mclapply x6          0.9s
  calcNormFactors TMM 200 x 8              serial               0.0s  2 gated
```

The engine column reports what actually ran, not what was asked for. A call under a size gate
says `serial`, because `identical()` output at serial speed otherwise looks exactly like a real
parallel run. If a pinned copy of upstream code stands down after an upstream change, the line
says `[match_quantiles stood down]`. Where the process's peak memory can be read (`/proc` on
Linux, the ps package elsewhere), the line also shows `peak N GB` after the elapsed time.

Warnings and messages raised inside workers reach the caller once each, in job order, as the
original emits them, on every backend, so `suppressWarnings()` and `suppressMessages()` work on
them as usual. When one companion calls another, as `removeBatchEffect_parallel()` calls
`lmFit_parallel()`, the memory guard still runs once per companion call you make, so it warns at
most once, and the timing line names the worker count the nested call dispatched with.

**Progress is on by default.** Every companion call prints a running count of its dispatches on
one line, redrawn in place up to 4 times a second (`  ComBat-seq 18,270 x 1,500: 42
dispatched`). `combat.timing` alone reports only once a call finishes, and ComBat-seq dispatches
its slow steps up to `2 * n_batch + 4` times per call, so on a cohort of hundreds of batches a
slow run would otherwise look the same as a stuck one. On macOS and Linux, a dispatch that really
runs in parallel also forks a small reporter process that draws a bar from the progress files
below while the main process waits on its workers. The count stands aside while the bar is up.
No reporter starts for a gated or serial dispatch, inside a worker, under
`options(combat.fork = FALSE)`, on Windows, or outside an interactive session or terminal; in a
knitr or Rscript log, set `combat.progress.dir` and follow the run with `rnaparallel_progress()`
from a second session.
Both lines clear themselves before the
`combat.timing` summary prints. `options(combat.progress = FALSE)` turns off the count and the
bar.

**Progress files.** Workers write them on every parallel dispatch, whatever `combat.progress`
says: each worker appends one line to its own file at the start and end of every chunk. No
option turns them off. Each write is one line per chunk, not per gene, so the cost does not show
up next to the compute itself even over hundreds of chunks and hours. With no directory set, each
parallel dispatch writes its files into its own subdirectory of
`file.path(tempdir(), "rnaparallel-progress")` and removes it when the dispatch returns, so only
the running session can find them, and only while that dispatch runs. They are the only view during a dispatch on Windows, where no reporter runs, and in an
RStudio Server session, where a forked process's console output does not reliably arrive.

To follow a long run from a second R session, set
`options(combat.progress.dir = "some/writable/path")` before the call. Then, from the second
session, while the run is still going:

```r
rnaparallel_progress("some/writable/path")
#> 47 done, 128 started, 0 stalled, 12.4m/chunk, ETA 16:42
```

`done`/`started`/`stalled` count chunks. `stalled` means started but never finished, which is
what a killed worker looks like. Seconds per chunk and the ETA need at least two finished chunks
to mean anything, and read `NA` until then.

A live bar reads the same directory, from that same second session:

```r
rnaparallel_progress("some/writable/path", watch = TRUE)
#> |===============================-------------------|  62%  ComBat-seq 40,609 x 9,493                79/128  ETA 16:42
```

Redraws every `interval` seconds (default 1), on one line. Returns about 2 seconds after every
started chunk has finished. Gives up after `stall_after` seconds (default 600) in which no chunk
starts or finishes, whether the run died or nobody is writing to `dir` at all. Each poll reads
only the bytes added since the last one, not the whole file again.

`rnaparallel_stale()` returns TRUE if the package was reinstalled while the current R session was
still running. The fix is to restart R.

## Memory

Forking copies nothing up front, so N workers do not cost N times the parent's memory at the
start. They cost whatever each one writes to as it runs, and for a row-split fit over a large
matrix that can be a real fraction of it. When the total exceeds what the machine has, the
kernel does not hand R a normal error: on a machine without swap, Linux kills the R process
outright, with no error message, no traceback, and `mclapply` reporting nothing. `future` only
says a future was interrupted. From outside, R just disappears.

Measured on a 40,609 x 9,493 matrix, 125 GB, no swap:

| workers | parent memory in use | outcome |
|---|---|---|
| 16 | 50 GB | killed |
| 8 | 100.8 GB | killed, 0 GB free at the fork |
| 4 | 23 GB | killed, 23 -> 111 GB in 30s |
| 2 | 23 GB | survived, 102 GB peak |

**The memory guard** (the internal `rp_mem_cap()`, not exported) runs once per companion call,
before its first dispatch. It reads how much memory is available and the main process's own
memory use. When the requested worker count would need more than 80% of what is available, it
lowers the worker count instead and warns with all three numbers. With the default
`combat.mem.divergence = 1`:

```
rnaparallel: 16 workers need ~80 GB on top of a 5 GB parent and only 10 GB is
available, which ends in an out-of-memory kill or heavy swapping, not an R error.
Using 1 instead. Set options(combat.mem.divergence=) if this workload dirties
less, or options(combat.mem.guard=FALSE) to disable.
```

It reads `/proc` on Linux and the ps package elsewhere. Without ps off Linux, or wherever neither
source answers, it reads `NA` and leaves the worker count unchanged, since it never blocks what
it cannot measure. `combat.mem.divergence` (default 1) is the fraction of the parent's memory
each worker is expected to modify, and so copy. The default of 1 assumes the worst case, a full
copy. Lower it for a workload known to modify less, for example a per-column trimmed mean.
`combat.mem.guard = FALSE` disables the whole check.

`combat.mem.chunk.cells` (unset by default) caps the cells in one chunk. Set to N, it raises the
chunk count past what `chunks` or `workers` asked for until no chunk holds more than N cells
(genes times samples), as long as N cells fit the smallest chunk a companion allows: one gene
row, one sample column for `calcNormFactors_parallel()`, or two rows for `lmFit_parallel()`. It
never lowers the chunk count, reads no memory, and does not change the output. It bounds each
chunk, not how many run at once; the worker count, and so the guard above, decides that.
ComBat-seq's dispersion jobs across batches are not chunks and hold one batch each regardless.

The guard does not read memory again between dispatches. A call whose main process grows after
it starts forks its later workers from a larger parent than the guard measured. ComBat-seq does
this: it builds its full-size working matrices (fitted means, dispersions, adjusted counts) after
entry, before its later dispatches. Size `workers` for a large ComBat-seq run with that growth in
mind.

If the main R process is killed, its forked workers keep running and keep their share of memory
until the machine restarts, since they were reparented to the operating system's init process,
not shut down. A normal kill signal does not stop them either, because R installs its own
handler and the worker is blocked mid-computation. Each forked worker therefore checks, at the
start of every chunk, whether its parent is still the master, and exits if not. Reading the
parent's process id needs the ps package (in Suggests) on R builds without `Sys.getppid()`,
R 4.4.2 among them; without it the check is skipped and orphaned workers keep running. Socket
workers skip the check and exit when their connection closes.

**`rnaparallel_set_mem_limit()`** sets `R_MAX_VSIZE`, the most vector memory one R process will
let itself use, for every future session that reads the `.Renviron` it writes. R checks this
limit every time it asks for memory, so a process that hits it fails with a normal R error
instead of a kernel kill:

```r
rnaparallel_set_mem_limit(dry_run = TRUE)   # see the computed value, write nothing
rnaparallel_set_mem_limit()                 # write it to ~/.Renviron
```

It reads total RAM (`/proc/meminfo` on Linux, PowerShell on Windows, the ps package elsewhere),
multiplies it by `fraction` (default 0.5), and rounds down to the largest of
8/16/32/64/128/256/512/1024 GB at or below that. Under 8 GB it writes the exact value in Mb.
It writes to `R_ENVIRON_USER` when that is set, otherwise `~/.Renviron`, or to a `path` passed
explicitly, which is also how to test this without touching a real file. The write only takes
effect on the next R session, since `.Renviron` is read once at startup. R reads `R_ENVIRON_USER`
when it is set, otherwise the first of `./.Renviron` and `~/.Renviron` it finds, so a project
with its own `.Renviron` needs the line there too. On Windows, `path.expand("~")` resolves via
`USERPROFILE`, not the `HOME` environment variable, which matters if you were planning to
redirect it by setting `HOME`.

The limit is per process. Each forked worker inherits it and counts its own allocations against
it, so N workers that each stay under it can still exhaust RAM together, which is the fork kill
above. It turns one process's overshoot into an R error. It does not bound the total across
workers; only the memory guard, or fewer workers, does that.

## Full self-check

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
[run_example.R](inst/examples/run_example.R) repeats the ComBat-seq comparison on 20,000 genes by
480 samples, where the speedup shows. The rendered reports compare all five companions at cohort
scale.

## How it works

`ComBat_seq_parallel()` runs the original `sva::ComBat_seq` code without changing it. It calls
the original with six symbols rebound in a child of its own environment: `match_quantiles`,
`estimateGLMTagwiseDisp`, `glmFit`, `glmFit.default`, `sapply` and `lapply`. The `sapply` and
`lapply` shims pick out the one call of each that runs a dispersion estimate per batch and
dispatch it across batches; every other `sapply` and `lapply` in the body goes to base R. Shares
below are of serial run time, measured on 10,000 genes by 500 samples with `group = NULL`.

| internal step | reached through | share of serial run time | work divided by | why the output is unchanged |
|---|---|---|---|---|
| `match_quantiles` | `match_quantiles` | 66.5% | gene rows | each output cell reads only its own gene |
| `estimateGLMTagwiseDisp` | `lapply` | 14.0% | batches | each batch's estimate reads only that batch's columns and draws no random numbers |
| `estimateGLMCommonDisp` | `sapply` | 13.1% | batches | the same; it is not split by rows, because it maximizes a sum over genes and a split sum can move the result in the last bit |
| `glmFit`, `glmFit.default` | `glmFit`, `glmFit.default` | remainder | gene rows, except one-way designs | offset, dispersion, weights and start values all slice with the rows |

The shares move with batch count: on 6,000 genes by 400 samples in 10 batches, also with
`group = NULL`, tagwise dispersion was 60% of the run and `match_quantiles` 7%.

Two row splits are rebound but often do not run. The `glmFit` split refuses a one-way design,
one whose columns are exactly the levels of one factor: edgeR fits it with a one-group kernel
whose result, for a gene that does not converge cleanly, depends on which other genes share the
matrix. With `group = NULL` the design is batch indicators only, so the fit runs whole; it splits
when `group` or `covar_mod` makes the design more than one-way. The tagwise row split runs only
with `prior.df = 0` and a design that is not one-way. ComBat-seq passes `prior.df = 0`. With
`group` alone or no covariates its per-batch design is one-way, so the tagwise estimate is
dispatched across batches and never split by rows. So in ComBat-seq the tagwise row split runs
only when `covar_mod` makes a batch's own design more than one-way, that is, gives it more
distinct rows than columns. A continuous covariate usually does. A single categorical covariate
with `group = NULL` does not, because an intercept plus its indicator columns is still a one-way
layout. Where the split runs, it runs in parallel only when the dispatch across batches runs in
the master.

Everything else in the body resolves to sva's (or the sourced upstream copy's) and edgeR's own
code, including `vec2mat`, `estimateGLMCommonDisp` itself, `getOffset` and `monte_carlo_int_NB`.
`monte_carlo_int_NB`, which runs when `shrink = TRUE`, stays serial because each batch's
`sample()` draws depend on the random-number state the previous batch left.

`calcNormFactors_parallel()` runs edgeR's `normLibSizes` (`calcNormFactors` in older releases)
and splits its per-sample loop by sample column. Each factor depends only on its own column and
on quantities computed once from all samples, such as TMM's reference column. For
`method = "TMM"` it also gives `.calcFactorTMM` a one-slot cache for `rank`, so each column's two
rank vectors are computed once instead of twice, with or without workers. It wraps edgeR's
default method only, which the `DGEList` method also ends in; any other class with its own edgeR
method runs different code, so the companion refuses it.

`lmFit_parallel()` rebinds limma's `lm.series` and `gls.series` and splits them by gene row.
`duplicateCorrelation_parallel()` runs limma's `duplicateCorrelation` on row blocks with
`mixedModel2Fit` rebound to a copy whose `La.svd` calls are memoised, and with `factor` and
`model.matrix` memoised in a child of `duplicateCorrelation`'s own environment. The pooled tail,
a trimmed mean over every gene, runs once on the concatenated result.

`removeBatchEffect_parallel()` does no splitting of its own. limma's `removeBatchEffect` is
one `lmFit` call and one matrix product. The companion rebinds that `lmFit` to `lmFit_parallel()`
and gets all of its speedup there. The matrix product stays whole, because a split product can
round differently once BLAS is threaded. The original's cost grows with sample count and with
batch count, since the batch columns join the design of that one fit.

Chunks are handed to workers in interleaved order (row 1 to worker A, row 2 to worker B, and so
on), so a matrix sorted by expression level does not leave one worker idle while the others
finish. If a worker dies, a chunk comes back twice, or a result comes back short, the run stops
rather than silently returning a wrong answer.

`lmFit_parallel()` refuses two configurations rather than running them serially:
`method = "robust"`, which fits through `mrlm`, and `ndups >= 2`, which makes `unwrapdups`
reshape the matrix so a gene no longer occupies one row. `duplicateCorrelation_parallel()`
requires `block`, unlike the original. With `block = NULL` the original pairs rows through
`unwrapdups`, and `if (spacing == "topbottom") spacing <- nrow(M)/2` branches on the row count of
whatever block it is handed, so a row split would pair different genes and raise nothing.

**Splitting a matrix by rows changes what limma sees in three ways.** Each one returns a wrong
answer with no warning. The package checks all three before it ever splits, reproduced against
limma 3.62.2:

| trap | wrong answer | guard |
|---|---|---|
| `asMatrixWeights` reads the row count of `block`, not of the full matrix, to decide what kind of weights it is looking at | per-array weights get read as per-gene weights when `block`'s height happens to match the sample count | weights are expanded against the whole matrix once, before any split |
| limma's `lm.series` and `gls.series` pick one of two fitting methods based on whether every expression value is finite after unusable weights are set to `NA` (at or below 0 in `lm.series`, below 1e-15 in `gls.series`), and whether the weights are absent or per-array | a chunk that happens to be all-finite can take a faster method the full matrix would not have qualified for | the finite check runs once on the whole matrix, before splitting |
| `stats::lm.fit` silently turns a one-column response into a plain vector | a one-gene chunk swaps a matrix function (`colMeans`) for a vector function (`mean`), changing the shape of the result | every chunk keeps at least two genes |

Two pieces of upstream code are transcribed rather than called. `match_quantiles`, 66.5% of
ComBat-seq's serial time, runs as a row-vectorized transcription of the sva 3.54.0 body, used
only while the installed sva's `match_quantiles` deparses to the pinned text byte for byte. On
any difference each slice runs sva's own function instead, and the timing line says
`[match_quantiles stood down]`. The tagwise row split computes the three default statements of
edgeR 4.4.2's `estimateGLMTagwiseDisp.default` (`offset`, `span` and `AveLogCPM`) itself, once
for the whole matrix, so every chunk sees the same values. They are pinned in
`.tagwise_defaults_pinned` and compared with `identical()`, as parsed statements, against the
installed edgeR's method; on any difference the split is skipped and edgeR's own function runs
whole. See [License](README.md#license).

## Why these five, and not the rest

A function gets a companion only when it passes three tests.

1. **It is slow at cohort scale.**
2. **It has no parallel option of its own.** A scan of the formals of every exported function
   found no `BPPARAM`, `parallel`, `cores`, `mc.cores`, `nthreads` or `workers` argument in
   limma 3.62.2 or edgeR 4.4.2, and neither package imports BiocParallel. DESeq2 1.46.0 takes
   `BPPARAM` in `DESeq`, `results` and `lfcShrink`. In sva 3.54.0, `ComBat` takes `BPPARAM` and
   `ComBat_seq` does not.
3. **It splits by gene rows or sample columns without changing the answer**, bit for bit.

Serial time for each step, on one matrix of 20,000 genes by 300 samples (3 batches, 2 groups,
150 subjects; 17,649 genes after `filterByExpr`), Apple M3, R 4.4.2 with its reference BLAS,
each original at its defaults, best of three for steps under 5 s, made with
[tools/profile_pipelines.R](tools/profile_pipelines.R):

| step | serial | companion | reason |
|---|---:|---|---|
| `limma::voomWithQualityWeights` | 906.5 s | no | sequential: each gene starts from the sample weights the previous gene left |
| `limma::duplicateCorrelation` | 348.9 s | yes | one REML fit per gene |
| `limma::lmFit`, block and correlation, voom weights | 60.0 s | yes | one generalized least squares fit per gene |
| `sva::ComBat_seq` | 32.7 s | yes | quantile match per gene, dispersions per batch |
| `DESeq2::DESeq` | 32.0 s | no | already parallel through `BPPARAM` |
| `edgeR::estimateDisp` | 15.2 s | no | built and removed, not `identical()` |
| `edgeR::glmQLFit` | 4.1 s | no | built and removed, not `identical()` |
| `edgeR::normLibSizes`, TMM | 2.2 s | yes | one factor per sample; 10.3 s on the TCGA matrix |
| `DESeq2::vst` | 0.91 s | no | not timed at cohort scale |
| `limma::lmFit`, voom weights | 0.90 s | yes | one weighted fit per gene; 4.4 s on the TCGA matrix |
| `sva::ComBat`, log-CPM, parametric | 0.55 s | no | already parallel through `BPPARAM`, and cheap |
| `limma::voom` | 0.47 s | no | its `lowess` trend's span is a fraction of the gene count, so only the arithmetic after it splits; a companion measured 0.99x |
| `edgeR::glmQLFTest` | 0.18 s | no | not timed at cohort scale |
| `limma::removeBatchEffect` | 0.11 s | yes | cheap at 3 batches; 4.5 s on the TCGA matrix with 54 plates |
| `DESeq2::results` | 0.06 s | no | cheap |
| `edgeR::filterByExpr` | 0.03 s | no | cheap |
| `limma::eBayes` | 0.013 s | no | cheap, and pools across genes |
| `limma::topTable` | 0.004 s | no | cheap |
| `limma::decideTests` | 0.004 s | no | cheap |
| `edgeR::topTags` | 0.002 s | no | cheap |
| `limma::contrasts.fit` | 0.001 s | no | cheap |

The rest fall into five groups.

**Pools across genes**, so a row split changes the answer: limma's `eBayes`, `squeezeVar`,
`fitFDist`, `treat`, `topTable` and `decideTests` (through `p.adjust`), `voom`'s `lowess` trend,
`arrayWeights` (every method), `normalizeBetweenArrays`, `normalizeQuantiles`, cyclic loess,
`camera`, `roast`, `fry` and `romer`; edgeR's `estimateDisp`, `estimateGLMCommonDisp`,
`estimateGLMTagwiseDisp` with `prior.df > 0` and `glmQLFit` (it squeezes the quasi-likelihood
dispersions); sva's `sva`, `svaseq`, `num.sv` and `fsva`; DESeq2's `estimateDispersionsFit`.

**Not timed at cohort scale**, so no measurement shows that a fork would pay: `vst` (0.91 s) and
`glmQLFTest` (0.18 s) were timed only on the 300-sample matrix above, and `glmLRT`, `exactTest`
and `cpm` not at all. `contrasts.fit`, `topTable`, `decideTests`, `topTags`, `filterByExpr` and
`results` each took under 0.1 s here. `glmQLFTest`, `glmLRT` and `exactTest` also fit through
edgeR's GLM kernels, the same kind of kernel the two removals below rest on.

**Already parallel.** DESeq2's `DESeq`, `results` and `lfcShrink` take `BPPARAM`. With
`fitType = "glmGamPoi"` and `parallel = TRUE`, DESeq2 warns that this is not implemented and
drops the glmGamPoi estimator, so the result differs from a serial glmGamPoi run. `sva::ComBat`
takes `BPPARAM` and splits across batches (default `SerialParam()`). sva's `f.pvalue` is matrix
algebra that a multithreaded BLAS already spreads across cores.

**Built, measured and removed.** Each ran faster, but its output was not `identical()`, and a
small speedup is not worth a different answer. Both rest on an edgeR fitting kernel whose result
for a gene that does not converge cleanly depends on which other genes share its block.

| companion | speedup measured | what changed in the output |
|---|---|---|
| `edgeR::estimateDisp` | 1.4x | 19,999 of 20,000 tagwise dispersions changed, once one library was over-sequenced |
| `edgeR::glmQLFit` | 1.6x | 22 of 18,270 TCGA genes returned a different deviance and iteration count depending on which chunk they landed in, with no flag set |

**Sequential by design.** `voomWithQualityWeights`, with its default `method = "genebygene"`,
updates the sample weights inside its loop over genes: each gene starts from the weights the
previous gene produced, so gene i depends on genes 1 to i-1 and no split reproduces it.

**Using the companions together.** Each companion takes its original's arguments and returns its
original's object, so the next step reads it unchanged. One step needs care.
`voom(counts, design, block = b)` with `correlation` unset passes `correlation = NULL` to
limma's `lmFit`, whose blocked branch then runs a serial `duplicateCorrelation` on the whole
matrix inside `voom`. Run `duplicateCorrelation_parallel()` first and pass its consensus to both:

```r
v   <- voom(dge, design)
dc  <- duplicateCorrelation_parallel(v, design, block = subject)
v   <- voom(dge, design, block = subject, correlation = dc$consensus.correlation)
fit <- lmFit_parallel(v, design, block = subject, correlation = dc$consensus.correlation)
```

`lmFit_parallel()` given `block` and `correlation = NULL` resolves the consensus itself through
`duplicateCorrelation_parallel()`, on the full matrix, so that path is parallel as well. With
`correlation` left out, limma's `lmFit` stops and asks for it, and so does the companion.

## Cross-platform

Every companion returns `identical()` output on macOS, Linux, and Windows, across forked and
socket workers. The rendered reports ran with BLAS pinned to one thread: R's reference BLAS on
macOS and Windows, OpenBLAS on Linux. A multithreaded BLAS has not been tested. The path most
exposed to it is `lmFit_parallel()` with `block`, no probe weights and every value finite, which
splits limma's triangular solve (`backsolve`) by gene; the reports do not reach that path.

| | macOS | Linux | Windows |
|---|---|---|---|
| CPU | Apple M3 | 2 x Intel Xeon @ 2.80 GHz | Intel Core Ultra 9 185H |
| Physical cores | 8 (4 full-speed + 4 efficiency) | 16 | 16 (6 full-speed + 8 efficiency + 2 low-power) |
| Logical cores | 8 | 32 | 22 |
| Can fork a worker | yes | yes | **no** |
| Dispatch | forked children | forked children | separate processes over a socket |
| Backend | `mclapply` | `mclapply` | `future` (`multisession`) |
| Workers swept | 2, 4, 6, 8 | 2, 4, 8, 16 | 2, 4, 6, 8 |

Fork means the operating system clones the running R process, so a worker shares memory with the
parent and starts instantly. Windows cannot fork, so each worker is a fresh R process and the
data has to be sent to it over a socket, which costs time the other two platforms do not pay.

Run time in seconds on the TCGA cohort, 18,270 genes by 1,500 samples in 54 plates, the same
matrix as the README table. Each cell shows the original's time, the fastest companion time, and
the worker count that reached it.

| stage | macOS | Linux | Windows |
|---|---:|---:|---:|
| `sva::ComBat_seq` | 1,619.8s -> **298.6s** (8w) | 3,043.8s -> 325.9s (16w) | 2,573.3s -> 720.1s (6w) |
| `limma::duplicateCorrelation` | 663.7s -> 175.0s (6w) | 469.1s -> **64.9s** (16w) | 784.8s -> 379.3s (4w) |
| `edgeR::normLibSizes` (TMM) | 10.3s -> **1.5s** (8w) | 19.5s -> 4.7s (8w) | 26.5s -> 15.2s (2w) |
| `limma::lmFit` | 4.4s -> **1.5s** (8w) | 11.4s -> 3.4s (16w) | 10.1s -> 9.2s (4w) |
| `limma::removeBatchEffect` | 4.5s -> **1.5s** (6w) | 2.5s -> 1.8s (8w) | 4.4s -> 4.9s (2w) |

On Linux, set `OPENBLAS_NUM_THREADS=1` before starting R. Otherwise every forked worker starts
its own set of math threads and they compete for the same cores. This matters for large matrix
multiplication. It barely matters for `lmFit`, since
limma solves one small problem per gene rather than one large one.

Windows' original run time sits between the M3 and the Xeon (2,573.3s), so its lower speedup is
not because the starting point was already fast, it genuinely finishes slower. With no fork,
every worker is a full separate process, and only 6 of the 16 cores are full-speed, so
ComBat-seq peaks at 3.57x at 6 workers and gets worse past that. `lmFit` and `removeBatchEffect`
on Windows (1.10x and 0.90x) run at close to the original's speed, not faster, by design: both
size gates are set so these two never split on Windows, because sending the data to a fresh
worker there costs more than the split saves. Left to split anyway, `lmFit` measured 0.14x at
21.6M cells and never caught up to the original, which is why the gate exists.

## Tuning internals

**Backends:** `"mclapply"` (forks, default except as below), `"future"`, `"BiocParallel"`,
`"foreach"`, `"serial"`, or any custom `function(idx, f, workers)`. All return identical results.
Forking wins by default because a forked worker reads the parent's matrix directly, with no copy.
A socket worker has to receive its own copy of every chunk, which on limma's `lmFit` row split
measured slower than not parallelizing at all.

**Windows cannot fork**, so backend choice is the whole decision there. `mclapply` runs serially
on Windows and says so once, and so does `BiocParallel`, whose `MulticoreParam` cannot fork there
either. Of the two backends that run workers on Windows, `"foreach"` gets worse as workers are
added: on the simulated ComBat-seq matrix in the Windows report (10,000 genes by 1,000 samples,
10 batches) it fell from 1.65x at 2 workers to 0.84x at 8, because doParallel's cluster form
sends every task over a socket with its closure and data. `"future"` keeps one set of workers
running across every step and reached 3.50x at 6 workers on the same matrix. Use `"future"` on
Windows.

```r
library(future); plan(multisession, workers = 6)   # on Windows the package then picks "future"
```

On Windows the backend is chosen per call: when `future` and `future.apply` are installed and the
active plan has more than one worker, the default is `"future"`; otherwise it is `mclapply`, with
a one-time notice that it is running serially. Everywhere else the default is `mclapply`
whatever plan is set; pass `parallel_backend = "future"` to use a plan there. Setting
`options(combat.backend = "future")`, or any other backend name, replaces this default for every
call on every platform. The package never sets a plan itself; a plan you set is yours to manage.

**Workers** default to `min(8, detectCores() - 2)`. On Windows that is capped again at the
performance-core count, because a socket worker on an efficiency core also costs a serialized
copy; on macOS and Linux it is not.

**Nesting is blocked on every backend.** Calling a companion from inside another parallel call
runs that inner call serially instead of opening a second worker pool underneath the first.
Inside your own fork-based loop (`mclapply`, `future` with a multicore plan) the default backend
also runs serially, while your own socket or callr workers are not marked and dispatch as usual.
Spend the worker budget inside one loop iteration, not spread thin across the whole loop.

**Size gates decide when an input is too small to bother splitting.** Nine options set the
cell or gene count below which a call runs serially instead:

| option | default | counts | what it gates |
|---|---:|---|---|
| `combat.min.cells` | 20,000 | cells | ComBat-seq quantile match; `lmFit_parallel()` with probe weights or a voom `EList`, together with `combat.min.wt.genes` |
| `combat.min.disp.cells` | 30,000 | cells | tagwise dispersion row split, which ComBat-seq reaches only when `covar_mod` makes a batch's design more than one-way |
| `combat.min.glm.cells` | 100,000 | cells | ComBat-seq GLM fit row split |
| `combat.min.batch.cells` | 20,000 | cells, whole filtered matrix | ComBat-seq common and tagwise dispersion across batches |
| `combat.min.ls.cells` | 6e6 | cells | `lmFit_parallel()` without probe weights |
| `combat.min.norm.cells` | 2e5 | cells | `calcNormFactors_parallel()` TMM and TMMwsp column loops |
| `combat.min.order.cells` | 4e6 | cells | `calcNormFactors_parallel()` RLE and upperquartile column loops, and the quantile loop that picks TMM's reference column |
| `combat.min.dupcor.cells` | 5,000 | cells | `duplicateCorrelation_parallel()` |
| `combat.min.wt.genes` | 2,000 | genes, not cells | `lmFit_parallel()` with probe weights or a voom `EList` |

`removeBatchEffect_parallel()` takes whichever `lmFit_parallel()` gate its inner call reaches:
`combat.min.ls.cells` without weights, `combat.min.cells` and `combat.min.wt.genes` when weights
pass through `...`.

An unset gate also depends on whether the backend copies the payload to its workers instead of
inheriting it through a fork. Copying backends are every backend on Windows, `future` without a
multicore plan, `foreach` unless a doParallel or doMC registration drives `mclapply`, `serial`,
and any backend under `options(combat.fork = FALSE)`; a custom backend function counts as
inheriting wherever fork exists. On a copying backend, `lmFit_parallel()`'s two gates,
`combat.min.ls.cells` and its own reading of `combat.min.cells`, close to `Inf`, so neither
lmFit branch splits: over sockets the split never caught up to the original (0.14x at 21.6M
cells, 0.24x at 60M, measured on Windows). `combat.min.norm.cells` rises to 2e6, since TMM pays
over sockets only above about 2 million cells (1.05x at 1.8M, 1.58x at 21.6M), and
`combat.min.order.cells` rises to 4e7. ComBat-seq's gates, including its own reading of
`combat.min.cells`, do not move.

Any gate can be overridden explicitly, and an option you set always wins; the output is
`identical()` either way, since the gate only decides speed, not correctness.
`options(combat.fork = FALSE)` forces every call to run serially regardless of size.
`combat_cluster_stop()` releases any cached worker pool.

**`calcNormFactors_parallel()`** wraps `normLibSizes` on current edgeR, or `calcNormFactors` on
older versions. `normLibSizes` errors on negative counts where the older name returned factors
with a warning instead.
