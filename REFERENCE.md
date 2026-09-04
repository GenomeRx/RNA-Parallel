# Reference

Deep-dive material split out of the README: internals, full cross-platform data, and tuning
detail. Read this when the README's summary isn't enough.

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
`sva::match_quantiles`, 66.5% of ComBat-seq's serial time. See [License](README.md#license).

## Cross-platform

Every arm on macOS, Linux, and Windows returns `identical()` output: fork and socket dispatch,
threaded and single-threaded BLAS.

| | macOS | Linux | Windows |
|---|---|---|---|
| CPU | Apple M3 | 2 x Intel Xeon @ 2.80 GHz | Intel Core Ultra 9 185H |
| Physical cores | 8 (4P+4E) | 16 | 16 (6P+8E+2 low-power) |
| Logical | 8 | 32 | 22 |
| NUMA nodes | 1 | 2 | 1 |
| BLAS | Accelerate | OpenBLAS 0.3.8 | reference, pinned to 1 thread |
| `fork()` | yes | yes | **no** |
| Dispatch | forked children | forked children | PSOCK processes |
| Backend | `mclapply` | `mclapply` | `future` (`multisession`) |
| Workers swept | 2, 4, 6, 8 | 2, 4, 8, 16 | 2, 4, 6, 8 |

Absolute wall clock at cohort scale (orig → best, worker count):

| stage | macOS | Linux | Windows |
|---|---:|---:|---:|
| ComBat-seq | 1,619.8s → **298.6s** (8w) | 3,043.8s → 325.9s (16w) | 2,573.3s → 720.1s (6w) |
| duplicateCorrelation | 663.7s → 175.0s (6w) | 469.1s → **64.9s** (16w) | 784.8s → 379.3s (4w) |
| calcNormFactors (TMM) | 10.3s → **1.5s** (8w) | 19.5s → 4.7s (8w) | 26.5s → 15.2s (2w) |
| lmFit | 4.4s → **1.5s** (8w) | 11.4s → 3.4s (16w) | 10.1s → 9.2s (4w) |
| removeBatchEffect | 4.5s → **1.5s** (6w) | 2.5s → 1.8s (8w) | 4.4s → 4.9s (2w) |

On Linux, pin `OPENBLAS_NUM_THREADS` before R starts, or every forked worker opens its own thread
pool (unpinning moves DGEMM 3.6x, 61→216 GFLOPS; `lmFit` barely changes since limma solves a
small QR per gene).

Windows' baseline sits between the M3 and the Xeon (2,573.3s), so its low ratio isn't a
short-baseline effect, it's just slower to finish. No `fork()` means a worker is a full copied
process, and only 6 of 16 cores are performance cores, so ComBat-seq peaks at 3.57x at 6 workers
and turns over past the P-core count. `lmFit`/`removeBatchEffect` on Windows (1.10x/0.90x) are
parity, not speedups: both size gates close where the payload gets copied, so the companion runs
the original plus wrapper overhead. Left to split, `lmFit` measured 0.14x at 21.6M cells and
never reached parity, which is why the gate exists.

## When a companion is worth reaching for

Measured on the M3, at the default worker count. Ratio is companion against original.

| companion | smallest measured | crossover | large |
|---|---|---:|---:|
| `duplicateCorrelation_parallel()` | 1.68x at 200 x 12 | pays at every size measured | 9.52x at 3,000 x 100 |
| `calcNormFactors_parallel()` | 1.44x at 2,000 x 20 | pays at every size measured | 6.26x at 20,000 x 500 |
| `ComBat_seq_parallel()` | 0.69x at 300 x 20 | about 1,000 genes | 5.37x on the cohort |
| `lmFit_parallel()`, voom or probe weights | 0.59x at 1,000 x 24 | about 4,000 genes | 2.79x at 60,000 x 48 |
| `lmFit_parallel()`, no probe weights | 1.00x | 6M cells | splits only above the gate |
| `removeBatchEffect_parallel()` | 0.92x at 20,000 x 50 | 6M cells | 1.91x at 20,000 x 500 |

Three groups: **`duplicateCorrelation`/`calcNormFactors` pay unconditionally** (one gene = one
REML fit; TMM's `rank` hoist wins with zero workers). **`ComBat-seq` and weighted `lmFit` have a
floor** below which fork overhead beats the savings. **`lmFit` without probe weights, and
`removeBatchEffect`, are parity until the input is large**: limma vectorises every gene into one
`lm.fit` (milliseconds), and `removeBatchEffect` is just that call plus a BLAS product.

Choosing wrong costs little: under-gate is one original call plus ~0.3-0.6ms overhead, real only
in a loop over thousands of small units, since gates decide per-call and can't see the loop.

## Tuning internals

**Backends:** `"mclapply"` (forks, default), `"future"`, `"BiocParallel"`, `"foreach"`,
`"serial"`, or any `function(idx, f, workers)`, all return identical results. Forking wins by
default because a forked worker reads the matrix copy-on-write; socket backends re-serialise per
chunk and measured slower than not parallelising.

**Windows can't fork**, so backend choice is the whole decision. `mclapply` runs serially (says
so once); `BiocParallel` substitutes a serial param. Of the two that actually parallelize:
`"foreach"` gets *worse* with more workers (1.18x→0.28x at 2→16w, cached cluster rebuilds on
every width change and ComBat-seq alternates widths per dispatch); `"future"` holds one
`multisession` pool across all of them and scales normally.

```r
library(future); plan(multisession, workers = 6)   # package picks "future" for you
```

Backend resolves per call: an active plan selects `"future"`; no plan leaves `mclapply` (one-time
serial notice). The package never sets a plan itself; a caller's plan is theirs.

**Nesting is blocked on every backend.** A dispatch already inside one of this package's workers
runs serially (`mc.allow.recursive = FALSE` on fork; equivalent guard over PSOCK, which used to
spawn workers² processes without it). Your own loop is unaffected. Spend the worker budget
*inside* a loop body, not across it. Inverting one 15-cohort screen measured 245s → 81.6s.

**Size gates** are why small inputs show no speedup. Nine options set the cell/gene count below
which a call runs serially: `combat.min.cells` (20,000), `combat.min.disp.cells` (30,000),
`combat.min.glm.cells` (100,000), `combat.min.ls.cells` (6e6), `combat.min.norm.cells` (2e5),
`combat.min.order.cells` (4e6), `combat.min.dupcor.cells` (5,000), `combat.min.batch.cells`
(20,000), `combat.min.wt.genes` (2,000, counted in genes not cells).

Two move without fork, in opposite directions: `lmFit` never reaches parity on Windows (0.14x at
21.6M cells, 0.24x at 60M) so `combat.min.ls.cells` **closes** there (covers
`removeBatchEffect_parallel()` too). TMM pays once big enough (1.05x at 1.8M cells → 1.58x at
21.6M) so `combat.min.norm.cells` **rises** to 2e6. Override either explicitly; output is
`identical()` regardless. `options(combat.fork = FALSE)` forces serial everywhere;
`combat_cluster_stop()` releases cached clusters.

**`calcNormFactors_parallel()`** wraps `normLibSizes` on current edgeR, `calcNormFactors` on
older ones. `normLibSizes` errors on negative counts where the old name returned NaN-warned
factors.

## Seeing what is running

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

The engine column reports what ran, not what was asked for: a call under a size gate says
`serial`, since `identical()` output at serial pace otherwise looks exactly like a fork. If a
pinned original excerpt stands down after an upstream change, the line says
`[match_quantiles stood down]`.

**Progress ticks by default.** Every parallel call prints a live count on one line, overwritten
in place with a carriage return (`  ComBat-seq 18,270 x 1,500: 42 dispatched`), throttled to
4/sec. `combat.timing` alone reports once a call finishes; ComBat-seq dispatches its hot paths up
to `2 * n_batch + 3` times per call, so a cohort of hundreds of batches used to show nothing on
screen between "computing" and the final line, and a slow run looked identical to a stuck one.
The line clears itself before the `combat.timing` summary prints. Set
`options(combat.progress = FALSE)` for silence.

**File progress, for a call that blocks for hours.** The console tick above only fires between
dispatches in the master process; once the master calls into `mclapply`/`future`/`BiocParallel`/
`foreach` it blocks synchronously until every chunk returns, and nothing in the master can print
during that block. Workers cannot fix this by writing to the master's console either: a forked
child's stdout is not reliably multiplexed back to an RStudio Server session, and PSOCK/
BiocParallel workers share no console with the master at all.

Set `options(combat.progress.dir = "some/writable/path")` before the call. Each worker appends
a TSV line, one file per worker PID, at the start and end of every chunk it runs. From a
SEPARATE R session, while the run is still going:

```r
rnaparallel_progress("some/writable/path")
#> 47 done, 128 started, 0 stalled, 12.4m/chunk, ETA 16:42
```

`done`/`started`/`stalled` count chunks (`stalled` = started, never finished, which is what a
killed worker looks like). The seconds-per-chunk and ETA need at least two finished chunks to
mean anything and read `NA` until then. Off unless a directory is set; every write is one line
per chunk, not per gene, so the cost is immaterial next to the compute itself even at hundreds
of chunks over hours.

**A live bar**, this package's own format, from the same directory, from that same SEPARATE
session:

```r
rnaparallel_progress("some/writable/path", watch = TRUE)
#> |==================================================| 62%  ComBat-seq 40,609 x 9,493  79/128  ETA 16:42
```

Redraws every `interval` seconds (default 1) on one line, overwritten in place, same mechanism
as the console tick above but drawn from the WATCHING process rather than the one blocked
inside the parallel call, since only the watcher is free to keep redrawing. Stops on its own
once `started` stops growing for `stall_after` seconds (default 600), whether that means the
run finished or nobody is writing to `dir` at all. Each poll only reads the bytes appended
since the last one, not the whole file again: verified against a real run, 32 of 34 polls
either skipped a file that had not changed or only read its new bytes.

`rnaparallel_stale()` returns TRUE if the package was reinstalled under a running session. The
fix is to restart R.

## Memory

Forking is copy-on-write, so N workers do not cost N times the parent process. They cost
whatever each one writes to, and for a row-split fit over a large matrix that is a real
fraction of it. When the total exceeds what the machine has, the kernel does not hand R an
allocation error: on a box without swap it SIGKILLs the process outright, with no condition to
catch, no traceback, and `mclapply` reporting nothing. `future` says only that a future was
interrupted. From outside, R vanishes.

Measured on a 40,609 x 9,493 matrix, 125 GB, no swap:

| workers | parent RSS | outcome |
|---|---|---|
| 16 | 50 GB | killed |
| 8 | 100.8 GB | killed, 0 GB free at the fork |
| 4 | 23 GB | killed, 23 -> 111 GB in 30s |
| 2 | 23 GB | survived, 102 GB peak |

**`rp_mem_cap()`** runs automatically before every dispatch, reading `MemAvailable` and the
caller's own RSS from `/proc`. When the requested worker count would need more than 80% of
what is available, it degrades to a smaller count instead and warns with all three numbers:

```
rnaparallel: 16 workers need ~40 GB on top of a 5 GB parent and only 10 GB is
available, which on a machine without swap is a kernel kill, not an R error.
Using 3 instead. Set options(combat.mem.divergence=) if this workload dirties
less, or options(combat.mem.guard=FALSE) to disable.
```

`NA` off Linux (or anywhere `/proc` is missing) means proceed unchanged; the guard never blocks
what it cannot measure. `combat.mem.divergence` (default 1, assume a worker can dirty the whole
parent) is the fraction of the parent each worker is assumed to dirty, workload-dependent, not
a constant; lower it explicitly for a workload known to dirty less, e.g. a per-column trimmed
mean. `combat.mem.guard = FALSE` disables the whole check.

A fork whose master was killed is reparented to init and keeps running, holding its share of
the matrix for as long as the machine stays up; `SIGTERM` does not clear it, because R installs
a handler and the worker is blocked mid-computation. Every worker now checks
`Sys.getppid() == 1` and exits on its own if its master is gone.

**`rnaparallel_set_mem_limit()`** is a second, independent net against the same failure.
`rp_mem_cap()` degrades the worker count based on a live reading before a fork; this instead
sets `R_MAX_VSIZE`, R's own vector-heap ceiling checked on every allocation, to a fixed value
for every future session:

```r
rnaparallel_set_mem_limit(dry_run = TRUE)   # see the computed value, write nothing
rnaparallel_set_mem_limit()                 # write it to ~/.Renviron
```

Reads total RAM, halves it by default (`fraction =`), and rounds to the nearest of
8/16/32/64/128/256/512/1024 GB. Writes to `~/.Renviron` (or a `path` passed explicitly, which
is also how to test this without touching a real file); the write only takes effect on the
NEXT R session, since `.Renviron` is read once at startup. On Windows, `path.expand("~")`
resolves via `USERPROFILE`, not the `HOME` environment variable, which matters if you were
planning to redirect it by setting `HOME`.

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
[run_example.R](inst/examples/run_example.R) runs the same check at a size where speedup shows.
