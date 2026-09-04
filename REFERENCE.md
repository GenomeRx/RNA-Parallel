# Reference

Deep detail split out of the README: internals, full cross-platform numbers, and tuning. Read
this when the README's summary is not enough.

## Glossary

- **companion**: a function such as `lmFit_parallel()` that returns the same output as the
  original function, faster.
- **worker**: one extra R process that does part of the work.
- **chunk**: the block of gene rows (or sample columns) one worker handles at a time.
- **dispatch**: handing chunks out to workers and collecting the results.
- **fork**: on macOS and Linux, the system clones the running R process. The clone shares memory
  with the parent, so no data is copied to start a worker.
- **backend**: the mechanism used to run workers (`mclapply`, `future`, `BiocParallel`, `foreach`,
  or serial).
- **size gate**: an input-size threshold. Below it, the companion runs the original code
  unchanged, since forking a worker costs more than the split saves at that size.
- **cells**: genes times samples. Every size gate counts in this unit.
- **parity**: running at 1.0x, the same speed as the original.

## When a companion is worth reaching for

Measured on an Apple M3, at the default worker count. Ratio is companion speed against the
original function, same input, same output.

| companion | smallest size tested (genes x samples) | size where it turns faster | largest size tested |
|---|---|---:|---:|
| `duplicateCorrelation_parallel()` | 1.68x at 200 x 12 | faster at every size tested | 9.52x at 3,000 x 100 |
| `calcNormFactors_parallel()` | 1.44x at 2,000 x 20 | faster at every size tested | 6.26x at 20,000 x 500 |
| `ComBat_seq_parallel()` | 0.69x at 300 x 20 | about 1,000 genes | 5.37x on the full cohort |
| `lmFit_parallel()`, voom or probe weights | 0.59x at 1,000 x 24 | about 4,000 genes | 2.79x at 60,000 x 48 |
| `lmFit_parallel()`, no probe weights | 1.00x | 6M cells | splits only above the gate |
| `removeBatchEffect_parallel()` | 0.92x at 20,000 x 50 | 6M cells | 1.91x at 20,000 x 500 |

Three groups. `duplicateCorrelation_parallel()` and `calcNormFactors_parallel()` are faster at
every size: `duplicateCorrelation` fits a separate model per gene, so there is always work to
divide, and `calcNormFactors_parallel()` computes the gene ranking once instead of once per
sample, which is faster even with no workers. `ComBat_seq_parallel()` and `lmFit_parallel()` with
probe weights have a floor below which starting a worker costs more than it saves.
`lmFit_parallel()` without probe weights, and `removeBatchEffect_parallel()`, run at the
original's speed until the input is large: limma turns every gene into one `lm.fit()` call that
already takes milliseconds, and `removeBatchEffect` is just that call plus one matrix product.

Using a companion on data that is too small costs one original call plus about 0.3 to 0.6 ms.
That only adds up if it runs thousands of times in a loop, because each call checks its own size
gate and cannot see the loop around it.

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

The engine column reports what actually ran, not what was asked for. A call under a size gate
says `serial`, because `identical()` output at serial speed otherwise looks exactly like a real
parallel run. If a pinned copy of upstream code stands down after an upstream change, the line
says `[match_quantiles stood down]`.

**A live counter by default.** Every parallel call prints a running count on one line, redrawn in
place (`  ComBat-seq 18,270 x 1,500: 42 dispatched`), updated up to 4 times a second.
`combat.timing` alone reports once a call finishes, but ComBat-seq dispatches its slow internal
steps up to `2 * n_batch + 3` times per call, so a cohort of hundreds of batches used to print
nothing between "computing" and the final line, and a slow run looked identical to a stuck one.
The line clears itself before the `combat.timing` summary prints. Set `options(combat.progress =
FALSE)` for silence.

**File progress, for a call that blocks for hours.** The console counter above only updates
between rounds of parallel work. During a round, the main R process is waiting on its workers and
cannot print anything, and the workers cannot print to your console either: a forked worker's
output does not reliably reach an RStudio Server session, and socket workers share no console
with the main process at all.

Set `options(combat.progress.dir = "some/writable/path")` before the call. Each worker appends a
line to its own file at the start and end of every chunk it runs. From a second R session, while
the run is still going:

```r
rnaparallel_progress("some/writable/path")
#> 47 done, 128 started, 0 stalled, 12.4m/chunk, ETA 16:42
```

`done`/`started`/`stalled` count chunks. `stalled` means started but never finished, which is
what a killed worker looks like. Seconds per chunk and the ETA need at least two finished chunks
to mean anything, and read `NA` until then. This is off unless a directory is set. Each write is
one line per chunk, not per gene, so the cost does not show up next to the compute itself even
over hundreds of chunks and hours.

A live progress bar reads the same directory, from that same second session:

```r
rnaparallel_progress("some/writable/path", watch = TRUE)
#> |==================================================| 62%  ComBat-seq 40,609 x 9,493  79/128  ETA 16:42
```

Redraws every `interval` seconds (default 1), on one line. Stops on its own once no new chunk has
started for `stall_after` seconds (default 600), whether that means the run finished or nobody is
writing to `dir` at all. Each poll reads only the bytes added since the last one, not the whole
file again.

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

**`rp_mem_cap()`** runs automatically before every dispatch, reading how much memory is
available and the caller's own memory use. When the requested worker count would need more than
80% of what is available, it lowers the worker count instead and warns with all three numbers:

```
rnaparallel: 16 workers need ~40 GB on top of a 5 GB parent and only 10 GB is
available, which on a machine without swap is a kernel kill, not an R error.
Using 3 instead. Set options(combat.mem.divergence=) if this workload dirties
less, or options(combat.mem.guard=FALSE) to disable.
```

Reads `NA` off Linux, or anywhere the system does not expose this reading: the guard leaves
worker count unchanged, since it never blocks what it cannot measure. `combat.mem.divergence`
(default 1) is the fraction of the parent's memory each worker is expected to modify, and so
copy. The default of 1 assumes the worst case, a full copy. Lower it for a workload known to
modify less, for example a per-column trimmed mean. `combat.mem.guard = FALSE` disables the whole
check.

If the main R process is killed, its workers keep running and keep their share of memory until
the machine restarts, since they were reparented to the operating system's init process, not shut
down. A normal kill signal does not stop them either, because R installs its own handler and the
worker is blocked mid-computation. Every worker now checks whether its recorded master process is
still alive and exits on its own if it is not.

**`rnaparallel_set_mem_limit()`** is a second, independent safeguard against the same failure.
`rp_mem_cap()` lowers the worker count based on a live reading before a fork; this instead sets
`R_MAX_VSIZE`, the most memory R will let itself use, for every future session. R checks this
limit every time it asks for memory, so hitting it fails with a normal R error instead of a
kernel kill:

```r
rnaparallel_set_mem_limit(dry_run = TRUE)   # see the computed value, write nothing
rnaparallel_set_mem_limit()                 # write it to ~/.Renviron
```

Reads total RAM, halves it by default (`fraction =`), and rounds to the nearest of
8/16/32/64/128/256/512/1024 GB. Writes to `~/.Renviron` (or a `path` passed explicitly, which is
also how to test this without touching a real file). The write only takes effect on the next R
session, since `.Renviron` is read once at startup. On Windows, `path.expand("~")` resolves via
`USERPROFILE`, not the `HOME` environment variable, which matters if you were planning to
redirect it by setting `HOME`.

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
[run_example.R](inst/examples/run_example.R) runs the same check at a size where the speedup
shows.

## How it works

`ComBat_seq_parallel()` runs the original `sva::ComBat_seq` code without changing it. It only
swaps in parallel versions of six internal functions that the original calls. Shares below are
of serial run time, measured on 10,000 genes by 500 samples.

| internal function | share of serial run time | work divided by | why the output is unchanged |
|---|---|---|---|
| `match_quantiles` | 66.5% | gene rows | each output cell only reads its own gene |
| `estimateGLMTagwiseDisp` | 14.0% | gene rows | safe only when `prior.df` is 0; any other value runs the original edgeR code on the whole matrix |
| `estimateGLMCommonDisp` | 13.1% | batches | each batch's estimate does not depend on the other batches, and uses no random numbers |
| `glmFit`, `glmFit.default` | remainder | gene rows | offset, dispersion, weights, and start values all slice with the rows |
| `monte_carlo_int_NB` | small | not split | each draw depends on the previous batch's state |

`calcNormFactors_parallel()` splits by sample column. `lmFit_parallel()` and
`duplicateCorrelation_parallel()` split by gene row.

`removeBatchEffect_parallel()` does no splitting of its own. The original `removeBatchEffect` is
one `lmFit` call plus some arithmetic. The companion swaps that one `lmFit` call for
`lmFit_parallel()` and gets all of its speedup from there. This is worth doing because the
original, serial `lmFit` inside `removeBatchEffect` gets slower than linearly as sample count
grows: measured 8.5 seconds at 948 samples, 44.2 seconds at 3,000 samples, and at 9,493 samples
the plain original did not finish inside a 10-minute benchmark cutoff at all.

Chunks are handed to workers in interleaved order (row 1 to worker A, row 2 to worker B, and so
on), so a matrix sorted by expression level does not leave one worker idle while the others
finish. If a worker dies, a chunk comes back twice, or a result comes back short, the run stops
rather than silently returning a wrong answer.

**Some functions stay serial, so the output stays identical.** Some combine information across
every gene at once: `eBayes`, `squeezeVar`, `fitFDist`, `arrayWeights`, `normalizeBetweenArrays`,
`normalizeQuantiles`. Others scale with the total gene count rather than with any one gene:
`voom`'s `lowess` span, `p.adjust`, `topTable`, `decideTests`. Splitting either kind by rows would
change the answer, so they run unmodified. `contrasts.fit` takes 0.001 seconds, less than the
cost of starting one worker, so there is nothing to gain there either. `lmFit(method = "robust")`
and `ndups >= 2` reshape rows in a way that would break a row split and error instead, which is
why `duplicateCorrelation_parallel()` requires `block` to be set explicitly.

**Two companions were built, measured, and removed.** Each ran faster, but its output was not
`identical()`, and a small speedup is not worth a different answer.

| companion | speedup measured | what changed in the output |
|---|---|---|
| `estimateDisp` | 1.4x | 19,999 of 20,000 tagwise dispersions changed, once one library was over-sequenced |
| `glmQLFit` | 1.6x | 22 of 18,270 TCGA genes stopped converging, depending on which chunk a gene landed in |

**Splitting a matrix by rows changes what limma sees in three ways.** Each one returns a wrong
answer with no warning. The package checks all three before it ever splits, reproduced against
limma 3.62.2:

| trap | wrong answer | guard |
|---|---|---|
| `asMatrixWeights` reads the row count of `block`, not of the full matrix, to decide what kind of weights it is looking at | per-array weights get read as per-gene weights when `block`'s height happens to match the sample count | weights are expanded against the whole matrix once, before any split |
| limma picks one of two fitting methods based on whether every weight in the matrix is finite | a chunk that happens to be all-finite can take a faster method the full matrix would not have qualified for | the finite check runs once on the whole matrix, before splitting |
| `stats::lm.fit` silently turns a one-column response into a plain vector | a one-gene chunk swaps a matrix function (`colMeans`) for a vector function (`mean`), changing the shape of the result | every chunk keeps at least two genes |

One internal function is kept as an exact copy of the original source and stands down if that
source changes upstream: `sva::match_quantiles`, which is 66.5% of ComBat-seq's serial time. See
[License](README.md#license).

## Cross-platform

Every companion returns `identical()` output on macOS, Linux, and Windows, across forked and
socket workers, and across threaded and single-threaded BLAS.

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

Run time in seconds on the full cohort (40,609 genes by 9,493 samples). Each cell shows the
original's time, the fastest companion time, and the worker count that reached it.

| stage | macOS | Linux | Windows |
|---|---:|---:|---:|
| ComBat-seq | 1,619.8s -> **298.6s** (8w) | 3,043.8s -> 325.9s (16w) | 2,573.3s -> 720.1s (6w) |
| duplicateCorrelation | 663.7s -> 175.0s (6w) | 469.1s -> **64.9s** (16w) | 784.8s -> 379.3s (4w) |
| calcNormFactors (TMM) | 10.3s -> **1.5s** (8w) | 19.5s -> 4.7s (8w) | 26.5s -> 15.2s (2w) |
| lmFit | 4.4s -> **1.5s** (8w) | 11.4s -> 3.4s (16w) | 10.1s -> 9.2s (4w) |
| removeBatchEffect | 4.5s -> **1.5s** (6w) | 2.5s -> 1.8s (8w) | 4.4s -> 4.9s (2w) |

On Linux, set `OPENBLAS_NUM_THREADS=1` before starting R. Otherwise every forked worker starts
its own set of math threads and they compete for the same cores. This matters for large matrix
multiplication (3.6x difference measured, 61 to 216 GFLOPS). It barely matters for `lmFit`, since
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

**Backends:** `"mclapply"` (forks, default), `"future"`, `"BiocParallel"`, `"foreach"`,
`"serial"`, or any custom `function(idx, f, workers)`. All return identical results. Forking wins
by default because a forked worker reads the parent's matrix directly, with no copy. A socket
worker has to receive its own copy of every chunk, which measured slower than not parallelizing
at all.

**Windows cannot fork**, so backend choice is the whole decision there. `mclapply` runs serially
on Windows and says so once; `BiocParallel` substitutes a serial version of itself. Of the two
backends that actually run in parallel on Windows: `"foreach"` gets worse as more workers are
added (1.18x at 2 workers, 0.28x at 16), because it rebuilds its worker pool every time the chunk
count changes, and ComBat-seq changes that count on every step. `"future"` keeps one worker pool
running across every step and scales the way you would expect. Use `"future"` on Windows.

```r
library(future); plan(multisession, workers = 6)   # package picks "future" for you
```

Backend is chosen per call: an active `future` plan selects `"future"`; no plan leaves
`mclapply` (with a one-time notice that it is running serially, on Windows). The package never
sets a plan itself; a plan you set is yours to manage.

**Nesting is blocked on every backend.** Calling a companion from inside another parallel call
runs that inner call serially instead of opening a second worker pool underneath the first.
Your own R loops are unaffected; only nested companion calls are blocked. Spend the worker budget
inside one loop iteration, not spread thin across the whole loop. Inverting a 15-cohort screen
this way measured 245s to 81.6s.

**Size gates decide when an input is too small to bother splitting.** Nine options set the
cell or gene count below which a call runs serially instead:

| option | default | counts | companion it gates |
|---|---:|---|---|
| `combat.min.cells` | 20,000 | cells | ComBat-seq (general) |
| `combat.min.disp.cells` | 30,000 | cells | ComBat-seq dispersion step |
| `combat.min.glm.cells` | 100,000 | cells | ComBat-seq GLM fit step |
| `combat.min.ls.cells` | 6e6 | cells | `lmFit_parallel()`, `removeBatchEffect_parallel()` |
| `combat.min.norm.cells` | 2e5 | cells | `calcNormFactors_parallel()` |
| `combat.min.order.cells` | 4e6 | cells | ComBat-seq ordering step |
| `combat.min.dupcor.cells` | 5,000 | cells | `duplicateCorrelation_parallel()` |
| `combat.min.batch.cells` | 20,000 | cells | ComBat-seq batch step |
| `combat.min.wt.genes` | 2,000 | genes, not cells | `lmFit_parallel()` with probe weights |

On Windows, two of these move from their default, in opposite directions. `combat.min.ls.cells`
is raised so `lmFit_parallel()` and `removeBatchEffect_parallel()` never split on Windows, since
`lmFit` never caught up to the original there (0.14x at 21.6M cells, 0.24x at 60M).
`combat.min.norm.cells` rises to 2e6, since TMM normalization only pays off on Windows above
roughly 2 million cells (1.05x at 1.8M cells, 1.58x at 21.6M).

Any gate can be overridden explicitly; the output is `identical()` either way, since the gate
only decides speed, not correctness. `options(combat.fork = FALSE)` forces every call to run
serially regardless of size. `combat_cluster_stop()` releases any cached worker pool.

**`calcNormFactors_parallel()`** wraps `normLibSizes` on current edgeR, or `calcNormFactors` on
older versions. `normLibSizes` errors on negative counts where the older name returned factors
with a warning instead.
