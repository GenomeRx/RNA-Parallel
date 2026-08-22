# rnaparallel 0.4.4

`duplicateCorrelation_parallel` no longer holds any statmod source. The old lift transcribed the
leading statements of `mixedModel2Fit` so the design-invariant QR and SVD could be computed once,
and needed roughly a hundred lines of body-text pinning to notice if statmod ever drifted, because
a lift that broke inside limma's `tryCatch(..., error = nafun)` would have returned NA for every
gene rather than an error. The vendor's own `mixedModel2Fit` object now runs unchanged and only
the primitive it repeats is shadowed: `La.svd(QtZ, nu = mq, nv = 0)` is called once per gene with
byte-identical arguments, measured at 300 calls across 300 genes and 47.5% of the vendor's wall
clock. The memo keys on a bitwise `identical()` of the actual arguments, so it needs no gate: a
statmod that stops calling `La.svd` produces a cache that never hits, which is slow rather than
wrong.

What it costs, measured on the same cohort with only this code changed: nothing on the arm this
package benchmarks, and some speed on the arm it does not. The TCGA sweep moves from 1.97x, 3.14x
and 4.34x to 1.96x, 3.18x and 4.45x at 2, 4 and 8 workers, inside the 3.2% the two baseline
readings differ by. That arm passes a `voom` result, and `mixedModel2Fit` scales `X` by each
gene's own weights before it reaches `La.svd`, so the argument is not design-invariant there:
counted, 300 genes give 300 calls and 300 distinct argument sets, a 0% hit rate. The same fact is
why the old lift refused to install on weighted input at all. On unweighted input the argument is
invariant, 300 calls give 1 distinct set, and the memo reaches 5.12x serial where the transcription
reached 8.23x. That is the whole trade: an input almost nobody passes gets slower, the published
numbers do not move, and the only GPL-derived code in an MIT package is gone.

Removing that excerpt removes the only GPL-2 | GPL-3 derived code from an MIT package. The licence
section, README and inline notices follow.

Smaller: the reachability gate covered three of the five symbols the dupcor rebind installs, and
now covers `factor` and `model.matrix` too. `lmFit_parallel`'s documentation listed `EListRaw` as
an accepted class, which limma rejects outright, and omitted `PLMset` and `marrayNorm`, which work.

# rnaparallel 0.4.3

`combat.min.dupcor.genes` is now `combat.min.dupcor.cells`, counted in gene-by-array cells rather
than genes, default 5000 in place of 150. The old name is no longer read, so a script setting it
gets the default instead of the value it asked for. The gene-count gate sent real work serial:
one gene is one REML fit and the fit's cost grows with the arrays too, so a hundred genes across
a thousand arrays measured 4.97x once it was allowed to fork, and the old gate refused it.

`duplicateCorrelation_parallel` now forks on that cells axis. Everything measured at or above
5000 cells gains; below roughly 3000 wins and losses interleave, because serial cost tracks the
number of fits and the per-fit array cost rather than their product, so the gate sits above the
mixed band rather than at a crossover.

Speed and memory, all `identical()` to the vendor. Each parallel `glmFit` worker returns the six
fields the parent reads instead of a whole `DGEGLM`; edgeR also ships back full counts, dispersion
and offset slices that were discarded on arrival, about 540 MB a dispatch on a cohort of this
size. `lmFit_parallel`'s finiteness guard prefilters with `rowSums` before the exact scan, and its
array-weight uniformity check now reads the rows each block leads with rather than the whole
matrix. That path previously measured 0.98x, slower than the function it wraps.

Fixed in the harness, not the package: `tools/exactproof.R` ran at production defaults, where its
fixtures fall under the size gates, so checks covering the parallel path executed the serial
fallback and still reported green. It now zeroes the gates the way the test suite does, and its
chunk loop no longer drops the default arm to `c()` eating a `NULL`. The count went from 70 checks
to 85, all passing. One limma regression test named for the gls punch rule was refused earlier by
the array-weight check and never reached what it tested; its weight moved off a block-leading row.

Documentation: the README said the repository contains no copy of ComBat-seq. It holds two gated
vendor excerpts, the `sva::match_quantiles` body and the head of `statmod::mixedModel2Fit`, each
pinned as text so the companion stands down when upstream changes. The README and LICENSE now
say so, the licence section names both origins, and each block carries its own notice in source. Sub-linear scaling was explained
by workers being stranded on efficiency cores; measured, eight forked children given identical
work finish within 1.08x of each other, so that mechanism is gone and the concurrency cost it should
have named is stated instead.

# rnaparallel 0.4.2

A fresh-eyes review pass over the whole package, every finding reproduced before it was fixed.

Two mattered. `data.frame` counts, which `sva::ComBat_seq` accepts, crashed the quantile-match
gate because `is.finite` has no data.frame method; a non-matrix input now goes to the vendor
whole. And the per-batch tagwise dispatch checked its results against the unfiltered row count,
while ComBat-seq drops genes that are all zero within a batch before that stage, so on any
sparse input the dominant stage was computed in parallel, discarded, and recomputed serially
with no signal. The check now reads the vendor's own filtered matrix. Both carry regression
tests.

Smaller: garbage `combat.fork` and `combat.min.*` values refuse loudly at the entry point
instead of silently changing behaviour; concurrent PSOCK cluster creation retries on a port
collision; the `future` backend's Windows claim in the README now matches the code; documented
gate defaults match the code; the duplicated entry-point validation, call-head walker and
row-subset helpers are each one definition; and the rebind count is documented as six, naming
the per-batch tagwise `lapply`.

# rnaparallel 0.4.1

## Fork thresholds retuned against measured crossovers

Every gate that decides whether a dispatch is worth a fork was re-measured, and several sat above
the point where forking starts to pay, which made the companion slower than the vendor on small
inputs. `combat.min.disp.cells` 5e4 to 3e4, `combat.min.norm.cells` 5e5 to 2e5,
`combat.min.order.cells` 5e6 to 4e6, and `glmFit` gained its own `combat.min.glm.cells` at 1e5
where it had been inheriting a 2e4 threshold that forked fits up to 3.4x slower than running them
whole. The dispatch layer also stopped rebuilding a `BiocParallel` param on every call (212 ms to
60 ms) and stopped spawning a process to count cores on every call.


## Correctness fix in the ComBat-seq companion

`ComBat_seq_parallel` returned output that depended on `chunks` when the design was
batch-only and some genes' fits did not converge. `edgeR`'s one-group kernel is not a pure
function of the gene it is fitting in that case: the same gene's coefficient changes with
which other genes share the matrix. Measured on 240 genes by 12 samples with one library
over-sequenced 1000x, 1, 2 and 3 cells differed at chunks 2, 4 and 8. Both `glmFit` and the
tagwise dispersion now detect the layout and call the vendor whole.

The speed lost to that gate was recovered on a different axis rather than accepted. ComBat-seq
computes the tagwise dispersion once per batch, and batches are independent with no shared state
and no RNG, so those calls are now dispatched across batches instead of across gene rows. That
never reaches the kernel a row split does. It is also where the time is: profiled at 6,000 genes
by 400 samples across 10 batches, the tagwise dispersion was 60% of the corrected run and the
quantile match was 7%.

A batch-only design measures 3.49x at 2 workers, 5.33x at 4 and 5.97x at 8, so the exact
implementation is now well past the inexact one it replaced (3.35x at 8). A design with
covariates through `covar_mod` never entered the gate and is unaffected; the TCGA figure of 4.99x
stood, and the rendered report now measures 5.01x. Every fixture in the suite used benign
counts with near-equal library sizes, which is why 285 tests stayed green over this; the
over-sequenced fixture is now a regression test.

Parallel companions for limma and edgeR, built the same way as the ComBat-seq one: the vendor
function is called unchanged on each block, with at most a symbol rebound in a child of its own
environment, and output is `identical()` to a serial run.

| companion | vendor | axis | speedup, rendered TCGA at 8 workers |
|---|---|---|---:|
| `calcNormFactors_parallel()` | `edgeR::calcNormFactors` | sample columns | 7.64x |
| `duplicateCorrelation_parallel()` | `limma::duplicateCorrelation` | gene rows | 4.29x |
| `lmFit_parallel()` | `limma::lmFit` | gene rows | 4.03x |

`calcNormFactors_parallel()` carries a change that needs no workers at all. TMM evaluates
`rank(logR)` and `rank(absE)` twice each per sample; computing them once is bit-identical and
worth 1.58x by itself. The split is by sample column rather than gene row, because every
normalisation method ranks or takes a median across genes.

## Deliberately not parallelised

`voom` measured 0.99x, since its `lowess` trend takes `f = span` as a fraction of the gene count
and only the post-trend arithmetic splits. `contrasts.fit` is exactly splittable and the vendor
takes 0.001 s. `eBayes`, `fitFDist` and `topTable` pool across all genes by construction. None of
them ships as a `workers` argument that does nothing.

`estimateDisp` was built and then deleted, which is worth recording because it looked like the
best target on the list. `adjustedProfileLik` is 83.6% of it and reads only its own gene's row, so
the split is sound on paper and it measured 1.4x. It was not `identical()` at stock defaults once
one library was heavily over-sequenced: 19,999 of 20,000 tagwise dispersions moved. The mechanism
was not pinned down, and a 1.4x gain is not worth a companion whose output depends on `chunks`.
An exactness suite for this package needs a fixture with a 1000x library-size outlier; benign
`rnbinom` matrices with near-equal library sizes never reach the branch that fails.

`glmQLFit` was built, measured at 1.6x, and deleted after the report's own assert caught it on
real TCGA data: `mglmLevenberg` records a deviance and an iteration count whose values depend on
which genes share the block when a fit does not converge cleanly, 22 of 18,270 genes, with no
flag set on any of them. Every small fixture had passed. `glmQLFTest` is 0.13 s and was not
built.

Every entry point now reaps its own fork children on exit, success or error. A child wedged in a
signal-unsafe state survives `mclapply`'s cleanup, and survivors accumulating across a long
session has crashed a 24 GB machine. Children that predate the call, a user's own `future` or
`mclapply` workers, are snapshotted on entry and spared.

`trend.method = "none"` errors inside `edgeR::estimateDisp` 4.4.2 whenever any gene falls below
`min.row.sum`, because `m0` is built with `nrow(y)` rows while `l0` has `sum(sel)`. That is an
upstream bug, reproduced against the untouched vendor, and is noted here only so nobody attributes
it to this package.

## Three ways a limma row split goes wrong

All three were reproduced against limma 3.62.2 and all three return a wrong answer quietly.

`asMatrixWeights` dispatches on the block's row count and tests its gene branch before its array
branch, so a per-array weight vector is read as per-gene weights by any block whose row count
equals the sample count. Coefficients moved by 0.325. Weights are expanded once against the full
matrix before anything is split.

`NoProbeWts` is an AND-reduction over every cell selecting between two numerically different
algorithms that return different component sets. An all-finite block flips to the fast path while
the whole matrix took the slow one, and 114 of 400 `sigma` differed. The branch is now proved
stable before splitting, and the run falls back to serial when it cannot be.

`stats::lm.fit` drops a one-column response to a vector, so a one-gene block swaps `colMeans` for
`mean` and loses the gene name: 4,080 of 16,000 one-gene blocks were not `identical()`. Chunks now
carry at least two genes.

# rnaparallel 0.3.0

Initial public release. `rnaparallel` collects parallel companions for RNA-seq tools, each one
returning output `identical()` to the function it replaces. ComBat-seq is the first; a limma
companion is in progress.

## What it does

Runs `sva::ComBat_seq` unmodified and parallelises its hot paths. Five symbols are rebound in a
child of the backend's own environment: `glmFit`, `glmFit.default`, `match_quantiles` and
`estimateGLMTagwiseDisp` split by gene row, and `sapply` dispatches the per-batch common
dispersion across batches. Output is `identical()` to a serial run at every stage, not merely
close to it.

`estimateGLMCommonDisp` cannot be split by gene row: it optimises over a sum across all genes, so
splitting changes floating-point accumulation order and can move the argmax in the last place.
Batches are independent of each other and carry no RNG, so dispatching across batches is exact
where dispatching across rows would not be. The Monte Carlo integration stays serial, since its
draws depend on the state the previous batch left.

## Measured

TCGA HNSC, LUAD and LUSC, 18,270 protein-coding genes by 1,500 primary tumours across 54
sequencing plates, on an Apple M3 with 4 performance and 4 efficiency cores:

| implementation | seconds | speedup |
|---|---:|---:|
| `sva::ComBat_seq` | 1656.1 | 1.00x |
| `ComBat_seq_parallel`, 4 workers | 467.7 | 3.54x |
| `ComBat_seq_parallel`, 8 workers | 330.4 | 5.01x |

Simulated counts at 10,000 genes by 1,000 patients across 10 batches, `group` supplied, each
arm in a fresh session: 3.72x at four workers and 4.43x at eight. Every arm returns a matrix
`identical()` to the original.

## Guardrails

`identical()` alone cannot tell a working parallel layer from a dead one: a wrapper that forwards
to `sva::ComBat_seq` satisfies every equivalence test. `combat_backend()` therefore collects the
heads of calls whose head is a bare symbol and refuses to run when a rebind target is missing, and
the suite counts dispatches through the public entry point. `all.names()` is not sufficient for
that check, since it flattens `edgeR::glmFit` to include `glmFit`.

Chunks are interleaved rather than contiguous, because gene order is not random and contiguous
blocks of an expression-sorted matrix leave workers idle. Each chunk carries a tag and is
reassembled by it, so a backend returning results out of order cannot bind genes to the wrong
rows. A dead worker, a duplicate chunk, a chunk of the wrong height, or a result count below the
request halts the run rather than being patched.

`workers` is run as asked, bounded by `detectCores()` and by the two-core limit `R CMD check`
imposes. On a hybrid CPU the efficiency cores contribute a fraction of a performance core, so
asking for more workers than performance cores can be slower on some workloads; the package says
so once per session rather than reducing the number.

Below `combat.min.cells` a row-split dispatch runs serially, since the fork costs more than the
work. That gate selects which process computes, never the result. The common dispersion has no
such gate, each batch being one whole estimate.

Windows cannot fork, so the default `mclapply` backend and `BiocParallel`'s `MulticoreParam` fall
back to serial there and say so once per session. `parallel_backend = "foreach"` uses a PSOCK
cluster and `"future"` a multisession plan, both of which run in parallel on Windows.

230 assertions pass under `devtools::test()`, 222 under `R CMD check`, and `R CMD check` is clean.
