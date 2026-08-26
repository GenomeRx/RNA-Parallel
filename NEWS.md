# rnaparallel 0.4.9

Diagnostic fix; the companions are unchanged from 0.4.8.

- `rnaparallel_stale()` now compares the loaded namespace version against the installed one, so it
  answers even from a session whose namespace predates the function. Call it as
  `tryCatch(rnaparallel_stale(), error = function(e) TRUE)` — the error is the answer.
- The roxygen sync test ran `roxygenise()` in-process, replacing the namespace mid-suite. It runs
  in a subprocess now.

# rnaparallel 0.4.8

- **New: `removeBatchEffect_parallel()`.** Rebinds the single `lmFit` call inside
  `limma::removeBatchEffect` and runs the vendor unchanged. 1.70x / 2.83x / 3.37x at 2 / 4 / 8
  workers on TCGA, `identical()` throughout.
- A backend that falls back to serial is no longer reported as parallel in the timing line.
  `mclapply` on Windows and a `future` plan resolving in one process were both counted as forks.
- A failed dispatch says when the cause is environmental: N of N chunks failing with the same
  error now says so, and names a reinstall when `rnaparallel_stale()` is `TRUE`. Failures carry an
  `rnaparallel_dispatch_error` class with `stage`, `n_chunks`, `n_errored` and `n_died`, so
  callers can branch on structure instead of matching message text.
- `?workers` documents the cost of exceeding your performance core count; `?chunks` notes that
  `chunks = workers` is redundant. The README covers nesting a companion inside a forking loop.

# rnaparallel 0.4.7

- **New: `rnaparallel_stale()`.** Reinstalling the package under a running session corrupts its
  lazy-load database, which surfaces inside forked workers as what looks like a data fault. This
  reports it; the fix is to restart R.
- **New: progress reporting.** `options(combat.timing = TRUE)` prints one elapsed line per call,
  `combat.quiet` silences the vendor's own output, `combat.timing.min` hides fast steps, and every
  companion takes `label`. The engine column reports what *ran*: a call under a size gate says
  `serial`, because a companion returning `identical()` output at serial pace otherwise looks
  exactly like one that forked.

# rnaparallel 0.4.5

- **Behaviour change: `workers` defaults to `min(8, detectCores() - 2)`,** not a flat 4. Not a
  safety ceiling — six workers alongside a second forking R session has kernel-panicked a 24 GB
  machine, and nothing here can see that session.
- **Behaviour change: `calcNormFactors_parallel()` wraps `normLibSizes`** on current edgeR. Same
  function, but a matrix with a negative cell now **errors** where it previously returned
  NaN-warned factors; the NA-count and `lib.size` message texts also changed.
- Killing a worker no longer orphans the rest, and could previously segfault. Reached only via
  `parallel_backend = "foreach"`; the default `mclapply` was never affected.
- A `foreach` backend you registered yourself is no longer throttled to `workers`, and a caller's
  `doRNG` registration can no longer move the master random stream.
- Cached-cluster reuse is faster, and a machine with no `ps` on `PATH` no longer fails outright.
- `duplicateCorrelation_parallel()` holds no statmod source, removing the only GPL-derived code
  from an MIT package, and preserves limma's outcome for zero-row input.

# rnaparallel 0.4.3

- **Breaking: `combat.min.dupcor.genes` is now `combat.min.dupcor.cells`,** counted in
  gene-by-array cells, default 5000. The old name is no longer read, so a script setting it gets
  the default instead of the value it asked for.
- `duplicateCorrelation_parallel()` forks on that cells axis; the gene-count gate had been sending
  real work serial.
- Lower memory: a parallel `glmFit` worker returns the six fields the parent reads rather than a
  whole `DGEGLM`, about 540 MB a dispatch on a cohort of this size.
- `lmFit_parallel()`'s array-weight check reads block-leading rows rather than the whole matrix;
  that path had measured 0.98x, slower than the function it wraps.

# rnaparallel 0.4.2

- `data.frame` counts, which `sva::ComBat_seq` accepts, crashed the quantile-match gate. A
  non-matrix input now goes to the vendor whole.
- On sparse input the per-batch tagwise dispatch checked results against the unfiltered row count,
  so the dominant stage was computed in parallel, discarded, and recomputed serially with no
  signal. It now reads the vendor's own filtered matrix.
- Garbage `combat.fork` and `combat.min.*` values refuse loudly instead of silently changing
  behaviour.

# rnaparallel 0.4.1

- **Correctness: `ComBat_seq_parallel()` returned output that depended on `chunks`** when the
  design was batch-only and some fits did not converge. edgeR's one-group kernel is not a pure
  function of the gene it fits in that case. Both `glmFit` and the tagwise dispersion now detect
  that layout and call the vendor whole.
- Fork thresholds retuned against measured crossovers, several of which had made the companion
  slower than the vendor on small inputs. `glmFit` gained its own `combat.min.glm.cells`.
- Dispatch overhead cut: no `BiocParallel` param rebuild and no core-count process per call.

# rnaparallel 0.3.0

Initial public release. Runs `sva::ComBat_seq` unmodified with its hot paths rebound in a child of
the backend's own environment, returning output `identical()` to a serial run rather than merely
close to it. 5.01x at eight workers on TCGA, 18,270 genes by 1,500 tumours across 54 plates.

`estimateGLMCommonDisp` is dispatched across batches rather than gene rows, because it optimises
over a sum across all genes and a row split would change floating-point accumulation order. The
Monte Carlo integration stays serial, since its draws depend on the previous batch's state.
