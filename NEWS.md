# rnaparallel 0.5.0

Windows is a supported platform, the socket backends work, and the companions no longer pay for
work they throw away.

- **Three companions could not run on the `future` backend at all.** A closure is serialised with
  its defining environment, and every dispatched job in the package was built in a frame holding
  far more than its body reads. Under `plan(multisession)` on a 54.9 MiB matrix,
  `lmFit_parallel()` exported 664.29 MiB of globals and `future` refused it outright;
  `calcNormFactors_parallel()` exported 595.97 MiB. Both were hard errors where the original
  returns a value, on the backend the new adaptive default selects. Every job is now built
  against an environment holding only what it reads: measured payload per task, lmFit
  13,209,456 B to 809,288 B, TMM 16,624,160 B to 1,987,841 B, `duplicateCorrelation` 5,961,556 B
  to 853,219 B. Two unforced promises were part of it — `combat_parallel_lapply()` never forced
  the function it dispatches, so each task carried the caller's frame and, through its own
  unforced promises, the original's frames and the entry point's raw inputs.

- **The size gates were keyed to the operating system when they are about the payload.** They are
  break-evens for a worker that INHERITS the matrix, and Windows is a sufficient condition for
  copying rather than a necessary one. It is not even about `fork()`: `foreach` builds a FORK
  cluster on Unix and is still slow, because doParallel's cluster form serialises every task
  whatever its nodes were made with. Measured on 300,000 x 24 at four workers, every arm
  `identical()`: `limma::lmFit` 0.150 s, the companion on `mclapply` 0.111 s, on `foreach`
  0.628 s. That last is 0.24x, and `removeBatchEffect_parallel()` inherited it. The gates now
  resolve against the backend that will actually run, and `foreach` is asked rather than assumed —
  `doParallelMC` drives `mclapply` and copies nothing, `doParallelSNOW` does not.

- **Both `lmFit` branches were reading `combat.min.ls.cells`.** Raising the weightless gate also
  switched off the voom branch's split, silently. The test suite could not see it: every
  `combat.min.*` is set to 0 for the suite, which returns from the same first line for both
  branches.

- **`lmFit_parallel()` was slower than `limma::lmFit` on array-weighted input below its own
  gate**, because it scanned the whole matrix for branch stability before reaching the gate that
  sends the call serial anyway. 20,000 x 24 went from 0.73x to 1.09x, 60,000 x 48 from 0.65x to
  1.02x. That scan also answered one TRUE or FALSE with five full-size allocations, now a min/max:
  200,000 x 48 measured 40.85 ms to 27.04 ms.

- **ComBat-seq stopped paying for scans it discards.** The tagwise separability gate ran two
  full-slice scans before a design test that vetoes it on every batch of every run; testing the
  design first is 1.66 ms to 0.033 ms on an 18,270 x 28 slice and stops charging each worker
  ~6 MB of transient allocation per batch. `match_quantiles`, the dominant stage, is vectorised
  over the whole slice instead of per row: 0.508 s to 0.370 s, checked against
  `sva::match_quantiles` on 600 randomised cases. Interleaved chunking uses a strided sequence
  rather than `split()`, and a single chunk no longer copies every bound field twice.

- **Two more size gates, for the two regimes that were still slower than the original.** limma's
  weighted branch is an interpreted per-gene loop whose cost barely moves with array count, so
  genes amortise the fork and cells were the wrong unit: at a fixed 1,000 genes the split behaved
  the same at 8, 24 and 48 arrays while the cell count crossed the gate that decided it. New
  `combat.min.wt.genes` (2,000) takes that branch from 0.57x-0.75x below it to parity, and 1.52x
  above. ComBat-seq's two across-batch dispatches carried no floor at all, on the argument that a
  whole-matrix estimate per batch is always worth a fork; true at cohort scale and false at 300
  genes, where they measured 0.69x. New `combat.min.batch.cells` (20,000) takes that to 1.20x —
  faster than the original while dispatching nothing, because the vectorised quantile match is a
  serial win.

- **A lean environment must be parented at the package, not at the global environment.** A
  dispatched closure is rebuilt against an environment holding only what its body reads, and that
  environment's parent decides how the body's remaining calls resolve. Parented at `globalenv()`,
  a user binding named `vapply` reached them: measured, it moved `calcNormFactors_parallel`
  8.59e-06 off `edgeR::normLibSizes` with no error and no warning, on the serial backend, on every
  platform, while leaving the original untouched. Six sites, now parented where those calls resolved
  before the closures were leaned.

- The quantile-match gate refused non-finite inputs but not a negative `old_mu`, an `old_phi`
  shorter than the matrix, or a mismatched `dim`. The original errors on all three; the vectorised
  form returned a plausible half-matched matrix. Unreachable from ComBat-seq, whose fitted values
  are non-negative, but the gate's contract is its own. Its finiteness tests are now reductions
  rather than full-matrix logicals, and `combat_row_order` inverts a permutation by scatter
  instead of sorting it.

- **New: each companion's help page says when it is worth reaching for.** Two of the five lose to
  their original below a floor and nothing said so. `duplicateCorrelation` and `calcNormFactors` pay
  unconditionally; ComBat-seq has a floor near a thousand genes; the weighted `lmFit` branch one
  near four thousand; the unweighted branch and `removeBatchEffect` are parity until the matrix is
  large. Under a gate a companion is one original call plus 0.3 to 0.6 ms.

- **The package stopped warning about its own default.** `min(8, detectCores() - 2)` is 6 on an
  M3 with 4 performance cores, so every fresh session opened by saying that might be slower than
  a number it had declined to pick — contradicting the 5.37x at eight workers this package
  publishes for that chip. There were two performance-core routines disagreeing (8 against 4);
  there is one now, and the message fires only where the payload is copied, which is where it is
  true.

- A cached cluster is no longer torn down when a NARROWER pool is asked for. ComBat-seq alternates
  widths within one run, and each change rebuilt the pool: measured on PSOCK, 20 calls alternating
  3/6 workers took 6299 ms against 22 ms at a fixed width. `calcNormFactors_parallel()` was the one
  companion still defaulting to a literal `"mclapply"` while its own documentation said otherwise,
  and the Windows serial notice recommended `foreach`, which this release measures at 0.24x.

Windows is a supported platform, and the dispatch layer no longer nests.

- **Nested dispatch is blocked on every backend, not just fork.** `ComBat_seq_parallel()`
  dispatches the tagwise loop across batches and ships the original closure, whose environment
  still carries the rebound `estimateGLMTagwiseDisp`; inside the worker that symbol dispatched
  again over gene rows. The guard against this was `mc.allow.recursive = FALSE`, an argument to
  `parallel::mclapply`, so it covered the fork branch alone. Windows runs `mclapply` serially and
  reaches its workers only through `foreach`/PSOCK, which had no equivalent — so the one platform
  that needed the guard was the one platform without it. Measured there: `workers = 2L` produced
  two outer workers and four nested ones, which extrapolates to 272 processes at a 16-worker arm,
  each a fresh R process with edgeR and limma loaded. `combat_parallel_lapply()` now marks the
  worker process it dispatches into and takes the serial path when it finds itself already inside
  one. A caller's own parallel loop is unaffected, since nothing marks their workers.

- **The default backend is chosen adaptively where there is no fork().** `mclapply` cannot fork on
  Windows, so it is correct, serial, and says so once. `future` is the only backend that runs real
  workers there without the caller registering a cluster, and it measured 2.92x on the cohort
  against `foreach`'s 0.28x — but it needs a plan, and this package will not set one. Defaulting to
  it unconditionally would hand most callers a warning per dispatch and no speedup. So it is taken
  exactly when a plan is already active: `plan(multisession, workers = 6)` is now the whole step,
  with no `options(combat.backend=)` to discover. Resolved per call, and unchanged on every
  platform that can fork.

- **The default worker count is capped at PERFORMANCE cores where there is no fork().** On a hybrid
  CPU those are not the core count: an Ultra 185H reports 16 physical cores, 6 of which are
  performance cores and 10 efficiency. Performance cores are read from the topology rather than a
  original table — only they carry SMT, so logical minus physical gives 6 — degrading to the physical
  count on uniform machines. The cap is deliberately not applied where fork() exists: macOS
  measures 5.37x at eight workers on a chip with four performance cores, because a forked worker on
  an efficiency core still adds throughput, while a socket worker also costs a serialised copy. On
  Windows the default moves from 8 to 6, which is where the measured ComBat-seq curve peaks.

- **The TMM column split gate rises without fork(), where the least-squares gate closes.** Same
  problem, opposite answers, which is why both were measured rather than assumed. TMM does pay once
  it is large enough — 1.05x at 1.8M cells against 1.58x at 21.6M — so `combat.min.norm.cells`
  moves to 2e6, an order of magnitude above the fork break-even, instead of closing. Treating it
  like `lmFit` would have discarded a real 1.6x.

- **The least-squares row split no longer runs where it cannot pay.** `combat.min.ls.cells`
  defaults to 6e6 cells, which is a *fork* break-even: a forked worker starts almost free and
  reads the matrix through copy-on-write. Without `fork()` every chunk is serialised to its own
  R process, and `lm.fit` is cheap enough per cell that the transfer is never repaid. Measured
  over PSOCK with the gate forced open, 1,200 samples: 0.14x at 21.6M cells and 0.24x at 60M,
  improving with size and never approaching parity, and worse with every worker added. On the
  TCGA cohort it measured 0.50x, and `removeBatchEffect_parallel()` inherits it because the
  original rebinds a single `lmFit` call -- so two companions were slower than the functions they
  wrap. No threshold rescues that, so on a platform without fork the fast branch now stays
  whole: both measure at parity instead of 0.50x and 0.38x. An explicit `combat.min.ls.cells`
  still reaches the split, and output is `identical()` either way -- this decides who computes,
  never what is computed.

- **New: the Windows verification report.** `inst/examples/RNA_Parallel_windows.Rmd` carries the
  same sections in the same order as the macOS and Linux reports, and adds a backend sweep,
  because on a platform with no `fork()` the backend decides whether anything runs in parallel
  at all. `render_windows.ps1` pins the BLAS thread count before R starts and finds pandoc from
  a standalone install as well as from RStudio.

- A test that reaches `mclapply` through a custom executor now skips on Windows rather than
  failing there, matching the guard its sibling tests already carried.

# rnaparallel 0.4.9

Diagnostic fix; the companions are unchanged from 0.4.8.

- `rnaparallel_stale()` now compares the loaded namespace version against the installed one, so it
  answers even from a session whose namespace predates the function. Call it as
  `tryCatch(rnaparallel_stale(), error = function(e) TRUE)` — the error is the answer.
- The roxygen sync test ran `roxygenise()` in-process, replacing the namespace mid-suite. It runs
  in a subprocess now.

# rnaparallel 0.4.8

- **New: `removeBatchEffect_parallel()`.** Rebinds the single `lmFit` call inside
  `limma::removeBatchEffect` and runs the original unchanged. 1.70x / 2.83x / 3.37x at 2 / 4 / 8
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
  `combat.quiet` silences the original's own output, `combat.timing.min` hides fast steps, and every
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
  non-matrix input now goes to the original whole.
- On sparse input the per-batch tagwise dispatch checked results against the unfiltered row count,
  so the dominant stage was computed in parallel, discarded, and recomputed serially with no
  signal. It now reads the original's own filtered matrix.
- Garbage `combat.fork` and `combat.min.*` values refuse loudly instead of silently changing
  behaviour.

# rnaparallel 0.4.1

- **Correctness: `ComBat_seq_parallel()` returned output that depended on `chunks`** when the
  design was batch-only and some fits did not converge. edgeR's one-group kernel is not a pure
  function of the gene it fits in that case. Both `glmFit` and the tagwise dispersion now detect
  that layout and call the original whole.
- Fork thresholds retuned against measured crossovers, several of which had made the companion
  slower than the original on small inputs. `glmFit` gained its own `combat.min.glm.cells`.
- Dispatch overhead cut: no `BiocParallel` param rebuild and no core-count process per call.

# rnaparallel 0.3.0

Initial public release. Runs `sva::ComBat_seq` unmodified with its hot paths rebound in a child of
the backend's own environment, returning output `identical()` to a serial run rather than merely
close to it. 5.01x at eight workers on TCGA, 18,270 genes by 1,500 tumours across 54 plates.

`estimateGLMCommonDisp` is dispatched across batches rather than gene rows, because it optimises
over a sum across all genes and a row split would change floating-point accumulation order. The
Monte Carlo integration stays serial, since its draws depend on the previous batch's state.
