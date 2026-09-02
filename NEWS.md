# rnaparallel 0.5.0

Windows is a supported platform, the socket backends work, and the companions no longer pay for
work they throw away.

- **`match_quantiles_parallel()`: bind+reorder fused into one scatter, one allocation not two.**
  `do.call(rbind, parts)` (allocate the whole output) then `m[ord, , drop = FALSE]` (allocate
  it again to undo the interleaving) is now `m[idx[[k]], ] <- parts[[k]]` per chunk into a
  single preallocated matrix -- the permute copy is gone entirely, and it was always paid on
  the real parallel path since interleaved chunks are never already sorted. Safe specifically
  here because `match_quantiles_rows()` already strips dimnames from every chunk, so there is
  no rowname parity to reconstruct across the scatter; `glmFit_rows_parallel`'s equivalent
  `bind()` carries real gene names and is left on the rbind+permute path pending that separate
  verification. `identical()` to the pre-fusion path across the existing 204-case suite,
  including the randomised match_quantiles equivalence tests.

- **Verbose/progress format made consistent across all five companions.** A second Fable
  review, this time of `verbose.R` and its five call sites, found the timing line, watch-mode
  bar, and default label wording had each drifted independently:
  - One shared label width (`RP_LABEL_WIDTH`, 40 chars) now backs the timing line, the
    watch-mode bar's stage field, and the console tick. Previously the timing line truncated
    at 34 and the watch bar at 28, both narrower than `duplicateCorrelation`'s own default
    label (36 chars), so those two surfaces silently cut it mid-number
    (`...18,270 x 1,50`) while the untruncated console tick did not.
  - `calcNormFactors_parallel()`'s default label now reads `calcNormFactors TMM 12,000 x 700`
    instead of bare `TMM 12,000 x 700`, matching the other four companions' "own function name
    first" convention (`lmFit ...`, `ComBat-seq ...`, `duplicateCorrelation ...`,
    `removeBatchEffect ...`) so a shared `combat.progress.dir` across a multi-stage script
    groups TSV rows under a name that actually maps back to the call that produced them.
  - Fixed two `@param label` doc examples (`lmFit_parallel()`, `duplicateCorrelation_parallel()`)
    that had been copy-pasted from ComBat-seq's own doc and showed the wrong default string.
  - `removeBatchEffect_parallel()` now carries the same "timing and quieting are on.exit
    hooks" comment the other four companions already had at their own `rp_step_begin()` call.

  Deferred to a follow-up (real findings, larger/riskier changes): unifying worker
  warning/message replay behaviour (currently only `duplicateCorrelation_parallel()` collects
  and replays distinct child conditions once; the same PSOCK/multisession swallowing risk
  applies to the other four); and having a nested re-entry skip re-running `rp_mem_cap()` a
  second time within one user-facing call, which can currently print a second degrade warning
  and/or have the timing line's engine column understate the actual worker count used by an
  inner dispatch.

- **Fable review pass: dispatch overhead trimmed on every call, socket serialization halved on
  the edgeR RLE path, one fewer allocation on ComBat-seq's quantile matcher.** Three concrete
  fixes from a full-package memory review, each verified with the existing 204-case identical()
  suite (0 new failures) and a full `R CMD check` (0 ERROR/WARNING):
  - `combat_parallel_lapply()` used to build its chunk-tagging machinery (`idx_tagged`,
    `f_tagged`, `untag`) before checking whether the call was even going to dispatch. Every
    early-return path (nested-worker re-entry, `combat.fork = FALSE`, a degenerate single
    chunk/worker) paid a full duplicate of `idx` plus a closure allocation for tagging it then
    threw away unused. Moved the tagging build to right before the two branches that actually
    consume it.
  - `rp_apply_shim()` (the edgeR RLE column-parallel path) stored the matrix twice per socket
    task: once directly, once inside `FUN`'s own captured frame (`.calcFactorRLE`'s `data`).
    R's serializer dedups repeated objects within one environment chain, not across two, so
    each task wrote the whole matrix twice and each worker unserialized two copies. Now reads
    `data` back out of `FUN`'s own environment instead of storing a second reference, halving
    per-task payload on the exact path the largest inputs take (falls back to the old
    double-store if a future edgeR release renames the internal variable).
  - `match_quantiles_rows()` allocated a fresh logical `NA` matrix and then promoted it via
    `out[] <- counts_sub`, paying for a type that was always going to become `counts_sub`'s own
    storage mode. Now copies `counts_sub` directly and strips dimnames, one allocation instead
    of two, same output.

- **`rp_mem_cap()` refuses to fork past what the machine can hold.** Forking N workers off a
  large parent is copy-on-write, so it does not cost N times the parent, but a row-split fit
  writes a real fraction of it, and a machine without swap does not return an allocation error
  when that fraction runs out: the kernel SIGKILLs the process, with no condition to catch, no
  traceback, and `mclapply` reporting nothing. Measured on a 40,609 x 9,493 matrix, 125 GB, no
  swap: 16 workers off a 50 GB parent died, 8 off 100 GB died, 4 off 23 GB died at 111 GB, 2 off
  23 GB survived. Reads `MemAvailable` and the caller's own RSS from `/proc` on Linux, or via
  the `ps` package (`ps_system_memory()`/`ps_memory_info()`) on Windows and macOS, before every
  dispatch, and degrades the worker count instead, warning with the three numbers so a run that
  cannot proceed at full concurrency says why instead of vanishing. NA when neither source is
  available means proceed unchanged; the guard never blocks what it cannot measure.
  `combat.mem.divergence` (default 1: assume a worker can dirty the whole parent) is the
  fraction of the parent each worker is assumed to dirty, workload-dependent so it is an
  option; `combat.mem.guard = FALSE` disables the whole check. The default was raised from
  an earlier 0.25 after checking it against the PR's own measured numbers: 4 workers off a
  23 GB parent that actually grew to 111 GB (a real divergence of about 0.96) computed only
  23 GB needed at 0.25 and would have proceeded unwarned into the same kill. Lower it
  explicitly for a workload known to dirty less, e.g. a per-column trimmed mean.

- **`rnaparallel_set_mem_limit()`** is a second, independent net against the same SIGKILL:
  `R_MAX_VSIZE` is R's own vector-heap ceiling, checked on every allocation, that turns an
  overshoot into a catchable `cannot allocate vector of size X` error rather than leaving the
  kernel to silently kill the process. `rp_mem_cap()` above degrades the worker count based
  on a live reading before a fork; this sets a fixed session-wide ceiling instead, in
  `~/.Renviron` (or a `path` you pass), by reading total RAM, halving it, and rounding to the
  nearest of 8/16/32/64/128/256/512/1024 GB. `dry_run = TRUE` shows the computed value without
  writing anything; the write only takes effect on the NEXT R session, since `.Renviron` is
  read once at startup.

- **A fork whose master died now exits on its own.** Reparented to init, an orphaned fork kept
  running and holding its share of the matrix for as long as the machine stayed up; two runs
  left 111 GB and 116 GB stranded that way, immune to `SIGTERM` because R installs a handler
  and the worker is blocked mid-computation. Every worker now checks `Sys.getppid() == 1` and
  exits if its master is gone. `Sys.getppid()` is not a base R function on every build (it does
  not exist at all on the R 4.6.1 UCRT Windows build this package is tested on), so the check
  now goes through `rp_getppid()`, which falls back to `ps::ps_ppid()` when installed and
  returns `NA` (skip the check) otherwise, rather than crashing every single dispatch on a
  build that lacks it.

- **`rnaparallel_progress(dir, watch = TRUE)` now renders a live `|====------|` bar**, matching
  `data.table::fread()`'s own style: `|==================================================| 62%
  ComBat-seq 40,609 x 9,493 79/128 ETA 16:42`. Same mechanism as before, watched from a SEPARATE
  process/terminal while the running session is blocked inside its parallel call.
  Redraws now TAIL each worker's file rather than re-reading it whole every poll:
  verified on a real run, 32 of 34 polls against a growing file skipped or only read the
  new bytes, the other 2 being the first sighting of each worker's file. Matters most on a
  long `watch = TRUE` session against a long-running dispatch, where the read cost used to
  grow with elapsed time even though each poll only cares about the handful of lines
  written since the last one. Stall detection now tracks a "done" row the same as a "start"
  row, not `started` alone: a run whose chunks all dispatched early (`preschedule = TRUE`) but
  whose "done" rows are still trickling in used to hit the stall exit mid-run with a misleading
  "no new chunks" message, even though chunks were actively finishing. A stall also now returns
  the last real progress summary instead of `NULL`, matching the function's own documented
  return contract, and includes the stalled-chunk count in the stall message when nonzero.
  `interval` and `stall_after` are validated (must be a single positive number) rather than
  reaching `Sys.sleep()` as a raw error or busy-polling on a bad value.

- **Peak RSS in the `combat.timing` line**, from `VmHWM` on Linux or `ps::ps_memory_info()`'s
  peak working set on Windows/macOS, so the number that decides whether a
  fit survives is visible in the log a run already produces: `pooled ComBat-seq    mclapply
  x16    30.3h    peak 51 GB`.

- **`combat.progress` is now on by default.** No option needs to be set: every parallel call
  ticks a live "N dispatched" line, overwritten in place, so you can always tell a slow run from
  a stuck one. `combat.timing` alone only reports after a call finishes, and ComBat-seq alone
  dispatches its hot paths up to `2 * n_batch + 3` times per call, so a cohort with hundreds of
  batches used to show nothing between "computing" and the final summary. Set
  `options(combat.progress = FALSE)` to go back to silent.

- **`options(combat.progress.dir = ...)` plus `rnaparallel_progress()`** cover the gap the
  console tick above cannot: once the master calls into the parallel backend it blocks until
  every chunk returns, sometimes for hours, and nothing can print from a blocked process. Set
  a directory and each worker appends a start/done line per chunk to its own file as it works;
  call `rnaparallel_progress(dir)` from a SEPARATE session to see chunks done, a mean
  seconds-per-chunk, and an ETA, while the run is still going. Built from a real pooled run
  (40,609 genes, 9,493 specimens, 464 batches, 16 workers) that ran 14h53m and then 13h more
  with no output at all, on either the console or the workers' own stdout, because forked and
  PSOCK workers cannot reliably reach an RStudio Server console. A file both sides can read
  survives all four backends the package supports. Off unless a directory is set.

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
