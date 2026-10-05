# RNA-Parallel

[![version](https://img.shields.io/badge/version-0.5.0-blue)](NEWS.md)
[![R](https://img.shields.io/badge/R-%E2%89%A5%204.1-blue)](DESCRIPTION)
[![license](https://img.shields.io/badge/license-MIT-blue)](LICENSE.md)

Parallel companions for five functions in sva, edgeR and limma. Each one calls the original
function unmodified and returns output `identical()` to it, bit for bit. Same arguments, same
defaults (exceptions under [Use](#use)), same result, faster at cohort scale.

Rendered analysis: [macOS](https://namstacks.github.io/RNA-Parallel/) ·
[Linux](https://namstacks.github.io/RNA-Parallel/linux.html) ·
[Windows](https://namstacks.github.io/RNA-Parallel/windows.html)

## Speedup

TCGA, 18,270 genes by 1,500 tumors. Every arm `identical()` to the original on all three
platforms. macOS: M3, 4P+4E. Linux: 2x Xeon, 16 cores. Windows: Ultra 9 185H, 6P+10E.

| companion | runs | macOS | Linux | Windows |
|---|---|---:|---:|---:|
| `ComBat_seq_parallel()` | `sva::ComBat_seq` | 5.43x @ 8w | **9.34x @ 16w** | 3.57x @ 6w |
| `calcNormFactors_parallel()` | `edgeR::normLibSizes` | **6.78x @ 8w** | 4.12x @ 8w | 1.74x @ 2w |
| `lmFit_parallel()` | `limma::lmFit` | 3.02x @ 8w | **3.37x @ 16w** | 1.10x @ 4w |
| `duplicateCorrelation_parallel()` | `limma::duplicateCorrelation` | 3.79x @ 6w | **7.22x @ 16w** | 2.07x @ 4w |
| `removeBatchEffect_parallel()` | `limma::removeBatchEffect` | **2.97x @ 6w** | 1.39x @ 8w | 0.90x @ 2w |

Bold = fastest platform per row. Windows has no `fork()` (a worker is a full copied process, and
only 6 of 16 cores are performance cores), which caps its scaling; `lmFit`/`removeBatchEffect`
there are parity, not speedups. Full breakdown and every wall-clock number in
[REFERENCE.md](REFERENCE.md#cross-platform).

Nothing here is reimplemented, apart from one helper transcribed from sva and three default
statements transcribed from edgeR (see [License](#license)). The original function runs, with
only its slow internal steps swapped for parallel versions in a child of its own
environment. `identical()` is checked directly, not a tolerance.

## Why these five

RNA-Parallel is a companion, not a replacement. It plugs into five steps of a limma, edgeR or
sva ComBat-seq pipeline, each chosen because all three hold:

1. **It is slow at cohort scale.** Most steps finish in seconds; a few take minutes.
2. **It has no parallel option of its own.** No exported limma or edgeR function takes a
   `BPPARAM`, `parallel` or worker argument, and neither does `sva::ComBat_seq`. DESeq2's
   `DESeq()`, `results()` and `lfcShrink()` and sva's `ComBat()` already take `BPPARAM`.
3. **It splits without changing the answer.** Steps that pool every gene into one estimate
   (eBayes priors, dispersion trends, voom's lowess curve, quantile normalization) do not.

Serial time, simulated 20,000 genes by 300 samples, Apple M3 (`tools/profile_pipelines.R`):

| step | serial | companion | why |
|---|---:|---|---|
| `limma::voomWithQualityWeights` | 906.5 s | no | each gene's update feeds the sample weights the next gene starts from, so the loop is sequential and no split reproduces it |
| `limma::duplicateCorrelation` | 348.9 s | yes | one REML fit per gene |
| `limma::lmFit`, block + correlation | 60.0 s | yes | one Cholesky per gene under voom weights |
| `sva::ComBat_seq` | 32.7 s | yes | quantile matching splits by gene rows; common and tagwise dispersion split by batch; `glmFit` splits by rows only when `group` or `covar_mod` is given |
| `DESeq2::DESeq` | 32.0 s | no | already parallel through `BPPARAM` |
| `edgeR::estimateDisp`, `edgeR::glmQLFit` | 15.2 s, 4.1 s | no | built and removed: output was not `identical()` |
| `edgeR::normLibSizes` (TMM) | 2.2 s | yes | one factor per sample, each against one reference column chosen once from all samples |
| `limma::lmFit` (voom weights), `limma::removeBatchEffect` | 0.9 s, 0.1 s | yes | one weighted fit per gene; one `lmFit` with the batch columns added |
| `DESeq2::vst`, `edgeR::glmQLFTest` | 0.9 s, 0.2 s | no | not timed at cohort scale |
| `limma::voom` | 0.5 s | no | its `lowess` trend's span is a fraction of the gene count, so only the arithmetic after it splits; a companion measured 0.99x |
| `eBayes`, `topTable`, `filterByExpr` and other quick steps | under 0.1 s | no | nothing to gain |

TMM, `lmFit` and `removeBatchEffect` are cheap here and grow with the cohort: on the Speedup
matrix (18,270 genes by 1,500 tumors) they take 10.3 s, 4.4 s and 4.5 s serial on the M3
([REFERENCE.md](REFERENCE.md#cross-platform)). `vst` and `glmQLFTest` were not timed at that
scale, so nothing yet shows a fork would pay for them. Each companion takes its original's
arguments and returns its original's object, so `voom`, `eBayes`, `topTable`, `glmQLFit` and
DESeq2 read the output unchanged. `voom(counts, design, block = b)` with `correlation` unset
runs a serial `duplicateCorrelation` inside `voom`, so run `duplicateCorrelation_parallel()`
first and pass its `consensus.correlation` to `voom` and to `lmFit_parallel()`;
`lmFit_parallel()` given `block` and `correlation = NULL` resolves it through
`duplicateCorrelation_parallel()` itself. Every function considered:
[REFERENCE.md](REFERENCE.md#why-these-five-and-not-the-rest).

## Install

The package needs `edgeR` and nothing else. Every other package below is optional, and
only for the companion you actually call.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("edgeR")

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("NamStacks/RNA-Parallel",
                        dependencies = c("Depends", "Imports"),
                        upgrade = "never")
```

Keep `upgrade = "never"`. Without it the installer offers to rebuild every outdated
package in your library, most of which have nothing to do with this one.

That is four packages on a clean library: `edgeR`, plus `limma`, `locfit` and `statmod`,
which `edgeR` and `limma` require themselves. `lattice` and the rest ship with R.

Only one companion asks for anything beyond that:

| Companion | Also install |
|---|---|
| `ComBat_seq_parallel()` | `sva` |
| the other four | nothing |

Backends work the same way. `"mclapply"` and `"serial"` need nothing; `"future"` needs `future`
and `future.apply`, `"BiocParallel"` needs `BiocParallel`, and `"foreach"` needs `foreach` and
`doParallel`. On macOS and Windows, install `ps` so the memory guard can read memory; Linux reads
`/proc` and needs nothing. Without a reading the guard leaves the worker count as it is. Where base
R has no `Sys.getppid()`, `ps` also lets a forked worker read its parent's process id and exit
once its master has died; without it, orphaned workers keep running.

## Use

Each companion takes its original's arguments in the same order with the same defaults, and adds
`workers`, `chunks`, `parallel_backend`, `backend` and `label`. All five in one pass:

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

Interface differences: `duplicateCorrelation_parallel` requires `block`, whereas the original
defaults it to `NULL`. `lmFit_parallel` and `removeBatchEffect_parallel` refuse
`method = "robust"` and `ndups >= 2`, which limma accepts; call limma directly for those.
`removeBatchEffect_parallel` declares `design = NULL`, current limma's default; left unset it is
not forwarded, so the installed limma's own default applies. `calcNormFactors_parallel` refuses any class other than `DGEList` that edgeR has its own method
for, such as a `SummarizedExperiment`; pass its count matrix instead. See
[REFERENCE.md](REFERENCE.md#how-it-works) for why.

**Check it yourself:** [REFERENCE.md](REFERENCE.md#full-self-check) has the full comparison
against all five originals, no download required. At cohort scale, the same checks run as
[rendered reports](https://namstacks.github.io/RNA-Parallel/) on all three platforms, sourced from
[inst/examples/](inst/examples/). [tests/](tests/testthat) covers every argument path, chunk
layout, backend, and dispatch count: 400+ assertions.

**See what's running:** on macOS and Linux every parallel dispatch draws a live `|====------|`
bar by default, with chunks done and an ETA, overwritten in place. It draws only in an interactive
session or on a terminal. A call that runs serially, every call on Windows, and every call whose
output goes to a log tick a running "N dispatched" count instead.
`options(combat.progress = FALSE)` silences both. `options(combat.timing = TRUE)` adds the real
engine per call once it finishes, `serial` vs `mclapply x6`, etc., and `label = "..."` names the
call in the progress and timing lines. The bar is drawn by a forked child, whose output may not
reach an RStudio Server console. To follow a logged run, or one on RStudio Server or Windows, set
`options(combat.progress.dir = "some/path")` and call
`rnaparallel_progress("some/path", watch = TRUE)` from a second session. See
[REFERENCE.md](REFERENCE.md#seeing-what-is-running).

**Memory:** forking a large matrix can exceed a machine's RAM even when the parent alone fits,
and on a box without swap the kernel SIGKILLs the process with no R error at all. The memory
guard degrades the worker count before that happens, using a live reading, once per companion
call, on by default; set `options(combat.mem.guard = FALSE)` to disable it.
`rnaparallel_set_mem_limit()` sets `R_MAX_VSIZE`, R's own allocation ceiling, to half the
machine's RAM rounded down to a whole tier (8, 16, 32 GB and up; an exact Mb value below 8 GB), so
one process's overshoot becomes a catchable error instead of a silent kill. The limit is per
process: each forked worker inherits it, so it does not bound the total across workers. On macOS
it reads total RAM through `ps`. See [REFERENCE.md](REFERENCE.md#memory).

## Tuning

| knob | default | change it when |
|---|---|---|
| `workers` | `min(8, detectCores() - 2)`, capped at performance cores on Windows | rarely; going past your performance-core count can be slower |
| `chunks` | `workers` | only to cut peak memory per worker |
| `parallel_backend` | `getOption("combat.backend")`, else `"mclapply"`; on Windows `"future"` when a multi-worker future plan is active | you cannot fork, or a cluster is already running |

Full backend, nesting, and size-gate detail in [REFERENCE.md](REFERENCE.md#tuning-internals).

## License

MIT for this companion, copyright GenomeRx 2026, in [LICENSE](LICENSE). The original packages are
called at run time from your own installation and none of them is redistributed here: sva is
Artistic-2.0, limma and edgeR are GPL (>= 2). Two exceptions, both marked in source in
`R/helper_seq_parallel.R`, and both used only while the installed original still matches them:
a row-vectorized transcription of `sva::match_quantiles` (derived from Artistic-2.0 code by
Zhang, Parmigiani, Johnson), and `.tagwise_defaults_pinned`, three one-line default formulas
(`offset`, `span`, `AveLogCPM`) taken from edgeR 4.4.2's `estimateGLMTagwiseDisp.default`,
which the tagwise row split computes once for the whole matrix. On any upstream change the
companion detects the difference and calls the original instead.

## Citation

Cite the method you used and this companion. `citation("rnaparallel")` prints every entry below.

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
Chen Y, Chen L, Lun ATL, Baldoni PL, Smyth GK (2025). edgeR v4: powerful differential analysis
of sequencing data with expanded functionality and improved support for small counts and larger
datasets. *Nucleic Acids Research* 53(2), gkaf018.
doi:[10.1093/nar/gkaf018](https://doi.org/10.1093/nar/gkaf018).
Package: <https://bioconductor.org/packages/release/bioc/html/edgeR.html>.

**This companion.** Nguyen N (2026). *rnaparallel: Exact Parallel Companions for sva ComBat-Seq, edgeR and limma.*
<https://github.com/NamStacks/RNA-Parallel>.
