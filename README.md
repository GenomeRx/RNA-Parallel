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

TCGA, 18,270 genes by 1,500 tumours. Every arm `identical()` to the original on all three
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

Nothing is reimplemented: the original function runs, called with hot paths rebound in a child of
its own environment. `identical()` is asserted, not a tolerance.

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
defaults it to `NULL`. See [REFERENCE.md](REFERENCE.md#how-it-works) for why.

**Check it yourself:** [REFERENCE.md](REFERENCE.md#full-self-check) has the full comparison
against all five originals, no download required. At cohort scale, the same checks run as
[rendered reports](https://genomerx.github.io/RNA-Parallel/) on all three platforms, sourced from
[inst/examples/](inst/examples/). [tests/](tests/testthat) covers every argument path, chunk
layout, backend, and dispatch count: 400+ assertions.

**See what's running:** every parallel call ticks a live "N dispatched" line by default,
overwritten in place, so a long ComBat-seq run against hundreds of batches never sits silent.
`options(combat.timing = TRUE)` adds the real engine per call once it finishes, `serial` vs
`mclapply x6`, etc. For a call that blocks inside the parallel backend for hours, set
`options(combat.progress.dir = "some/path")` and call `rnaparallel_progress(dir, watch = TRUE)`
from a SEPARATE session for a live `|====------|` bar, `data.table::fread()` style, with chunks
done and an ETA. See [REFERENCE.md](REFERENCE.md#seeing-what-is-running).

**Memory:** forking a large matrix can exceed a machine's RAM even when the parent alone fits,
and on a box without swap the kernel SIGKILLs the process with no R error at all. `rp_mem_cap()`
degrades the worker count before that happens, using a live reading, on by default; set
`options(combat.mem.guard = FALSE)` to disable it. `rnaparallel_set_mem_limit()` is a second,
independent net: it sets `R_MAX_VSIZE`, R's own allocation ceiling, to half the machine's RAM so
an overshoot becomes a catchable error instead of a silent kill. See
[REFERENCE.md](REFERENCE.md#memory).

## Tuning

| knob | default | change it when |
|---|---|---|
| `workers` | `min(8, detectCores() - 2)`, capped at performance cores without `fork()` | rarely; going past your performance-core count can be slower |
| `chunks` | `workers` | only to cut peak memory per worker |
| `parallel_backend` | `"mclapply"` | you cannot fork, or a cluster is already running |

Full backend, nesting, and size-gate detail in [REFERENCE.md](REFERENCE.md#tuning-internals).

## License

MIT for this companion, copyright GenomeRx 2026, in [LICENSE](LICENSE). The original packages are
called at run time from your own installation and none of them is redistributed here: sva is
Artistic-2.0, limma and edgeR are GPL (>= 2). One exception, marked in source:
`R/helper_seq_parallel.R` carries a row-vectorised transcription of `sva::match_quantiles` so the
companion can detect an upstream change and stand down (derived from Artistic-2.0 code by Zhang,
Parmigiani, Johnson).

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
Package: <https://bioconductor.org/packages/release/bioc/html/edgeR.html>.

**This companion.** Nguyen N (2026). *rnaparallel: Parallel Companions for RNA-Seq Tools.*
<https://github.com/GenomeRx/RNA-Parallel>.
