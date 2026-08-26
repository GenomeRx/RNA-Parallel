#!/usr/bin/env bash
# Render the Linux verification report for rnaparallel.
#
# OpenBLAS reads its thread count when the library is loaded, so OPENBLAS_NUM_THREADS
# has to be in the environment BEFORE R starts -- calling Sys.setenv() inside the R
# session has no effect (measured on this machine: 0.97x, i.e. none). Without this,
# every forked worker opens its own BLAS thread pool and the box is oversubscribed.
set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1          # harmless if MKL is absent
export VECLIB_MAXIMUM_THREADS=1   # harmless on Linux; keeps the script portable

# pandoc ships with RStudio Server here and is not on PATH.
if [ -z "${RSTUDIO_PANDOC:-}" ]; then
  for d in /usr/lib/rstudio-server/bin/quarto/bin/tools/x86_64 \
           /usr/lib/rstudio-server/bin/quarto/bin/tools \
           /usr/lib/rstudio/bin/quarto/bin/tools/x86_64; do
    [ -x "$d/pandoc" ] && export RSTUDIO_PANDOC="$d" && break
  done
fi
echo "RSTUDIO_PANDOC=${RSTUDIO_PANDOC:-<not found>}"
echo "OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"

IN="${1:-RNA_Parallel_linux.Rmd}"
OUT="${2:-${IN%.Rmd}.html}"

exec Rscript -e "rmarkdown::render('${IN}', output_file='${OUT}', quiet=FALSE)"
