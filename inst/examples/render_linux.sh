#!/usr/bin/env bash
# Render the Linux verification report for rnaparallel.
#
# OpenBLAS reads its thread count when the library is loaded, so OPENBLAS_NUM_THREADS
# has to be in the environment BEFORE R starts -- calling Sys.setenv() inside the R
# session has no effect (measured on this machine: 0.97x, i.e. none). Without this,
# every forked worker opens its own BLAS thread pool and the box is oversubscribed.
#
#   ./render_linux.sh
#   ./render_linux.sh RNA_Parallel_linux.Rmd ../../docs/linux.html
set -euo pipefail
cd "$(dirname "$0")"

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
# A machine with no RStudio Server at all is the other common case, and there pandoc is
# usually just on PATH. Without this the guard below would refuse a machine that can in fact
# render, which is worse than having no guard.
if [ -z "${RSTUDIO_PANDOC:-}" ]; then
  p="$(command -v pandoc || true)"
  [ -n "$p" ] && export RSTUDIO_PANDOC="$(dirname "$p")"
fi
echo "RSTUDIO_PANDOC=${RSTUDIO_PANDOC:-<not found>}"
echo "OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"

# Fail here rather than after the render. Pandoc is the LAST step, so without this the report
# computes the whole cohort -- over an hour -- and then dies with nothing written.
if [ -z "${RSTUDIO_PANDOC:-}" ] || [ ! -x "${RSTUDIO_PANDOC}/pandoc" ]; then
  echo "pandoc not found. rmarkdown needs it to write the HTML, and it is the last step of a" >&2
  echo "render that takes over an hour, so this stops now rather than at the end." >&2
  echo "Install it, or point RSTUDIO_PANDOC at a directory holding it." >&2
  exit 1
fi

IN="${1:-RNA_Parallel_linux.Rmd}"
OUT="${2:-${IN%.Rmd}.html}"

exec Rscript -e "rmarkdown::render('${IN}', output_file='${OUT}', quiet=FALSE)"
