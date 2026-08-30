#!/usr/bin/env bash
# Render the macOS verification report for rnaparallel.
#
# Accelerate is the system BLAS here and manages its own thread pool, so there is no
# OPENBLAS_NUM_THREADS to pin as there is on Linux. What this guards is the other case: a
# user-installed OpenBLAS or MKL build of R, where every forked worker would open its own
# thread pool on top of the fork and oversubscribe the machine. A BLAS reads its thread count
# when the library LOADS, so these have to be in the environment before R starts --
# Sys.setenv() inside the session is too late, which is why render_linux.sh exports rather
# than sets. All four are harmless where the library they name is absent.
#
#   ./render_macos.sh
#   ./render_macos.sh RNA_Parallel.Rmd ../../docs/index.html
set -euo pipefail
cd "$(dirname "$0")"

export VECLIB_MAXIMUM_THREADS=1   # Accelerate
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# pandoc ships inside RStudio.app and is not on PATH in a plain shell. A machine with no
# RStudio at all is the other common case, so Homebrew, an active conda environment and PATH
# all count. Set RSTUDIO_PANDOC yourself to override the search.
if [ -z "${RSTUDIO_PANDOC:-}" ]; then
  for d in "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools" \
           "/Applications/RStudio.app/Contents/Resources/app/bin/quarto/bin/tools" \
           "/Applications/RStudio.app/Contents/MacOS/quarto/bin/tools" \
           "/Applications/RStudio.app/Contents/MacOS/pandoc" \
           "/opt/homebrew/bin" "/usr/local/bin" "${CONDA_PREFIX:-}/bin"; do
    [ -n "$d" ] && [ -x "$d/pandoc" ] && export RSTUDIO_PANDOC="$d" && break
  done
fi
if [ -z "${RSTUDIO_PANDOC:-}" ]; then
  p="$(command -v pandoc || true)"
  [ -n "$p" ] && export RSTUDIO_PANDOC="$(dirname "$p")"
fi

# Fail here rather than after the render. Pandoc is the LAST step, so without this the report
# computes the whole cohort -- an hour or more -- and then dies with nothing written.
if [ -z "${RSTUDIO_PANDOC:-}" ] || [ ! -x "${RSTUDIO_PANDOC}/pandoc" ]; then
  echo "pandoc not found. rmarkdown needs it to write the HTML, and it is the last step of a" >&2
  echo "render that takes over an hour, so this stops now rather than at the end." >&2
  echo "Install it (brew install pandoc) or point RSTUDIO_PANDOC at a directory holding it." >&2
  exit 1
fi

IN="${1:-RNA_Parallel.Rmd}"
OUT="${2:-${IN%.Rmd}.html}"

echo "RSTUDIO_PANDOC        = ${RSTUDIO_PANDOC:-<not found>}"
echo "VECLIB_MAXIMUM_THREADS= $VECLIB_MAXIMUM_THREADS"
echo "performance cores     = $(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || echo '?')"
echo "rendering $IN -> $OUT"

exec Rscript -e "rmarkdown::render('${IN}', output_file='${OUT}', quiet=FALSE)"
