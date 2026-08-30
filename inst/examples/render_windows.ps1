# Render the Windows verification report for rnaparallel.
#
# A multi-threaded BLAS reads its thread count when the library loads, so these have to be in
# the environment BEFORE R starts. Sys.setenv() inside the session is too late, which is the
# same reason render_linux.sh exports rather than sets. Stock Windows R ships the single-threaded
# reference BLAS and none of this matters there. It matters on a user-installed OpenBLAS or MKL
# build, where every PSOCK worker would otherwise open its own thread pool on top of the cluster
# and oversubscribe the machine.
#
#   .\render_windows.ps1
#   .\render_windows.ps1 RNA_Parallel_windows.Rmd out.html

param(
  [string]$In  = "RNA_Parallel_windows.Rmd",
  [string]$Out = ""
)

$ErrorActionPreference = "Stop"
Set-Location -Path $PSScriptRoot

if (-not $Out) { $Out = [System.IO.Path]::ChangeExtension($In, "html") }
if (-not (Test-Path $In)) { throw "input not found: $In" }

$env:OPENBLAS_NUM_THREADS = "1"
$env:OMP_NUM_THREADS      = "1"
$env:MKL_NUM_THREADS      = "1"

# pandoc ships with RStudio and is not on PATH in a plain PowerShell session. A machine with no
# RStudio at all is the other common case here, so a standalone install and PATH both count.
if (-not $env:RSTUDIO_PANDOC) {
  $roots = @($env:ProgramFiles, ${env:ProgramFiles(x86)}, "$env:LOCALAPPDATA\Programs") |
           Where-Object { $_ }
  $tails = @(
    "RStudio\resources\app\bin\quarto\bin\tools",
    "RStudio\resources\app\bin\pandoc",
    "RStudio\bin\quarto\bin\tools",
    "RStudio\bin\pandoc"
  )
  foreach ($r in $roots) {
    foreach ($t in $tails) {
      $d = Join-Path $r $t
      if (Test-Path (Join-Path $d "pandoc.exe")) { $env:RSTUDIO_PANDOC = $d; break }
    }
    if ($env:RSTUDIO_PANDOC) { break }
  }
}
if (-not $env:RSTUDIO_PANDOC) {
  $onPath = (Get-Command pandoc.exe -ErrorAction SilentlyContinue).Source
  if ($onPath) {
    $env:RSTUDIO_PANDOC = Split-Path $onPath
  } else {
    foreach ($d in @("$env:LOCALAPPDATA\Pandoc", "$env:ProgramFiles\Pandoc")) {
      if (Test-Path (Join-Path $d "pandoc.exe")) { $env:RSTUDIO_PANDOC = $d; break }
    }
  }
}

# Rscript is not on PATH by default on Windows, unlike every other platform this package
# renders on. The registry knows where R is even when the shell does not.
$rscript = (Get-Command Rscript.exe -ErrorAction SilentlyContinue).Source
if (-not $rscript) {
  $key = Get-ItemProperty "HKLM:\SOFTWARE\R-core\R" -ErrorAction SilentlyContinue
  if ($key.InstallPath) {
    $cand = Join-Path $key.InstallPath "bin\x64\Rscript.exe"
    if (Test-Path $cand) { $rscript = $cand }
  }
}
if (-not $rscript) {
  $cand = Get-ChildItem "$env:ProgramFiles\R" -Filter Rscript.exe -Recurse -ErrorAction SilentlyContinue |
          Sort-Object FullName -Descending | Select-Object -First 1
  if ($cand) { $rscript = $cand.FullName }
}
if (-not $rscript) { throw "Rscript.exe not found. Add R's bin directory to PATH." }

Write-Host "Rscript        = $rscript"
Write-Host "RSTUDIO_PANDOC = $(if ($env:RSTUDIO_PANDOC) { $env:RSTUDIO_PANDOC } else { '<not found>' })"

# Fail here rather than after the render. Pandoc is the LAST step, so without this the report
# computes the whole cohort -- over an hour -- and then dies with nothing written.
if (-not $env:RSTUDIO_PANDOC -or -not (Test-Path (Join-Path $env:RSTUDIO_PANDOC "pandoc.exe"))) {
  throw ("pandoc not found. rmarkdown needs it to write the HTML, and it is the last step of " +
         "a render that takes over an hour, so this stops now rather than at the end. " +
         "Install it, or set RSTUDIO_PANDOC to a directory holding pandoc.exe.")
}
Write-Host "BLAS threads   = $env:OPENBLAS_NUM_THREADS"
Write-Host "rendering $In -> $Out"

& $rscript -e "rmarkdown::render('$In', output_file='$Out', quiet=FALSE)"
exit $LASTEXITCODE
