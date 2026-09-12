$ErrorActionPreference = "Stop"

$projectRoot = Split-Path -Parent $PSScriptRoot
$timestamp = Get-Date -Format "yyMMddHHmmss"
$quarto = (Get-Command quarto -ErrorAction Stop).Source
$source = "scripts/629_cdr3_per_bin_top20_heatmap_report.qmd"
$resultsDir = Join-Path $projectRoot "results"
$logsDir = Join-Path $projectRoot "logs"
$outputName = "CDR3_PER_BIN_TOP20_HEATMAP_REPORT_$timestamp.html"
$logPath = Join-Path $logsDir "cdr3_per_bin_top20_heatmap_report_render_$timestamp.log"
$commandPath = Join-Path $logsDir "cdr3_per_bin_top20_heatmap_report_command_$timestamp.txt"

New-Item -ItemType Directory -Force -Path $resultsDir, $logsDir | Out-Null

$commandText = @"
cd scripts
quarto render 629_cdr3_per_bin_top20_heatmap_report.qmd --output $outputName --output-dir ../results --no-clean
"@
Set-Content -LiteralPath $commandPath -Value $commandText -Encoding UTF8

Push-Location $PSScriptRoot
try {
  $savedErrorActionPreference = $ErrorActionPreference
  $ErrorActionPreference = "Continue"
  try {
    & $quarto render (Split-Path -Leaf $source) --output $outputName --output-dir ../results --no-clean *>&1 |
      Tee-Object -FilePath $logPath
    $renderExitCode = $LASTEXITCODE
  } finally {
    $ErrorActionPreference = $savedErrorActionPreference
  }
  if ($renderExitCode -ne 0) {
    throw "Quarto render failed with exit code $renderExitCode"
  }
} finally {
  Pop-Location
}

$outputPath = Join-Path $resultsDir $outputName
if (-not (Test-Path -LiteralPath $outputPath -PathType Leaf)) {
  throw "Expected report was not created: $outputPath"
}

$item = Get-Item -LiteralPath $outputPath
$hash = (Get-FileHash -LiteralPath $outputPath -Algorithm SHA256).Hash
Write-Output ("REPORT`t{0}`t{1}`t{2:yyyy-MM-ddTHH:mm:sszzz}`t{3}" -f `
  $item.FullName, $item.Length, $item.LastWriteTime, $hash)
