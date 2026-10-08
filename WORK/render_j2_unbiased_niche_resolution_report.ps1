$ErrorActionPreference = "Stop"

$timestamp = Get-Date -Format "yyMMddHHmmss"
$projectRoot = (Get-Location).Path
$logPath = Join-Path $projectRoot "logs/J2_UNBIASED_NICHE_RESOLUTION_REPORT_$timestamp.log"
$errorLogPath = Join-Path $projectRoot "logs/J2_UNBIASED_NICHE_RESOLUTION_REPORT_$timestamp.stderr.log"
$source = Join-Path $projectRoot "scripts/1017_j2_unbiased_niche_resolution_report.qmd"
$outputDir = Join-Path $projectRoot "results"
$generated = Join-Path $outputDir "1017_j2_unbiased_niche_resolution_report.html"
$final = Join-Path $outputDir "J2_UNBIASED_NICHE_RESOLUTION_REPORT_$timestamp.html"
$quarto = "D:/bb/RSTUDIO/resources/app/bin/quarto/bin/quarto.exe"

New-Item -ItemType Directory -Force -Path (Split-Path $logPath) | Out-Null
$render = Start-Process `
  -FilePath $quarto `
  -ArgumentList @(
    "render",
    "`"$source`"",
    "--output-dir",
    "`"$outputDir`"",
    "--no-clean"
  ) `
  -WorkingDirectory $projectRoot `
  -WindowStyle Hidden `
  -RedirectStandardOutput $logPath `
  -RedirectStandardError $errorLogPath `
  -Wait `
  -PassThru

if ($render.ExitCode -ne 0) {
  throw "Quarto render failed with exit code $($render.ExitCode); see $logPath and $errorLogPath"
}

Move-Item -LiteralPath $generated -Destination $final
"REPORT=$final" | Add-Content -LiteralPath $logPath
"COMPLETE" | Set-Content -LiteralPath "$final.complete"
Write-Output "REPORT=$final"
