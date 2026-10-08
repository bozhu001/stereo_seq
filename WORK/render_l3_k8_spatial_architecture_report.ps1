$ErrorActionPreference = "Stop"

$timestamp = Get-Date -Format "yyMMddHHmmss"
$projectRoot = (Get-Location).Path
$logPath = Join-Path $projectRoot "logs/L3_K8_SPATIAL_ARCHITECTURE_REPORT_$timestamp.log"
$errorLogPath = Join-Path $projectRoot "logs/L3_K8_SPATIAL_ARCHITECTURE_REPORT_$timestamp.stderr.log"
$source = Join-Path $projectRoot "scripts/1012_l3_k8_spatial_architecture_report.qmd"
$outputDir = Join-Path $projectRoot "results"
$generated = Join-Path $outputDir "1012_l3_k8_spatial_architecture_report.html"
$final = Join-Path $outputDir "L3_K8_SPATIAL_ARCHITECTURE_REPORT_$timestamp.html"

New-Item -ItemType Directory -Force -Path (Split-Path $logPath) | Out-Null
$quarto = (Get-Command quarto).Source
$render = Start-Process `
  -FilePath $quarto `
  -ArgumentList @("render", "`"$source`"", "--output-dir", "`"$outputDir`"") `
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
