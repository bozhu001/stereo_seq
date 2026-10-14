param(
    [string]$Source = "results\reanalysis\K6_FACTOR_INTERNAL_DISEASE_20261013_143952",
    [string]$Output = ""
)

$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
$python = Join-Path $projectRoot "envs\stereo_seq_py\Scripts\python.exe"
if ([string]::IsNullOrWhiteSpace($Output)) {
    $timestamp = Get-Date -Format "yyyyMMdd_HHmmss"
    $Output = "results\reanalysis\CROSS_CHIP_EDGER_STAGE1_CLOSEOUT_$timestamp"
}
$logTimestamp = Get-Date -Format "yyyyMMdd_HHmmss"
$log = Join-Path $projectRoot "logs\1116_crosschip_edger_closeout_$logTimestamp.log"

& $python (Join-Path $PSScriptRoot "1116_closeout_crosschip_edger_evidence.py") `
    --project-root $projectRoot `
    --source $Source `
    --output $Output 2>&1 | Tee-Object -FilePath $log

if ($LASTEXITCODE -ne 0) {
    throw "1116 cross-chip edgeR closeout failed with exit code $LASTEXITCODE"
}

Write-Output (Join-Path $projectRoot $Output)
Write-Output $log
