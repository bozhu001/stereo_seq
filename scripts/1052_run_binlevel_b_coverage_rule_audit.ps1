param(
    [Parameter(Mandatory = $true)]
    [string]$OutputDir
)

$ErrorActionPreference = "Stop"
$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$outputPath = Join-Path $projectRoot $OutputDir
New-Item -ItemType Directory -Path $outputPath -Force | Out-Null
New-Item -ItemType Directory -Path (Join-Path $outputPath "logs") -Force | Out-Null
$logPath = Join-Path $outputPath "logs/1052_binlevel_b_coverage_rule_audit.log"
$pythonWsl = "/home/bozhu/miniconda3/envs/stereo/bin/python"

function Convert-ToWslPath([string]$Path) {
    $converted = & wsl.exe -d Ubuntu -- wslpath -a $Path
    if ($LASTEXITCODE -ne 0) {
        throw "Failed to convert path to WSL: $Path"
    }
    return $converted.Trim()
}

$rootWsl = Convert-ToWslPath $projectRoot
$outputWsl = Convert-ToWslPath $outputPath
$scriptWsl = "$rootWsl/scripts/1052_audit_binlevel_b_coverage_rule.py"
$baseWsl = "$rootWsl/results/reanalysis/EXTENDED_IMMUNE_SIGNATURE_VALIDATION_20261010_192044"

Start-Transcript -Path $logPath -Force
try {
    & wsl.exe -d Ubuntu -- $pythonWsl $scriptWsl `
        --score-file "$baseWsl/work/06_CENTER_NEIGHBOR_MARKER_SCORES.tsv.gz" `
        --assignments "$baseWsl/work/01_FROZEN_K345_ASSIGNMENTS.tsv.gz" `
        --cluster-validation "$baseWsl/04_CLUSTER_LEVEL_EXTENDED_SIGNATURE_VALIDATION.tsv" `
        --output-dir $outputWsl
    if ($LASTEXITCODE -ne 0) {
        throw "B coverage rule audit failed"
    }
}
finally {
    Stop-Transcript
}
