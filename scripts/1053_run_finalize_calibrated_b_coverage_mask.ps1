param(
    [Parameter(Mandatory = $true)]
    [string]$OutputDir
)

$ErrorActionPreference = "Stop"
$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$outputPath = Join-Path $projectRoot $OutputDir
New-Item -ItemType Directory -Path $outputPath -Force | Out-Null
New-Item -ItemType Directory -Path (Join-Path $outputPath "logs") -Force | Out-Null
$logPath = Join-Path $outputPath "logs/1053_finalize_calibrated_b_coverage_mask.log"
$pythonWsl = "/home/bozhu/miniconda3/envs/stereo/bin/python"

function Convert-ToWslPath([string]$Path) {
    $converted = & wsl.exe -d Ubuntu -- wslpath -a $Path
    if ($LASTEXITCODE -ne 0) { throw "Failed to convert path to WSL: $Path" }
    return $converted.Trim()
}

$rootWsl = Convert-ToWslPath $projectRoot
$outputWsl = Convert-ToWslPath $outputPath
$scriptWsl = "$rootWsl/scripts/1053_finalize_calibrated_b_coverage_mask.py"
$auditWsl = "$rootWsl/results/reanalysis/B_LINEAGE_BIN_COVERAGE_RULE_AUDIT_20261010_212159"
$pilotWsl = "$rootWsl/results/reanalysis/B_LINEAGE_BIN_COVERAGE_PILOT_20261010_211348"
$baseWsl = "$rootWsl/results/reanalysis/EXTENDED_IMMUNE_SIGNATURE_VALIDATION_20261010_192044"

Start-Transcript -Path $logPath -Force
try {
    & wsl.exe -d Ubuntu -- $pythonWsl $scriptWsl `
        --panel-scores "$auditWsl/05_BIN_INFORMATIVE_PANEL_SCORES.tsv.gz" `
        --thresholds "$auditWsl/06_PATIENT_UMI_BACKGROUND_THRESHOLDS.tsv" `
        --pilot-mask "$pilotWsl/01_BIN_LEVEL_B_COVERAGE_MASK.tsv.gz" `
        --edges "$baseWsl/work/04_EXACT_RADIUS50_CANDIDATE_NEIGHBOR_EDGES.tsv.gz" `
        --output-dir $outputWsl
    if ($LASTEXITCODE -ne 0) { throw "Final calibrated B coverage mask failed" }
}
finally {
    Stop-Transcript
}
