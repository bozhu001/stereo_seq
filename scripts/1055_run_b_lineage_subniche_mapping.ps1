param(
    [Parameter(Mandatory = $true)]
    [string]$OutputDir
)

$ErrorActionPreference = "Stop"
$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$outputPath = Join-Path $projectRoot $OutputDir
New-Item -ItemType Directory -Path $outputPath -Force | Out-Null
New-Item -ItemType Directory -Path (Join-Path $outputPath "logs") -Force | Out-Null
$logPath = Join-Path $outputPath "logs/1055_b_lineage_subniche_mapping.log"
$pythonWsl = "/home/bozhu/miniconda3/envs/stereo/bin/python"

function Convert-ToWslPath([string]$Path) {
    $converted = & wsl.exe -d Ubuntu -- wslpath -a $Path
    if ($LASTEXITCODE -ne 0) { throw "Failed to convert path to WSL: $Path" }
    return $converted.Trim()
}

$rootWsl = Convert-ToWslPath $projectRoot
$outputWsl = Convert-ToWslPath $outputPath
$scriptWsl = "$rootWsl/scripts/1055_map_b_lineage_to_immune_subniches.py"
$maskWsl = "$rootWsl/results/reanalysis/B_LINEAGE_CALIBRATED_COVERAGE_MASK_GE5_ONLY_20261010_220617/01_FINAL_CALIBRATED_B_COVERAGE_MASK.tsv.gz"
$assignWsl = "$rootWsl/results/reanalysis/EXTENDED_IMMUNE_SIGNATURE_VALIDATION_20261010_192044/work/01_FROZEN_K345_ASSIGNMENTS.tsv.gz"

Start-Transcript -Path $logPath -Force
try {
    & wsl.exe -d Ubuntu -- $pythonWsl $scriptWsl `
        --mask $maskWsl `
        --assignments $assignWsl `
        --output-dir $outputWsl
    if ($LASTEXITCODE -ne 0) { throw "B-lineage subniche mapping failed" }
}
finally {
    Stop-Transcript
}
