param(
    [string]$OutputStamp = (Get-Date -Format "yyyyMMdd_HHmmss")
)

$ErrorActionPreference = "Continue"
$projectRoot = Split-Path -Parent $PSScriptRoot
$python = Join-Path $projectRoot "envs\stereo_seq_py\Scripts\python.exe"
$output = Join-Path $projectRoot "results\reanalysis\NON_B_EXTERNAL_MICROENVIRONMENT_ATLAS_$OutputStamp"
$log = Join-Path $projectRoot "logs\1108_non_b_external_microenvironment_atlas_$OutputStamp.log"

$previousNativePreference = $ErrorActionPreference
& $python (Join-Path $PSScriptRoot "1108_build_non_b_external_microenvironment_atlas.py") `
    --project-root $projectRoot `
    --output $output 2>&1 | Tee-Object -FilePath $log

$analysisExitCode = $LASTEXITCODE
$ErrorActionPreference = $previousNativePreference
if ($analysisExitCode -ne 0) {
    throw "1108 analysis failed with exit code $analysisExitCode"
}
