$ErrorActionPreference = "Stop"
$projectRoot = Split-Path -Parent $PSScriptRoot
$python = Join-Path $projectRoot "envs\stereo_seq_py\Scripts\python.exe"
$script = Join-Path $projectRoot "scripts\1163_k7_spatial_tissue_niche_characterization.py"
$output = Join-Path $projectRoot "results\reanalysis\K7_SPATIAL_TISSUE_NICHE_CHARACTERIZATION"
$log = Join-Path $projectRoot "logs\1163_k7_spatial_tissue_niche_characterization.log"

& $python $script --output-dir $output 2>&1 | Tee-Object -FilePath $log
if ($LASTEXITCODE -ne 0) {
    throw "Spatial niche characterization failed with exit code $LASTEXITCODE"
}
