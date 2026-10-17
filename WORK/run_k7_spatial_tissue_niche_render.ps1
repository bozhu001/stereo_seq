$ErrorActionPreference = "Stop"

$repo = Split-Path -Parent $PSScriptRoot
$python = Join-Path $repo "envs\stereo_seq_py\Scripts\python.exe"
$script = Join-Path $repo "scripts\1164_render_k7_spatial_tissue_niche_atlas.py"
$output = Join-Path $repo "results\reanalysis\K7_SPATIAL_TISSUE_NICHE_CHARACTERIZATION"
$log = Join-Path $repo "logs\1164_render_k7_spatial_tissue_niche_atlas.log"

& $python $script --output-dir $output *>&1 | Tee-Object -FilePath $log
if ($LASTEXITCODE -ne 0) {
    throw "Atlas rendering failed with exit code $LASTEXITCODE"
}
