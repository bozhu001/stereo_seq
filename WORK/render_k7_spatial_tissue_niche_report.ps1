$ErrorActionPreference = "Stop"

$repo = Split-Path -Parent $PSScriptRoot
$python = Join-Path $repo "envs\stereo_seq_py\Scripts\python.exe"
$script = Join-Path $repo "scripts\1165_render_k7_spatial_tissue_niche_report.py"
$analysis = Join-Path $repo "results\reanalysis\K7_SPATIAL_TISSUE_NICHE_CHARACTERIZATION"
$stamp = Get-Date -Format "yyMMddHHmmss"
$log = Join-Path $repo "logs\1165_k7_spatial_tissue_niche_report_$stamp.log"

& $python $script --analysis-dir $analysis *>&1 | Tee-Object -FilePath $log
if ($LASTEXITCODE -ne 0) {
    throw "HTML report rendering failed with exit code $LASTEXITCODE"
}
