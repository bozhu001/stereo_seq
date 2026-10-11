$ErrorActionPreference = "Stop"
Set-Location -LiteralPath "D:\Lung ST\codex agent\stereo_seq\stereo_seq"
$stamp = Get-Date -Format "yyMMddHHmmss"
$logPath = "logs\1070_k3_niche_report_$stamp.log"
& ".\envs\stereo_seq_py\Scripts\python.exe" `
  "scripts\1070_render_k3_niche_report.py" `
  --atlas-dir "results\reanalysis\ALL21_SHARED_UNBIASED_DISCRETE_NICHE_ATLAS_20261011_191652_attempt4" `
  --heatmap-dir "results\reanalysis\ALL21_K3_NICHE_HEATMAP_DISEASE_COMPARISON_20261011_192920_attempt2" `
  --adjacency-dir "results\reanalysis\ALL21_K3_DIFFERENTIAL_NICHE_ADJACENCY_20261011_195242" `
  --output-dir "results" *> $logPath
if ($LASTEXITCODE -ne 0) {
  throw "Report render failed; see $logPath"
}
