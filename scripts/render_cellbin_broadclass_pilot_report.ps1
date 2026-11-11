$ErrorActionPreference = "Stop"

$projectRoot = Split-Path -Parent $PSScriptRoot
$timestamp = Get-Date -Format "yyMMddHHmmss"
$qmd = Join-Path $PSScriptRoot "cellbin_broadclass_pilot_21samples.qmd"
$results = Join-Path $projectRoot "results"
$logs = Join-Path $projectRoot "logs"
$baseHtml = Join-Path $results "cellbin_broadclass_pilot_21samples.html"
$finalHtml = Join-Path $results "cellbin_broadclass_pilot_21samples_$timestamp.html"
$log = Join-Path $logs "cellbin_broadclass_pilot_report_$timestamp.log"
$status = Join-Path $results "cellbin_broadclass_pilot_report_latest_status.tsv"

$env:QUARTO_PYTHON = Join-Path $projectRoot "envs\stereo_seq_py\Scripts\python.exe"

try {
    @(
        "timestamp=$timestamp",
        "qmd=$qmd",
        "output_dir=$results",
        "quarto_python=$env:QUARTO_PYTHON",
        "command=quarto render <qmd> --output-dir <results> --execute-daemon 0 --no-clean --quiet"
    ) | Set-Content -LiteralPath $log -Encoding UTF8
    $previousErrorPreference = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    & quarto render $qmd --output-dir $results --execute-daemon 0 --no-clean --quiet 2>&1 |
        Tee-Object -FilePath $log -Append
    $exitCode = if ($null -eq $LASTEXITCODE) { 0 } else { $LASTEXITCODE }
    $ErrorActionPreference = $previousErrorPreference
    if ($exitCode -ne 0) {
        throw "Quarto exited with code $exitCode"
    }
    Move-Item -LiteralPath $baseHtml -Destination $finalHtml
    @(
        "field`tvalue",
        "status`tPASS",
        "timestamp`t$timestamp",
        "html`t$finalHtml",
        "log`t$log"
    ) | Set-Content -LiteralPath $status -Encoding UTF8
}
catch {
    $errorMessage = ($_ | Out-String).Trim().Replace("`t", " ").Replace("`r", " ").Replace("`n", " ")
    @(
        "field`tvalue",
        "status`tFAIL",
        "timestamp`t$timestamp",
        "html`t$finalHtml",
        "log`t$log",
        "error`t$errorMessage"
    ) | Set-Content -LiteralPath $status -Encoding UTF8
    throw
}
