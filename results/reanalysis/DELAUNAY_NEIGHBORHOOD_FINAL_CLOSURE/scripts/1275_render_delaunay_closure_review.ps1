param(
    [string]$Quarto = 'D:\bb\RSTUDIO\resources\app\bin\quarto\bin\quarto.cmd'
)

$ErrorActionPreference = 'Stop'
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$closureDir = Split-Path -Parent $scriptDir
$timestamp = Get-Date -Format 'yyMMddHHmmss'
$outputName = "DELAUNAY_NEIGHBORHOOD_FINAL_CLOSURE_REVIEW_$timestamp.html"
$logPath = Join-Path $closureDir "logs\1275_render_review_$timestamp.log"
$statusPath = Join-Path $closureDir 'logs\1275_render_latest_status.json'

Push-Location $scriptDir
try {
    # Quarto writes normal progress messages to stderr; do not promote them to
    # terminating PowerShell errors. The native exit code remains authoritative.
    $ErrorActionPreference = 'Continue'
    & $Quarto render '1274_delaunay_closure_review.qmd' --output-dir '..\report' --output $outputName *>&1 |
        Tee-Object -FilePath $logPath
    if ($LASTEXITCODE -ne 0) {
        throw "Quarto render failed with exit code $LASTEXITCODE"
    }
    $outputPath = Join-Path $closureDir "report\$outputName"
    $status = [ordered]@{
        status = 'PASS'
        rendered_at = (Get-Date).ToString('o')
        output_path = $outputPath
        bytes = (Get-Item -LiteralPath $outputPath).Length
        render_log = $logPath
    }
    $status | ConvertTo-Json | Set-Content -LiteralPath $statusPath -Encoding UTF8
}
finally {
    Pop-Location
}
