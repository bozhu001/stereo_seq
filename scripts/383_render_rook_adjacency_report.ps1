$ErrorActionPreference = "Stop"

# Render-only wrapper for scripts/382_rook_adjacency_sensitivity_report.qmd.
# It reads completed TSV/JSON/PNG artifacts; it does not run statistical or
# spatial analysis and does not call scripts/380 or scripts/381.
$projectRoot = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $projectRoot

$timestamp = Get-Date -Format "yyMMddHHmmss"
$outputDirectory = Join-Path $projectRoot (
    "results/reanalysis/" +
    "bin50_all21_patient_pseudobulk_disease_signal_audit_260822/" +
    "rook_adjacency_sensitivity_260822200142"
)
$logDirectory = Join-Path $projectRoot (
    "results/reanalysis/" +
    "bin50_all21_patient_pseudobulk_disease_signal_audit_260822/logs"
)
$reportName = "rook_adjacency_sensitivity_report_${timestamp}.html"
$reportPath = Join-Path $outputDirectory $reportName
$fallbackReportPath = Join-Path $projectRoot $reportName
$logPath = Join-Path $logDirectory (
    "rook_adjacency_sensitivity_report_${timestamp}.log"
)

New-Item -ItemType Directory -Force -Path $outputDirectory, $logDirectory |
    Out-Null

$renderArguments = @(
    "render",
    "scripts/382_rook_adjacency_sensitivity_report.qmd",
    "--output", $reportName,
    "--output-dir", $outputDirectory
)

$startedAt = Get-Date
@(
    "Render started: $($startedAt.ToString('o'))"
    "Project root: $projectRoot"
    "Command: quarto $($renderArguments -join ' ')"
    "Analysis recomputation: disabled (report-only render)"
) | Set-Content -LiteralPath $logPath -Encoding UTF8

$previousErrorActionPreference = $ErrorActionPreference
$ErrorActionPreference = "Continue"
$renderOutput = & quarto @renderArguments 2>&1
$renderExitCode = $LASTEXITCODE
$ErrorActionPreference = $previousErrorActionPreference
$renderOutput | Add-Content -LiteralPath $logPath -Encoding UTF8

$recoveredFromQuartoCleanupFailure = $false
if (
    -not (Test-Path -LiteralPath $reportPath -PathType Leaf) -and
    (Test-Path -LiteralPath $fallbackReportPath -PathType Leaf)
) {
    Move-Item -LiteralPath $fallbackReportPath -Destination $reportPath
    $recoveredFromQuartoCleanupFailure = $true
}

if (-not (Test-Path -LiteralPath $reportPath -PathType Leaf)) {
    throw (
        "Quarto render failed with exit code $renderExitCode and the report " +
        "is missing. Log: $logPath"
    )
}

if ($renderExitCode -ne 0 -and -not $recoveredFromQuartoCleanupFailure) {
    throw "Quarto render failed with exit code $renderExitCode. Log: $logPath"
}

$report = Get-Item -LiteralPath $reportPath
if ($report.Length -le 0) {
    throw "Rendered report is empty: $reportPath"
}

# Quarto can leave absolute-input images as relative links when its Windows
# cleanup step fails before output relocation. Resolve those links against the
# QMD directory and embed the already-created PNGs so the HTML is standalone.
$html = [System.IO.File]::ReadAllText($reportPath)
$imageMatches = [regex]::Matches($html, '<img[^>]+src="([^"]+)"')
$embeddedImageCount = 0
foreach ($imageMatch in $imageMatches) {
    $source = $imageMatch.Groups[1].Value
    if ($source -match '^(data:|https?:|#)') {
        continue
    }

    $relativeSource = $source.Replace(
        '/',
        [System.IO.Path]::DirectorySeparatorChar
    )
    $sourcePath = [System.IO.Path]::GetFullPath(
        (Join-Path $PSScriptRoot $relativeSource)
    )
    if (-not (Test-Path -LiteralPath $sourcePath -PathType Leaf)) {
        throw "Local report image is missing: $sourcePath"
    }

    $extension = [System.IO.Path]::GetExtension($sourcePath).ToLowerInvariant()
    $mimeType = switch ($extension) {
        '.png' { 'image/png' }
        '.jpg' { 'image/jpeg' }
        '.jpeg' { 'image/jpeg' }
        '.svg' { 'image/svg+xml' }
        default { throw "Unsupported report image extension: $extension" }
    }
    $base64 = [System.Convert]::ToBase64String(
        [System.IO.File]::ReadAllBytes($sourcePath)
    )
    $dataUri = "data:${mimeType};base64,${base64}"
    $html = $html.Replace("src=`"${source}`"", "src=`"${dataUri}`"")
    $embeddedImageCount += 1
}
[System.IO.File]::WriteAllText(
    $reportPath,
    $html,
    [System.Text.UTF8Encoding]::new($false)
)
$report = Get-Item -LiteralPath $reportPath

$finishedAt = Get-Date
@(
    "Render finished: $($finishedAt.ToString('o'))"
    "Elapsed seconds: $([math]::Round(($finishedAt - $startedAt).TotalSeconds, 3))"
    "Output: $reportPath"
    "Output bytes: $($report.Length)"
    "Recovered after Quarto cleanup-only failure: $recoveredFromQuartoCleanupFailure"
    "Locally embedded images: $embeddedImageCount"
) | Add-Content -LiteralPath $logPath -Encoding UTF8

Write-Output "REPORT_PATH=$reportPath"
Write-Output "LOG_PATH=$logPath"
