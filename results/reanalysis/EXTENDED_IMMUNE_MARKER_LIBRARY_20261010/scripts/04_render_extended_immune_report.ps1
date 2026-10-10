param()

$ErrorActionPreference = "Stop"

$outputDir = Split-Path -Parent $PSScriptRoot
$source = Join-Path $PSScriptRoot "03_extended_immune_marker_report.qmd"
$reports = Join-Path $outputDir "reports"
$logs = Join-Path $outputDir "logs"
$quarto = "D:\Positron\resources\app\quarto\bin\quarto.exe"
$stamp = Get-Date -Format "yyMMddHHmmss"
$htmlName = "EXTENDED_IMMUNE_MARKER_LIBRARY_REPORT_$stamp.html"
$htmlPath = Join-Path $reports $htmlName
$logPath = Join-Path $logs "04_render_extended_immune_report_$stamp.log"
$statusPath = Join-Path $logs "04_render_extended_immune_report_latest.tsv"

New-Item -ItemType Directory -Path $reports -Force | Out-Null
New-Item -ItemType Directory -Path $logs -Force | Out-Null

"COMMAND`t& '$quarto' render '$source' --output '$htmlName' --output-dir '$reports'" |
    Set-Content -LiteralPath $logPath -Encoding utf8
"QUARTO_VERSION`t$(& $quarto --version)" |
    Add-Content -LiteralPath $logPath -Encoding utf8
"POWERSHELL_VERSION`t$($PSVersionTable.PSVersion)" |
    Add-Content -LiteralPath $logPath -Encoding utf8

try {
    $stdoutPath = "$logPath.stdout"
    $stderrPath = "$logPath.stderr"
    $arguments = @(
        "render",
        "`"$source`"",
        "--output",
        "`"$htmlName`"",
        "--output-dir",
        "`"$reports`""
    )
    $startArgs = @{
        FilePath = $quarto
        ArgumentList = $arguments
        Wait = $true
        PassThru = $true
        NoNewWindow = $true
        WorkingDirectory = $PSScriptRoot
        RedirectStandardOutput = $stdoutPath
        RedirectStandardError = $stderrPath
    }
    $process = Start-Process @startArgs
    if (Test-Path -LiteralPath $stdoutPath) {
        Get-Content -LiteralPath $stdoutPath | Add-Content -LiteralPath $logPath
        Remove-Item -LiteralPath $stdoutPath -Force
    }
    if (Test-Path -LiteralPath $stderrPath) {
        Get-Content -LiteralPath $stderrPath | Add-Content -LiteralPath $logPath
        Remove-Item -LiteralPath $stderrPath -Force
    }
    if ($process.ExitCode -ne 0) {
        throw "Quarto exited with code $($process.ExitCode)"
    }
    if (-not (Test-Path -LiteralPath $htmlPath)) {
        throw "Expected report was not created: $htmlPath"
    }
    $size = (Get-Item -LiteralPath $htmlPath).Length
    if ($size -le 0) {
        throw "Rendered report is empty: $htmlPath"
    }
    "status`treport_path`tlog_path`tsize_bytes`ttimestamp" |
        Set-Content -LiteralPath $statusPath -Encoding utf8
    "PASS`t$htmlPath`t$logPath`t$size`t$stamp" |
        Add-Content -LiteralPath $statusPath -Encoding utf8
    "REPORT_COMPLETE=$htmlPath"
    "SESSION_INFO=Quarto $(& $quarto --version); PowerShell $($PSVersionTable.PSVersion)"
} catch {
    "status`treport_path`tlog_path`tsize_bytes`ttimestamp" |
        Set-Content -LiteralPath $statusPath -Encoding utf8
    "FAIL`t$htmlPath`t$logPath`t0`t$stamp" |
        Add-Content -LiteralPath $statusPath -Encoding utf8
    $_ | Out-String | Add-Content -LiteralPath $logPath -Encoding utf8
    throw
}
