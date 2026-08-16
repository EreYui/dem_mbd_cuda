[CmdletBinding()]
param(
    [Parameter(Position = 0)]
    [string]$Case = "DATA",

    [ValidateRange(1, 256)]
    [int]$Processors = 4,

    [ValidateSet("CUDA", "CPU")]
    [string]$Engine = "CUDA",

    [ValidateRange(0, 32)]
    [int]$Device = 0,

    [switch]$OneStep,

    [string]$EdemExe = "C:\Program Files\Altair\2024\EDEM\bin\edem.exe"
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$CaseRoot = Join-Path (Join-Path $RepoRoot "Data") $Case
$EdemDir = Join-Path $CaseRoot "edem"
$Deck = Join-Path $EdemDir "DATA_cube_drop.dem"
$ManifestPath = Join-Path $EdemDir "case_manifest.json"

if (-not (Test-Path -LiteralPath $EdemExe -PathType Leaf)) {
    throw "EDEM executable not found: $EdemExe"
}
if (-not (Test-Path -LiteralPath $Deck -PathType Leaf)) {
    throw "EDEM deck not found: $Deck`nRun .\edem\generate_case.ps1 $Case first."
}
if (-not (Test-Path -LiteralPath $ManifestPath -PathType Leaf)) {
    throw "Case manifest not found: $ManifestPath"
}

$Manifest = Get-Content -LiteralPath $ManifestPath -Raw | ConvertFrom-Json
$TimeStep = [double]$Manifest.simulation.time_step_s
$WriteInterval = [double]$Manifest.simulation.write_interval_s
$RunTime = [double]$Manifest.simulation.total_time_s
$EngineCode = if ($Engine -eq "CUDA") { "2" } else { "0" }
$Arguments = @(
    "--console",
    "-i", $Deck,
    "--rewind",
    "--time-step", $TimeStep.ToString("R", [Globalization.CultureInfo]::InvariantCulture),
    "--write-out", $WriteInterval.ToString("R", [Globalization.CultureInfo]::InvariantCulture),
    "--processors", $Processors,
    "--engine", $EngineCode,
    "--file-logger", (Join-Path $EdemDir "edem_run.log")
)
if ($Engine -eq "CUDA") {
    $Arguments += @("--device", $Device)
}
if ($OneStep) {
    $Arguments += @("--steps-total", "1")
}
else {
    $Arguments += @("--run-time", $RunTime.ToString("R", [Globalization.CultureInfo]::InvariantCulture))
}

& $EdemExe @Arguments
if ($LASTEXITCODE -ne 0) {
    throw "EDEM exited with code $LASTEXITCODE. See $(Join-Path $EdemDir 'edem_run.log')"
}
