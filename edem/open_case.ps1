[CmdletBinding()]
param(
    [Parameter(Position = 0)]
    [string]$Case = "DATA",

    [string]$EdemExe = "C:\Program Files\Altair\2024\EDEM\bin\edem.exe"
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Deck = Join-Path (Join-Path (Join-Path $RepoRoot "Data") $Case) "edem\DATA_cube_drop.dem"

if (-not (Test-Path -LiteralPath $EdemExe -PathType Leaf)) {
    throw "EDEM executable not found: $EdemExe"
}
if (-not (Test-Path -LiteralPath $Deck -PathType Leaf)) {
    throw "EDEM deck not found: $Deck`nRun .\edem\generate_case.ps1 $Case first."
}

Start-Process -FilePath $EdemExe -ArgumentList @("-i", ('"{0}"' -f $Deck))
