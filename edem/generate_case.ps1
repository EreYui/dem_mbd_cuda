[CmdletBinding()]
param(
    [Parameter(Position = 0)]
    [string]$Case = "DATA",

    [string]$EdemRoot = "C:\Program Files\Altair\2024\EDEM"
)

$ErrorActionPreference = "Stop"
$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$PythonExe = "C:\Program Files\Altair\2024\common\python\python3.8\win64\python.exe"
$Wheel = Join-Path $EdemRoot "EDEMpy\win64\edempy-1.5.0-cp38-cp38-win_amd64.whl"
$Runtime = Join-Path $env:TEMP "dem_mbd_cuda_edempy_1_5_0"
$PureLib = Join-Path $Runtime "edempy-1.5.0.data\purelib"
$InstalledLib = Join-Path $EdemRoot "python"

if (-not (Test-Path -LiteralPath $PythonExe -PathType Leaf)) {
    throw "Altair Python 3.8 not found: $PythonExe"
}
if (-not (Test-Path -LiteralPath $Wheel -PathType Leaf)) {
    throw "EDEMpy wheel not found: $Wheel"
}
if (-not (Test-Path -LiteralPath (Join-Path $InstalledLib "edempy\Deck.py") -PathType Leaf) -and
    -not (Test-Path -LiteralPath (Join-Path $PureLib "edempy") -PathType Container)) {
    New-Item -ItemType Directory -Path $Runtime -Force | Out-Null
    tar -xf $Wheel -C $Runtime
}

$PreviousPythonPath = $env:PYTHONPATH
try {
    if (Test-Path -LiteralPath (Join-Path $InstalledLib "edempy\Deck.py") -PathType Leaf) {
        $env:PYTHONPATH = $InstalledLib
    }
    else {
        $env:PYTHONPATH = $PureLib
    }
    & $PythonExe (Join-Path $PSScriptRoot "generate_case.py") `
        --repo $RepoRoot `
        --case $Case `
        --edem-root $EdemRoot
    if ($LASTEXITCODE -ne 0) {
        throw "EDEM case generator exited with code $LASTEXITCODE"
    }
}
finally {
    $env:PYTHONPATH = $PreviousPythonPath
}
