param([string]$Python = "python")
$ErrorActionPreference = "Stop"

& $Python -c 'import sys; assert (3,12) <= sys.version_info[:2] < (3,15), "Python 3.12-3.14 required"'
if ($LASTEXITCODE -ne 0) { throw "Unsupported Python interpreter" }
& $Python (Join-Path $PSScriptRoot "scripts/verify_release.py")
if ($LASTEXITCODE -ne 0) { throw "Release verification failed" }

$venv = Join-Path $PSScriptRoot ".venv"
if (-not (Test-Path (Join-Path $venv "Scripts/python.exe"))) {
    & $Python -m venv $venv
    if ($LASTEXITCODE -ne 0) { throw "Virtual environment creation failed" }
}
$runtime = Join-Path $venv "Scripts/python.exe"
& $runtime -m pip install "${PSScriptRoot}[excel,3d]"
if ($LASTEXITCODE -ne 0) { throw "Dependency installation failed" }
& $runtime -m pip check
if ($LASTEXITCODE -ne 0) { throw "Dependency validation failed" }
& $runtime -c 'import gena; from load import load_data; from libs import ephemeris; from plots.scene import SceneSpec; print("GENA-CODE", gena.__version__)'
if ($LASTEXITCODE -ne 0) { throw "Import verification failed" }
Write-Output "Installed. Python interpreter: $runtime"
