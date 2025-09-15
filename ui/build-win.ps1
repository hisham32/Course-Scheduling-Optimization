Param(
    [string]$Python = "",               # autodetect if empty
    [string]$Name = "CEE-Scheduler",
    [switch]$Console
)

$ErrorActionPreference = "Stop"
$Here = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $Here   # we are now in ...\ui

# --- Pick a Python ---
if (-not $Python) {
    if (Get-Command py -ErrorAction SilentlyContinue) {
        $Python = "py"
    } elseif (Get-Command python -ErrorAction SilentlyContinue) {
        $Python = "python"
    } else {
        throw "No Python interpreter found. Install Python 3.10+ and ensure 'py' or 'python' is on PATH."
    }
}

Write-Host "==> Creating venv using: $Python" -ForegroundColor Cyan
& $Python -m venv .venv

$venvActivate = ".\.venv\Scripts\Activate.ps1"
if (-not (Test-Path $venvActivate)) {
    throw "Failed to create venv at .venv"
}
. $venvActivate

Write-Host "==> Upgrading pip" -ForegroundColor Cyan
python -m pip install --upgrade pip

Write-Host "==> Installing UI requirements + PyInstaller" -ForegroundColor Cyan
pip install -r requirements.txt
pip install pyinstaller==6.10.0

# Clean previous builds
if (Test-Path dist) { Remove-Item -Recurse -Force dist }
if (Test-Path build) { Remove-Item -Recurse -Force build }
if (Test-Path "$Name.spec") { Remove-Item -Force "$Name.spec" }

$consoleFlag = "--noconsole"
if ($Console) { $consoleFlag = "" }

Write-Host "==> Building $Name.exe with PyInstaller" -ForegroundColor Cyan
# We are already in the ui/ folder, so point to main.py (NOT ui\main.py)
pyinstaller `
    main.py `
    --name "$Name" `
    --onefile `
    $consoleFlag `
    --clean `
    --collect-all PySide6 `
    --collect-data certifi `
    --hidden-import PySide6.QtCore `
    --hidden-import PySide6.QtGui `
    --hidden-import PySide6.QtWidgets

if (-not (Test-Path "dist\$Name.exe")) {
    throw "Build failed. Check PyInstaller output above."
}

Write-Host "==> Build complete:" -ForegroundColor Green
Write-Host "    $(Resolve-Path dist\$Name.exe)"
