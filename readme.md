# CEE Scheduler Microservice + Desktop UI

FastAPI backend + Docker + PySide6 UI.

## Backend (Docker)

1. Build:
```powershell
docker build -t cee-scheduler-api -f backend/Dockerfile .
```

2. Run:
```powershell
docker run --rm -it `
  -p 8000:8000 `
  -v "${PWD}\data:/data" `
  -e GRB_WLSACCESSID="YOUR_WLSACCESSID" `
  -e GRB_WLSSECRET="YOUR_WLSSECRET" `
  -e GRB_LICENSEID="YOUR_LICENSEID" `
  cee-scheduler-api
```

Check: http://localhost:8000/health  
Docs: http://localhost:8000/docs

## UI (PySide6)

Install and run locally:
```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r ui\requirements.txt
python ui\main.py
```

## Build a Windows .exe for the UI

From Windows PowerShell:
```powershell
cd ui
powershell -ExecutionPolicy Bypass -File .\build-win.ps1
# exe will be at: ui\dist\CEE-Scheduler.exe
```

Options:
```powershell
.\build-win.ps1 -Console                 # show console window
.\build-win.ps1 -Name "CEE-Scheduler-Dev"
.\build-win.ps1 -Python "C:\Python311\python.exe"
```
