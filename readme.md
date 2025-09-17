# CEE Scheduler Microservice + Desktop UI

This repository contains a complete course scheduling system built for the Department of Civil and Environmental Engineering at the University of Utah.
It uses an **Integer Programming model with Gurobi** to assign courses to time slots while respecting prerequisites, co-requisites, lab/discussion structures, and slot preferences.

Backend: **FastAPI + Docker**
UI: **PySide6 desktop app**

---

## Backend (Docker)

### Using Docker Compose (recommended)

1. Copy `.env.example` to `.env` and fill in your **Gurobi WLS credentials**:

   ```env
   GRB_WLSACCESSID=your-access-id
   GRB_WLSSECRET=your-secret
   GRB_LICENSEID=your-license-id
   ```

2. Start the API in production mode:

   ```powershell
   docker compose up api --build
   ```

   Or start in development mode (hot reload on code changes):

   ```powershell
   docker compose up api-dev --build
   ```

3. Check:

   * Health: [http://localhost:8000/health](http://localhost:8000/health)
   * API Docs: [http://localhost:8000/docs](http://localhost:8000/docs)

### Manual build/run (if you don't want compose)

```powershell
docker build -t cee-scheduler-api -f backend/Dockerfile .
docker run --rm -it `
  -p 8000:8000 `
  -v "${PWD}\data:/data" `
  -e GRB_WLSACCESSID="YOUR_WLSACCESSID" `
  -e GRB_WLSSECRET="YOUR_WLSSECRET" `
  -e GRB_LICENSEID="YOUR_LICENSEID" `
  cee-scheduler-api
```

---

## UI (PySide6)

### Run locally

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r ui\requirements.txt
python ui\main.py
```

### Build a Windows .exe

From PowerShell:

```powershell
cd ui
powershell -ExecutionPolicy Bypass -File .\build-win.ps1
# exe will be at: ui\dist\CEE-Scheduler.exe
```