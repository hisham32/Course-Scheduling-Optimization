import os
import json
import time
import base64
import threading
from fastapi import FastAPI, HTTPException
from fastapi.responses import StreamingResponse
from fastapi.middleware.cors import CORSMiddleware
from .models import OptimizeRequest, OptimizeResponse, ScheduleRow
from .scheduler import run_and_optionally_write_excel

app = FastAPI(title="CEE Scheduler Microservice", version="1.4.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], allow_credentials=True, allow_methods=["*"], allow_headers=["*"],
)

LOG_PATH = "/tmp/gurobi.log"


@app.get("/health")
def health():
    return {"status": "ok"}


@app.post("/optimize", response_model=OptimizeResponse)
def optimize(req: OptimizeRequest):
    try:
        if not req.output_excel_path and not req.return_excel_bytes:
            raise HTTPException(
                status_code=400,
                detail="You must either provide output_excel_path or set return_excel_bytes=True."
            )

        df_final, excel_bytes, log_text = run_and_optionally_write_excel(
            plan_semester=req.plan_semester.value,
            offerings_path=req.offerings_path,
            input_data_path=req.input_data_path,
            output_excel_path=req.output_excel_path,
        )
        rows = [ScheduleRow(**r) for r in df_final.to_dict(orient="records")]

        resp = OptimizeResponse(
            status="success",
            message="Optimization complete",
            schedule=rows,
            excel_path=req.output_excel_path,
            log=log_text,
        )
        if req.return_excel_bytes and excel_bytes is not None:
            resp.excel_base64 = base64.b64encode(excel_bytes).decode("utf-8")
        return resp

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/optimize/stream")
def optimize_stream(req: OptimizeRequest):
    """
    Streams Gurobi log lines via SSE while optimization runs.
    Emits:
      - {"type":"log","text":"..."} repeatedly
      - {"type":"done","status":"success","rows":[...], "excel_base64": "..."} once at end
    """
    # if not req.output_excel_path and not req.return_excel_bytes:
    #     raise HTTPException(
    #         status_code=400,
    #         detail="You must either provide output_excel_path or set return_excel_bytes=True."
    #     )

    # Truncate previous log BEFORE starting worker, to avoid tailing old content
    try:
        with open(LOG_PATH, "w", encoding="utf-8") as _f:
            pass
    except Exception:
        pass

    final = {"done": False, "ok": False, "rows": None, "excel_b64": None, "error": None}

    def worker():
        try:
            df_final, excel_bytes, _ = run_and_optionally_write_excel(
                plan_semester=req.plan_semester.value,
                offerings_path=req.offerings_path,
                input_data_path=req.input_data_path,
                output_excel_path=req.output_excel_path,
            )
            final["rows"] = df_final.to_dict(orient="records")
            if req.return_excel_bytes and excel_bytes:
                final["excel_b64"] = base64.b64encode(excel_bytes).decode("utf-8")
            final["ok"] = True
        except Exception as e:
            final["error"] = str(e)
            final["ok"] = False
        finally:
            final["done"] = True

    t = threading.Thread(target=worker, daemon=True)
    t.start()

    def sse_iter():
        # wait for log file to show up
        start = time.time()
        f = None
        while not os.path.exists(LOG_PATH) and not final["done"]:
            if time.time() - start > 5:
                break
            time.sleep(0.1)

        heartbeat_last = 0.0
        try:
            if os.path.exists(LOG_PATH):
                f = open(LOG_PATH, "r", encoding="utf-8", errors="ignore")

            while True:
                line = f.readline() if f else ""
                if line:
                    payload = json.dumps({"type": "log", "text": line.rstrip("\n")})
                    yield f"data: {payload}\n\n"
                else:
                    # periodic heartbeat to defeat buffering
                    now = time.time()
                    if now - heartbeat_last > 1.0:
                        yield ": ping\n\n"
                        heartbeat_last = now

                    if final["done"]:
                        break
                    time.sleep(0.05)
        finally:
            if f:
                f.close()

        # final event
        if not final["ok"]:
            payload = json.dumps({"type": "done", "status": "error", "message": final["error"] or "Unknown error"})
        else:
            payload = json.dumps({
                "type": "done",
                "status": "success",
                "rows": final["rows"] or [],
                "excel_base64": final["excel_b64"]
            })
        yield f"data: {payload}\n\n"

    headers = {
        "Cache-Control": "no-cache, no-transform",
        "X-Accel-Buffering": "no",   # nginx
        "Connection": "keep-alive",
    }
    return StreamingResponse(sse_iter(), media_type="text/event-stream", headers=headers)
