
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from .models import OptimizeRequest, OptimizeResponse, ScheduleRow
from .scheduler import run_and_optionally_write_excel
import base64

app = FastAPI(title="CEE Scheduler Microservice", version="1.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/health")
def health():
    return {"status": "ok"}

@app.post("/optimize", response_model=OptimizeResponse)
def optimize(req: OptimizeRequest):
    try:
        df_final, excel_bytes = run_and_optionally_write_excel(
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
            excel_path=req.output_excel_path
        )
        if req.return_excel_bytes and excel_bytes is not None:
            resp.excel_base64 = base64.b64encode(excel_bytes).decode("utf-8")
        return resp
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
