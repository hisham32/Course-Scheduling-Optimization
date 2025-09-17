from enum import Enum
from typing import List, Optional
from pydantic import BaseModel, Field


class SemesterEnum(str, Enum):
    Fall = "Fall"
    Spring = "Spring"


class OptimizeRequest(BaseModel):
    plan_semester: SemesterEnum = Field(
        default=SemesterEnum.Fall,
        description="Target semester"
    )
    offerings_path: str = Field(
        default="/data/Hx_Data_Offerings.xlsx",
        description="Path to Hx_Data_Offerings.xlsx (mounted inside container)"
    )
    input_data_path: str = Field(
        default="/data/input_data.xlsx",
        description="Path to input_data.xlsx (mounted inside container)"
    )
    # Optional; no default path
    output_excel_path: Optional[str] = Field(
        default=None,
        description="Where to write the Excel schedule (optional). If omitted, set return_excel_bytes=True to receive it in the response."
    )
    return_excel_bytes: bool = Field(
        default=False,
        description="If true, API returns the generated Excel file as base64 in the response"
    )


class ScheduleRow(BaseModel):
    Course: str
    Day: str
    Start: str
    End: str


class OptimizeResponse(BaseModel):
    status: str
    message: str
    schedule: List[ScheduleRow] = []
    excel_path: Optional[str] = None
    excel_base64: Optional[str] = None
    log: Optional[str] = None
