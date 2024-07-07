from fastapi import APIRouter
from schemas.PvtRequest import PvtRequest
from schemas.PvtResponse import PvtResponse
from calculations.pvt import calc_pvt


router = APIRouter(prefix="/pvt", tags=["PVT"])


@router.post("/calculator")
def calculate_pvt(params_in: PvtRequest) -> PvtResponse:
    data = params_in.model_dump()
    return calc_pvt(**data)
