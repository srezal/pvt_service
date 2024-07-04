from pydantic import BaseModel, confloat, Field


class PvtRequest(BaseModel):
    p_res: float = Field(title="Пластовое давление, атм")
    t_res: confloat(ge=10, le=500) = Field(title="Пластовая температура, C")
    gamma_oil: confloat(ge=0.6, le=1) = Field(title="Отн. плотность нефти")
    gamma_gas: confloat(ge=0.5, le=1) = Field(title="Отн. плотность газа")
    gamma_wat: confloat(ge=0.98, le=1.2) = Field(title="Отн. плотность воды")
    wct: confloat(ge=0, le=100) = Field(title="Обводненность, %")
    rp: confloat(ge=0) = Field(title="Газовый фактор, м3/т")