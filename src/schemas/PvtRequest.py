from pydantic import BaseModel, confloat, Field


class PvtRequest(BaseModel):
    P: float = Field(title="Пластовое давление, Па")
    T: confloat(ge=10, le=500) = Field(title="Пластовая температура, К")
    GammaOil: confloat(ge=0.6, le=1) = Field(title="Отн. плотность нефти, доли ед.")
    GammaGas: confloat(ge=0.5, le=1) = Field(title="Отн. плотность газа, доли ед.")
    GammaWat: confloat(ge=0.98, le=1.2) = Field(title="Отн. плотность воды, доли ед.")
    Wct: confloat(ge=0, le=100) = Field(title="Обводненность, доли ед.")
    Rp: confloat(ge=0) = Field(title="Газовый фактор, м3/т")
    QLiq: float = Field(title="Поверхностный объём жидкости, м3")