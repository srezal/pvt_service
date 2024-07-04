from pydantic import BaseModel, Field


class PvtResponse(BaseModel):
    QMix: float = Field(title="Расход смеси, м3")
    MuMix: float = Field(title="Вязкость смеси, сПа")
    RhoMix: float = Field(title="Плотность смеси, кг/м3")