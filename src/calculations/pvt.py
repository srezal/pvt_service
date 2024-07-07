from math import log10, e
from .utils import *

def calc_saturated_oil_volumetric_coeff(R_sb: float, GammaGas: float, GammaOil: float, T: float):
    """
    Вычисление объёмного коэффициента насыщенной нефти
    Parameters
    ----------
    :param R_sb: газовый фактор, м3/м3.
    :param GammaGas: отн. плотность газа, доли ед.
    :param GammaOil: отн. плотность нефти, доли ед.
    :param T: пластовая температура, К
    :return: объёмный коэффициент, безразмерный.
    """
    return 0.972 + (147 / 10**6) * (5.61458333333 * R_sb * (GammaGas / GammaOil)**0.5 + 2.25 * T - 574.5875)**1.175


def calc_compressibility(R_sb: float, T: float, GammaGas: float, GammaOil: float, P_bp: float):
    """
    Вычисление коэффициента сжимаемости
    Parameters
    ----------
    :param R_sb: газовый фактор, м3/м3.
    :param T: пластовая температура, F
    :param GammaGas: отн. плотность газа, доли ед.
    :param GammaOil: отн. плотность нефти, доли ед.
    :param P_bp: давление в точке насыщения, Psi
    :return: коэффициент сжимаемости, безразмерный.
    """
    R_sb = R_sb / 0.17810760667903522 # Перевод в куб. фут / баррель. В коде перевод есть, в методичке нет (
    A = (-1.433 + 5 * R_sb + 17.2 * T - 1.180 * GammaGas + 12.61 * gamma_oil_api(GammaOil)) * 10**(-5)
    return A / P_bp

def calc_unsaturated_oil_volumetric_coeff(R_sb: float, GammaGas: float, GammaOil: float, T: float, P):
    """
    Вычисление объёмного коэффициента ненасыщенной нефти
    Parameters
    ----------
    :param R_sb: газовый фактор, м3/м3.
    :param GammaGas: отн. плотность газа, доли ед.
    :param GammaOil: отн. плотность нефти, доли ед.
    :param T: пластовая температура, К
    :param P: пластовое давление, Па
    :return: объёмный коэффициент, безразмерный.
    """
    yg = 1.2254503 + (0.001638 * T) - (1.76875 / GammaOil)
    P_bp = (10**yg / (1.9243101395421235 * 10**(-6))) * (R_sb / GammaGas)**(1/1.2048192771084338) # Давление в точке насыщения, барах
    b_bpp = calc_saturated_oil_volumetric_coeff(R_sb, GammaGas, GammaOil, T)
    T = kelvin_to_fah(T)
    p_res_in_psi = pascal_to_psi(P)
    P_bp_in_psi = pascal_to_psi(bar_to_pascal(P_bp))
    delta_P = P_bp_in_psi - p_res_in_psi
    C = calc_compressibility(R_sb, T, GammaGas, GammaOil, P_bp_in_psi)
    return b_bpp * e**(C*delta_P)


def calc_mu_oil(T: float, GammaOil: float, R_s: float):
    """
    Вычисление вязкости нефти
    Parameters
    ----------
    :param T: пластовая температура, К
    :param GammaOil: отн. плотность нефти, доли ед.
    :param R_s: газосодержание смеси, м3/м3
    :return: вязкость нефти, сПа.
    """
    R_s = R_s / 0.17810760667903522 # Перевод в куб. фут / баррель
    T = kelvin_to_fah(T)
    if T > 70:
        D = T**(-1.163) * 10**(3.0324 - 0.02023 * gamma_oil_api(GammaOil))
        mu_dead = 10**D - 1
    else:
        D_70 = 70**(-1.163) * 10**(3.0324 - 0.02023 * gamma_oil_api(GammaOil))
        D_80 = 80**(-1.163) * 10**(3.0324 - 0.02023 * gamma_oil_api(GammaOil))
        mu_oil_80 = 10**D_80 - 1
        mu_oil_70 = 10**D_70 - 1
        L_7_8 = log10(80 / 70)
        L_mu = log10(mu_oil_70 / mu_oil_80)
        c = L_mu / L_7_8
        b = 70**c * mu_oil_70
        D = log10(b) - c * log10(T)
        mu_dead = 10**D
    mu_live = 10.715 * (R_s + 100)**(-0.515) * mu_dead**(5.44 * (R_s + 150)**(-0.338))
    return mu_live


def calc_pvt(P: float, T: float, GammaOil: float, GammaGas: float,
                GammaWat: float, Wct: float, Rp: float, QLiq: float):
    """
    Вычисление PVT-свойств смеси
    Parameters
    ----------
    :param P: пластовое давление, Па
    :param T: пластовая температура, К
    :param GammaOil: отн. плотность нефти, доли ед.
    :param GammaGas: отн. плотность газа, доли ед.
    :param GammaWat: отн. плотность воды, доли ед.
    :param Wct: обводненность нефти, доли ед.
    :param Rp: газовый фактор, м3/т
    :param QLiq: расход жидкости, м3/cут
    :return: dict: QMix - объём смеси, м3; RhoMix - плотность смеси, кг/м3; MuMix - вязкость смеси, сПа
    """

    # Расчёт поверхностных объёмов
    Rp = Rp / GammaOil # Перевод из м3/т в м3/м3

    V_liq = QLiq / (24 * 3600)
    V_wat = V_liq * Wct
    V_oil = V_liq * (1 - Wct)
    V_gas = V_oil * Rp

    # Расчёт газосодержания нефти
    yg = 1.2254503 + 0.001638 * T - (1.76875 / GammaOil)
    R_s = GammaGas * (1.9243101395421235 * 10**(-6) * P / 10**yg)**1.2048192771084338

    # Расчёт параметров нефти
    if Rp > R_s: # Весь объём газа не сможет раствориться в нефти, так как нефть насыщенная
        b_oil = calc_saturated_oil_volumetric_coeff(R_s, GammaGas, GammaOil, T)
    else: # Нефть недонасыщенная, поэтому весь объём газа сможет раствориться в ней
        b_oil = calc_unsaturated_oil_volumetric_coeff(Rp, GammaGas, GammaOil, T, P)
    V_cond_oil = V_oil * b_oil
    rho_oil = 1000 * (GammaOil + 1.2217 * R_s * GammaGas / 1000) / b_oil
    mu_oil = calc_mu_oil(T, GammaOil, Rp)

    # Расчёт параметров газа
    b_g = 350.958 * T / P
    rho_gas = (28.97 * GammaGas) / (24.04220577350111 * b_g)
    if Rp < R_s:
        V_cond_gas = 0
    else:
        V_cond_gas = b_g * V_oil * (Rp - R_s)
    b = 2.57 + (1914.5 / (1.8 * T)) + 0.275 * GammaGas
    exp = e**(b * (rho_gas / 1000)**(1.11 + 0.04 * b))
    mu_gas = 10**(-4) * (7.77 + 0.183 * GammaGas) * ((1.8 * T)**1.5 / (122.4 + 373.6 * GammaGas + 1.8 * T)) * exp

    # Расчёт параметров жидкости
    wct_cond = V_wat / (V_cond_oil + V_wat)
    V_cond_liq = V_cond_oil + V_wat
    rho_liq = rho_oil * (1 - wct_cond) + 1000 * wct_cond # 1000 - плотность воды
    mu_liq = mu_oil * (1 - wct_cond) + 1 * wct_cond # 1 сПа = вязкость воды

    # Расчёт параметров смеси
    V_cond_mix = V_cond_liq + V_cond_gas
    GF = V_cond_gas / V_cond_mix
    rho_mix = rho_liq * (1 - GF) + rho_gas * GF
    mu_mix = mu_liq * (1 - GF) + mu_gas * GF

    result = {
        "QMix": V_cond_mix,
        "MuMix": mu_mix,
        "RhoMix": rho_mix
    }

    return result

    
    