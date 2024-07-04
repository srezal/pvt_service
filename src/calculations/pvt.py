from math import log10, e

def pascal_to_psi(pressure_in_pascal: float):
    return pressure_in_pascal * 1.4504*10**(-4)


def atm_to_pascal(pressure_in_atm: float):
    return pressure_in_atm * 101325


def bar_to_pascal(pressurer_in_bar: float):
    return pressurer_in_bar * 10**5 


def celsius_to_kelvin(temp_in_celsius: float):
    return temp_in_celsius + 273.15


def kelvin_to_far(temp_in_kelvin: float):
    return 1.8 * (temp_in_kelvin - 273.15) + 32
    

def gamma_oil_api(gamma_oil: float):
    return 141.5 / gamma_oil - 131.5


def calc_saturated_oil_volumetric_coeff(R_sb: float, gamma_gas: float, gamma_oil: float, t_res: float):
    return 0.972 + 147 / 10**6 * (5.614583 * R_sb * (gamma_gas / gamma_oil)**0.5 + 2.25 * t_res - 574.5875)**1.175


def calc_unsaturated_oil_volumetric_coeff(yg: float, R_sb: float, gamma_gas: float, gamma_oil: float, t_res: float, p_res):
    P_bp = (10**yg / (1.9243101395421235 * 10**(-6))) * (R_sb / gamma_gas)**(1/1.2048192771084338) # В барах
    p_res = pascal_to_psi(p_res)
    P_bp = pascal_to_psi(bar_to_pascal(P_bp))
    b_bpp = calc_saturated_oil_volumetric_coeff(R_sb, gamma_gas, gamma_oil, t_res)
    t_res = kelvin_to_far(t_res)
    A = (-1.433 + 5 * R_sb + 17.2 * t_res - 1.180 * gamma_gas + 12.61 * gamma_oil_api(gamma_oil)) * 10**(-5)
    C = A / P_bp
    return b_bpp * e**(C*(P_bp - p_res)) # Возможно нужен перевод в другие ед.


def calc_mu_oil(t_res: float, gamma_oil: float, R_s: float): # Здесь вопрос
    R_s = R_s / 0.17810760667903522
    t_res = kelvin_to_far(t_res)
    if t_res > 70:
        D = t_res**(-1.163) * 10**(3.0324 - 0.02023 * gamma_oil_api(gamma_oil))
        mu_dead = 10**D - 1
        mu_live = 10.715 * (R_s + 100)**(-0.515) * mu_dead**(5.44 * (R_s + 150)**(-0.338))
        return mu_live
    else:
        D_70 = 70**(-1.163) * 10**(3.0324 - 0.02023 * gamma_oil_api(gamma_oil))
        D_80 = 80**(-1.163) * 10**(3.0324 - 0.02023 * gamma_oil_api(gamma_oil))
        mu_oil_80 = 10**D_80 - 1
        mu_oil_70 = 10**D_70 - 1
        L_7_8 = log10(80 / 70)
        L_mu = log10(mu_oil_70 / mu_oil_80)
        c = L_mu / L_7_8
        b = 70**c * mu_oil_70
        D = log10(b) - c * log10(t_res)
        mu_oil = 10**D
        return mu_oil


def calc_pvt(p_res: float, t_res: float, gamma_oil: float, gamma_gas: float,
                gamma_wat: float, wct: float, rp: float, q_liq: float):

    p_res = atm_to_pascal(p_res)
    t_res = celsius_to_kelvin(t_res)
    wct = wct / 100

    # Расчёт поверхностных объёмов
    V_liq = q_liq
    V_wat = V_liq * wct
    V_oil = V_liq * (1 - wct)
    V_gas = V_oil * rp

    # Расчёт газосодержания нефти
    yg = 1.2254503 + 0.001638 * t_res - 1.76875 / gamma_oil
    R_s = gamma_gas * (1.9243101395421235 * 10**(-6) * p_res / 10**yg)**1.2048192771084338

    # Расчёт параметров нефти
    if rp < R_s: # Весь объём газа растворится в нефти
        b_oil = calc_saturated_oil_volumetric_coeff(rp, gamma_gas, gamma_oil, t_res)
    else:
        b_oil = calc_unsaturated_oil_volumetric_coeff(yg, rp, gamma_gas, gamma_oil, t_res, p_res)
    V_cond_oil = V_oil * b_oil
    rho_oil = 1000 * (gamma_oil + 1.2217 * R_s * gamma_gas / 1000) / b_oil
    mu_oil = calc_mu_oil(t_res, gamma_oil, R_s)
    #print(f'V_cond_oil = {V_cond_oil}\nrho_oil = {rho_oil}\nmu_oil = {mu_oil}')

    # Расчёт параметров газа
    b_g = 350.958 * t_res / p_res
    rho_gas = (28.97 * gamma_gas) / (24.04220577350111 * b_g)
    if rp < R_s:
        V_cond_gas = 0
    else:
        V_cond_gas = b_g * V_gas * (rp - R_s)
    b = 2.57 + 1914.5 / (1.8 * t_res) + 0.275 * gamma_gas
    exp = e**(b * (rho_gas / 1000)**(1.11 + 0.04 * b))
    mu_gas = 10**(-4) * (7.77 + 0.183 * gamma_gas) * (1.8 * t_res)**1.5 / (122.4 + 373.6 * gamma_gas + 1.8 * t_res) * exp
    #print(f'V_cond_gas = {V_cond_gas}\nrho_gas = {rho_gas}\nmu_gas = {mu_gas}')


    # Расчёт параметров жидкости
    V_cond_liq = V_cond_oil + V_wat
    rho_liq = rho_oil * (1 - wct) + 1000 * wct # 1000 - плотность воды
    mu_liq = mu_oil * (1 - wct) + 1 * wct # 1 сПа = вязкость воды
    #print(f'V_cond_liq = {V_cond_liq}\nrho_liq = {rho_liq}\nmu_liq = {mu_liq}')

    # Расчёт параметров смеси
    V_cond_mix = V_cond_liq + V_cond_gas
    GF = V_cond_gas / (V_cond_liq + V_cond_mix)
    rho_mix = rho_liq * (1 - GF) + rho_gas * GF
    mu_mix = mu_liq * (1 - GF) + mu_gas * GF

    result = {
        "QMix": V_cond_mix,
        "MuMix": mu_mix,
        "RhoMix": rho_mix
    }

    return result

    
    