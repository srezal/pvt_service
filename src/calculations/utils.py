def pascal_to_psi(pressure_in_pascal: float):
    return pressure_in_pascal * 1.4504 * 10**(-4)


def atm_to_pascal(pressure_in_atm: float):
    return pressure_in_atm * 101325


def bar_to_pascal(pressure_in_bar: float):
    return pressure_in_bar * 10**5 


def celsius_to_kelvin(temp_in_celsius: float):
    return temp_in_celsius + 273.15


def kelvin_to_fah(temp_in_kelvin: float):
    return 1.8 * (temp_in_kelvin - 273.15) + 32
    

def gamma_oil_api(gamma_oil: float):
    return 141.5 / gamma_oil - 131.5