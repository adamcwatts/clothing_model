import numpy as np
import scipy


def absorption(rh):
    # rh is relative humidity array
    # if relative humidity is bound between 0 and 1, can vector-ize function

    number_of_elements = rh.shape[0]
    rh_tuple = (number_of_elements)

    gain = np.ones(rh_tuple)
    vect_func = np.vectorize(regain_function)
    # gain = vect_func(rh)

    for k in range(number_of_elements):

        if rh[k] <= 0:
            gain[k] = 0

        elif rh[k] <= 1:
            gain[k] = regain_function(rh[k])

        elif rh[k] > 1:
            gain[k] = regain_function(1)
    # gain = gain[rh <= 0] * 0
    #
    # gain = gain[rh < 1] * vect_func(rh[rh < 1])
    #
    # gain = gain[rh > 1] * vect_func(rh[rh > 1] == 1)

    return gain


def regain_function(rh):
    return 0.55 * rh * ((1.0 / (0.25 + rh)) + (1.0 / (1.25 - rh)))


def h_vap_calc(t_celsius, t_kelvin):
    try:  # If Celsius is given
        t_kelvin = t_celsius + 273.15
        enthalpy_vapor = 0.001 * (2.792 * 10 ^ 6 - 160 * t_kelvin - 3.43 * t_kelvin ** 2)

    except NameError:  # no celsius, use kelvin
        enthalpy_vapor = 0.001 * (2.792 * 10 ^ 6 - 160 * t_kelvin - 3.43 * t_kelvin ** 2)


def evap_res_to_diffusion_res(temp_c, R_ef):
    # INPUTS:
    # 1) Temperature of Hot Plate
    # 2) Evaporative Resistance [m ^ 2 Pa / W ] which reduces to [s / m]
    #
    # Output: Intrinsic evaporative Diffusion Resistance [s / m] Units: Second per meter

    M_W = 18.01528  # g / mol
    R = 8.3144598  # J / mol / K
    temp_k = temp_c + 273.15  # Celsius to Kelvin
    H_Vap = h_vap_calc(temp_c, temp_k)  # [J / g]

    diffusion_resistance = R_ef * ((M_W * H_Vap) / (R * temp_k))

    return diffusion_resistance


def fabric_paramters(fabric_dictionary):
    AIR_SPECIFIC_HEAT_CAPACITY = 1.007  # [J/ kg K]

    fabric_porosity = 1 - AIR_SPECIFIC_HEAT_CAPACITY  # estimated porosity of fabric
    p_fab_dry = fabric_dictionary['dry fiber density'] * fabric_porosity  # [kg/m^3]
    fabric_specific_heat_capacity = (fabric_dictionary['fiber specific heat'] * fabric_porosity) + (
                fabric_dictionary['air porosity'] * AIR_SPECIFIC_HEAT_CAPACITY)
    diffusion_resistence = evap_res_to_diffusion_res(35, fabric_dictionary['R_ef'])  # TODO CHECK WHY 35?
    D_WF = fabric_dictionary['fabric thickness']

    fabric_dictionary['fabric_porosity'] = fabric_porosity
    fabric_dictionary['dry fabric density'] = p_fab_dry
    fabric_dictionary['fabric specific heat capacity'] = fabric_specific_heat_capacity
    fabric_dictionary['D_WF'] = D_WF

    return fabric_dictionary
    # TODO write a return to include either fabric parameters dictionary or array

# VERIFY ABSORPTION FUNCTION

# relative_h = np.linspace(0, 1, 10)
# print(absorption(relative_h))
