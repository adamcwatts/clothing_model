import numpy as np
import pandas as pd
import MODEL_BC_IC
import math
import scipy.integrate as integrate
from scipy.optimize import fsolve
from time import time

MOLECULAR_WEIGHT_H2O = 18.01528  # [g / mol]
H_VAPORIZATION = 2418  # J / g


def absorption(rh: 'array fraction') -> float:
    # rh is relative humidity array
    # if relative humidity is bound between 0 and 1, can vectorize function

    MAX_REGAIN = 2.6400
    try:
        number_of_elements = rh.shape[0]
    except (AttributeError, IndexError):
        number_of_elements = 1
        rh = np.array([rh])

    gain = np.ones(number_of_elements)

    # g_zero = np.where(rh < 0)
    # g_ones = np.where(rh > 1)

    gain[np.where(rh <= 0)] = 0  # No negative relative humidity
    gain[np.where(rh > 1)] = 2.64  # (PRE-CALCULATED MAX VALUE FROM FUNCTION CALL)

    g_proper = np.where((rh > 0) & (rh < 1))  # collect indices for satisfying this condition
    gain[g_proper] = vectorized_regain(rh[g_proper])

    # when rh <= 0
    # gain[rh <= 0] = gain[rh <= 0] * 0

    # when 0 <= rh <= 1
    # gain[rh <= 1] = gain[rh <= 1] * vectorized_regain(rh[rh <= 1])
    #
    # # when rh > 1
    # gain[rh > 1] = gain[rh > 1] * vectorized_regain(1)

    return gain


def absorption_nv(rh):
    number_of_elements = rh.shape[0]
    gain = np.ones(number_of_elements)

    for k in range(number_of_elements):

        if rh[k] <= 0:
            gain[k] = 0

        elif rh[k] <= 1:
            gain[k] = regain_function(rh[k])

        elif rh[k] > 1:
            gain[k] = regain_function(1)
    return gain


def regain_function(rh: 'fraction') -> float:
    return 0.55 * rh * ((1.0 / (0.25 + rh)) + (1.0 / (1.25 - rh)))


def h_vap_calc(t_celsius: 'celsius', t_kelvin: 'Kelvin') -> 'Joules/gram':
    if t_celsius is not None:  # If Celsius is given
        t_kelvin = t_celsius + 273.15
        enthalpy_vapor = 0.001 * (2.792 * 10 ** 6 - 160 * t_kelvin - 3.43 * t_kelvin ** 2)

    else:  # no celsius, use kelvin
        enthalpy_vapor = 0.001 * (2.792 * 10 ** 6 - 160 * t_kelvin - 3.43 * t_kelvin ** 2)
    return enthalpy_vapor


def h_absorption(rh):
    # Heat of absorption in J/g as a function of RH
    return 195 * (1.0 - rh) * ((1.0 / (0.2 + rh)) + (1.0 / (1.05 - rh)))


def h_sorp_calc(previous_rh, current_rh):
    # total change in heat of Sorption where positive values implies giving off heat due to sorption while
    # negative implies taking in heat (cooling) due to loosing bound water [J/g]
    return integrate.quad(h_absorption, previous_rh, current_rh)[0]  # return value only, not error too


def evap_res_to_diffusion_res(t_celsius: 'celsius', R_ef: 'm ^ 2 Pa / W') -> 's/m':
    # INPUTS:
    # 1) Temperature of Hot Plate when test occurred to get R_ef or any temp for dynamic value
    # 2) Evaporative Resistance [m ^ 2 Pa / W ] which reduces to [s / m]
    #
    # Output: Intrinsic evaporate Diffusion Resistance [s / m] Units: Second per meter

    molecular_weight_h20 = 18.01528  # g / mol
    R_gas_constant = 8.3144598  # J / mol / K
    temp_k = t_celsius + 273.15  # Celsius to Kelvin
    enthalpy_vapor = h_vap_calc(t_celsius, temp_k)  # [J / g]

    diffusion_resistance = R_ef * ((molecular_weight_h20 * enthalpy_vapor) / (R_gas_constant * temp_k))

    return diffusion_resistance


def fabric_parameters(fabric_dictionary) -> 'updated fabric dictionary':
    AIR_SPECIFIC_HEAT_CAPACITY = 1.007  # [J/ kg K]

    fabric_porosity = 1 - fabric_dictionary['porosity of air in fabric']  # estimated porosity of fabric
    p_fab_dry = fabric_dictionary['dry fiber density'] * fabric_porosity  # [kg/m^3]
    fabric_specific_heat_capacity = (fabric_dictionary['fiber specific heat'] * fabric_porosity) + \
                                    (fabric_dictionary['porosity of air in fabric'] * AIR_SPECIFIC_HEAT_CAPACITY)
    diffusion_resistance = evap_res_to_diffusion_res(MODEL_BC_IC.FABRIC_IC_INPUT['initial clothing temp'],
                                                     fabric_dictionary['R_ef'])  # TODO CHECK WHY 35?
    diffusivity_water_though_fabric = fabric_dictionary['fabric thickness'] / diffusion_resistance

    fabric_dictionary['fabric porosity'] = fabric_porosity
    fabric_dictionary['dry fabric density [kg/m^3]'] = p_fab_dry
    fabric_dictionary['dry fabric density [g/m^3]'] = p_fab_dry * 1000
    fabric_dictionary['fabric specific heat capacity [J/ Kg K]'] = fabric_specific_heat_capacity
    fabric_dictionary['diffusivity of water though fabric [m^2 /s]'] = diffusivity_water_though_fabric

    return fabric_dictionary


def fabric_1D_meshing(fabric_dictionary, number_of_nodes,
                      fraction_spacing_of_elements=None) -> 'turns fabric_dictionary into fabric pd data frame':
    # Generates 1D element from inputs, such as thickness of fabric and the spacing of the finite elements
    # 1) Total_Thickness: Total thickness of the fabric in units of meters [m] (data type- double, 1x1 array)
    # 2) Number of elements: How many times to split up fabric e.g. 1mm fabric where 3=number of elements, 1/3mm element
    #   (data type- double, 1x1 array)

    #  OPTIONAL INPUT
    #
    #  3) Fraction_Spacing_of_elements: overrides input 2. Custom spacing fabric
    #   (data type- double, 1xn array)

    # node_name = []
    # node_length = []

    if fraction_spacing_of_elements is None:
        node_length = np.empty(number_of_nodes)  # create empty array
        size_between_nodes = fabric_dictionary['fabric thickness'] / number_of_nodes  # calculate uniform spacing
        node_length.fill(size_between_nodes)  # fill array with said thickness above

    elif fraction_spacing_of_elements is not None:
        tolerance = 10 ** -10
        if not math.isclose(1, sum(fraction_spacing_of_elements), abs_tol=tolerance):
            print(f'ERROR - Total sum of fraction_spacing_of_elements needs to be equal to 1 or less than {tolerance}')
            return

        else:
            node_length = fabric_dictionary['fabric thickness'] * fraction_spacing_of_elements
            number_of_nodes += 1

    # node_names = [f'Fabric Node {x}' for x in range(number_of_nodes)]  # list comprehension of node names
    # data = {'Thickness [m]': node_length}  # combines into dictionary
    # fabric_df = pd.DataFrame(data)  # dictionary to pandas data frame structure
    # fabric_df.index = node_names
    # return fabric_df

    return node_length


def fractional_spacing_generator(number_of_nodes, number_gradient_at_end):
    n = number_gradient_at_end

    if n == 0:
        element_ratio = np.ones(number_of_nodes) / number_of_nodes
        return element_ratio

    new_nodes = (number_of_nodes + 1) + (n - 1) * 2
    element_ratio = np.ones(new_nodes)
    element_ratio[n:-n] = (1 / number_of_nodes)
    if n == 1:
        value = (1 / number_of_nodes) * (1 / 2) ** 1
        element_ratio[0] = value
        element_ratio[-1] = value

    if n >= 2:
        # print((n - i), (1 / number_of_nodes) * (1 / 2) ** (n + 0))
        for i in range(1, n):
            # print((n - i), (-n + (i - 1)))

            element_ratio[n - i] = (1 / number_of_nodes) * (1 / 2) ** (i + 1)
            element_ratio[-n + (i - 1)] = (1 / number_of_nodes) * (1 / 2) ** (i + 1)

        repeat_value = (1 / number_of_nodes) * (1 / 2) ** (n)
        element_ratio[0] = repeat_value
        element_ratio[-1] = repeat_value

    # print('SUM:', sum(element_ratio), '\n')

    return element_ratio


def concentration_calc(t_celsius: 'celsius', rh: 'relative humidity', t_kelvin: 'Kelvin' = None) -> 'g/m^3':
    R = 8.3144598 * 10 ** 3  # [cm^3 kPa K^−1 mol^−1]
    molecular_weight_H2O = 18.01528  # [g / mol]
    # saturated_vp_vectorized = np.vectorize(saturated_vapor_pressure)

    if t_kelvin is None:
        c_mol = (saturated_vapor_pressure(t_celsius) * rh) / (R * (t_celsius + 273.15))  # [mol/cm^3]
        c_mol *= molecular_weight_H2O  # returns [g/cm^3]
        c_mol *= 10 ** 6  # [g/m^3]
    else:
        c_mol = (saturated_vapor_pressure(t_kelvin - 273.15) * rh) / (R * t_kelvin)  # [mol/cm^3]
        c_mol *= molecular_weight_H2O  # returns [g/cm^3]
        c_mol *= 10 ** 6  # [g/m^3]

    return c_mol


def saturated_vapor_pressure(t_celsius: 'celsius', number_of_nodes=None) -> 'kPa':
    #  HUANG - 2018 improved Saturated Vapor Pressure Formula

    if number_of_nodes is None:
        try:
            number_of_nodes = t_celsius.shape[0]
        except (IndexError, AttributeError):
            number_of_nodes = 1
            t_celsius = np.array([t_celsius])

    pressure_saturated = np.ones(number_of_nodes)

    pressure_saturated[t_celsius > 0] = pressure_saturated[t_celsius > 0] * vp_equation_greater_than_freezing(
        t_celsius[t_celsius > 0])

    pressure_saturated[t_celsius <= 0] = pressure_saturated[t_celsius <= 0] * vp_equation_less_than_freezing(
        t_celsius[t_celsius <= 0])

    return pressure_saturated


def sat_vapor_pressure_eq_greater_0(t_celsius: 'celsius') -> 'kPa':
    vapor_pressure = math.exp(34.494 - (4924.99 / (t_celsius + 237.1))) / (t_celsius + 105) ** 1.57  # [kPa]
    vapor_pressure /= 1000  # Takes Pa and converts to kPa
    return vapor_pressure


def sat_vapor_pressure_eq_less_0(t_celsius: 'celsius') -> 'kPa':
    vapor_pressure = math.exp(43.494 - (6545.8 / (t_celsius + 278))) / (t_celsius + 868) ** 2  # [kPa]
    vapor_pressure /= 1000  # Takes Pa and converts to kPa
    vapor_pressure = math.exp(43.494 - (6545.8 / (t_celsius + 278.0))) / (t_celsius + 868) ** 2.0  # [kPa]
    vapor_pressure /= 1000  # Takes Pa and converts to kPa
    return vapor_pressure


# REDUNDANT
def condensation_check(concentration: 'concentration in grams per m^3',
                       t_kelvin: 'Temperature in kelvin') -> 'array of concentration of  water condensate and vapor':
    node_count = concentration.shape[0]
    condensation_concentration = np.zeros(node_count)
    concentration_water_vapor = np.zeros(node_count)

    max_concentration = concentration_calc(None, np.ones(node_count), t_kelvin)

    # mask where concentration in air is greater than saturation point
    condensation_mask = np.where(concentration > max_concentration)
    inverse_mask = np.setdiff1d(np.where(concentration), condensation_mask)  # inverse mask!

    condensation_concentration[condensation_mask] = concentration[condensation_mask] - max_concentration[
        condensation_mask]
    concentration_water_vapor[condensation_mask] = max_concentration[condensation_mask]

    concentration_water_vapor[inverse_mask] = concentration[inverse_mask]

    return concentration_water_vapor, condensation_concentration


def wet_fabric_calc(fabric_df, environmental_rh) -> 'wet_fabric_df':  # TODO Dont use iloc use index name
    #  OUTPUT: multidimensional array where columns refer to node location: TODO: fix this
    #
    # Row (1) refers to density of fabric considering water regain
    # Row (2) refers to specific heat capacity considering water regain
    # Row (3) refers to thermal conductivity considering water regain
    # Row (4) refers to the thickness of the fabric element (could be adjusted
    # for swelling, but isn't)
    #
    # Row(5) Total thickness of the entire fabric layer
    # Row(6) Diffusitvity value

    wet_fabric_properties = {}

    WATER_SPECIFIC_HEAT_CAPACITY = 4179  # [J / (kg K)]
    # T. Bergman and A. Lavine, Fundamentals of Heat and Mass Transfer, 8th ed. Wiley.

    THERMAL_CONDUCTIVITY_WATER = 0.613  # (W/mK) FUNCTION OF TEMP????
    # T. Bergman and A. Lavine, Fundamentals of Heat and Mass Transfer, 8th ed. Wiley.

    FABRIC_THERMAL_CONDUCTIVITY = 0.042  # W/mk TODO volumetric average not Lotens value
    # W. A. Lotens, “Heat transfer from humans wearing clothing,” 1993.

    number_of_nodes = fabric_df.shape[0]
    # wet_fabric_params = np.array((number_of_nodes, 6))
    extracted_data = fabric_df.iloc[0]
    absorption_factor = absorption(environmental_rh.values)

    gamma = (extracted_data.iloc[6] * extracted_data.iloc[3] * absorption_factor) / extracted_data.iloc[6]
    #  (p_fab_dry*Regain*Absorption(Relative_Humidities)) / p_fab_dry
    # fractional density of water in fabric

    wet_fabric_properties['density fraction of water in fabric [kg/m^3]'] = gamma

    wet_fabric_density = extracted_data.iloc[6] + (extracted_data.iloc[6] * extracted_data.iloc[3] * absorption_factor)
    # corrected density including water in fabric

    wet_fabric_properties['wet fabric density [kg/m^3]'] = wet_fabric_density

    wet_fabric_specific_heat = extracted_data.iloc[8] * (1 - gamma) + gamma * WATER_SPECIFIC_HEAT_CAPACITY
    wet_fabric_properties['wet fabric specific heat [J/kg K]'] = wet_fabric_specific_heat

    wet_fabric_thermal_conductivity = FABRIC_THERMAL_CONDUCTIVITY * (1 - gamma) + gamma * THERMAL_CONDUCTIVITY_WATER
    wet_fabric_properties['wet fabric thermal conductivity [W/mk]'] = wet_fabric_thermal_conductivity

    wet_fabric = pd.DataFrame.from_dict(wet_fabric_properties)
    wet_fabric.index = fabric_df.index
    return wet_fabric


def relative_humidity_calc(concentration: 'grams / m^3 H20 in air', temp_kelvin: 'Kelvin') -> 'RH in fraction':
    R = 8.3144598 * 10 ** (-3)  # [m^3 kPa / K mol]
    molar_concentration = concentration / MOLECULAR_WEIGHT_H2O
    vapor_pressure = molar_concentration * R * temp_kelvin
    temperature_celsius = temp_kelvin - 273.15
    relative_humidity = vapor_pressure / saturated_vapor_pressure(temperature_celsius)

    return relative_humidity


def rh_equilibrium(fabric_dataframe, water_vapor_concentration: "grams / m^3 H20 in the air", temperature: 'Kelvin',
                   previous_rh: 'Fabrics previous RH state') -> \
        'RH equilibrium between fabric and air due to sorption process':
    def func_1(x):
        # change_in_concentration of water absorbed water by fabric due to change in RH
        return fabric_dataframe['dry fabric density [g/m^3]'].array[0] * \
               fabric_dataframe['regain'].array[0] * (absorption(x) - absorption(previous_rh))

    def func_2(current_water_vapor_concentration, x):
        # change_in_RH in the air due to the fabric absorbing water from the air
        return relative_humidity_calc(current_water_vapor_concentration - func_1(x), temperature)

    def func_3(x):
        return func_2(water_vapor_concentration, x) - x

    guess = np.ones(water_vapor_concentration.shape[0]) * previous_rh
    # guess = np.array([0.6])
    rh_solution = fsolve(func_3, guess, maxfev=25)
    sorption = func_1(rh_solution)

    # print(rh_solution)
    return rh_solution, sorption


def q_condensation(condensation: "array of condensation [g/m^3] ", fabric_thickness: ' in [m]',
                   time_step: 'in seconds') -> 'Condensation flux in [W/m^2]':
    #  INPUTS:
    #  condensation vector [1, n]  in [g/m^3]
    #  fabric_thickness [1, n] for entire fabric discretization in [m]

    # OUTPUTS:
    # q_condensation array of heat evolved in [J/s m^2] or [W/m^2]

    return (condensation * H_VAPORIZATION * fabric_thickness) / time_step;


def q_evaporation(rh: 'array of relative humidity in fraction', condensation: "array of condensation [g/m^3] ",
                  fabric_thickness: ' in [m]',
                  time_step: 'in seconds') -> 'q_evaporation and concentration array':
    #  INPUTS:
    #  Relative Humidity vector [1, n]  in [fraction]
    #  condensation vector [1, n]  in [g/m^3]
    #  fabric_thickness [1, n] for entire fabric discrete in [m]
    # time step of ODE function [s]

    # OUTPUTS:
    # q_evaporation array of heat evolved in [J/s m^2] or [W/m^2]
    # condensation array corrected to 0 when evaporation takes place

    nodes = rh.shape[0]
    q_evap = np.zeros(nodes)

    evap_mask = np.where((rh < 1.0) & (condensation > 0))

    q_evap[evap_mask] = condensation[evap_mask] * H_VAPORIZATION * fabric_thickness[evap_mask] / time_step
    condensation[evap_mask] = 0

    return q_evap, condensation


def condensation_checker(rh_array: 'relative humidity array', concentration: 'concentration in grams per m^3',
                         temp: 'temperature array in Kelvin',
                         condensation: 'pre_loaded array of zeros') -> 'corrected RH array, updated concentration array,' \
                                                                       ' condensation array':
    rh_mask = rh_array > 1  # where RH > 1
    # condensation = np.zeros(rh_array.shape[0])

    if np.any(rh_mask):  # only run function if at least 1 node has RH > 1
        saturated_vapor_concentration = concentration_calc(None, 1, temp)  # new vapor in air at RH 1
        condensation[rh_mask] = concentration[rh_mask] - saturated_vapor_concentration[
            rh_mask]  # condensation (simulated - saturation)
        concentration[rh_mask] = saturated_vapor_concentration[rh_mask]  # updated air to saturation point
        rh_array[rh_mask] = 1  # update RH to saturation point

    return rh_array, concentration, condensation


def condensation_class(ode_class, index) -> 'corrected RH array, updated concentration array':
    rh_array = ode_class.relative_humidity_post_absorption_history[index, :]
    concentration = ode_class.fiber_water_concentration_history[index, :]
    condensation = ode_class.condensation_concentration_history[index, :]
    temp = ode_class.temps[index - 1:]

    rh_mask = rh_array > 1  # where RH > 1
    # condensation = np.zeros(rh_array.shape[0])

    if np.any(rh_mask):  # only run function if at least 1 node has RH > 1
        saturated_vapor_concentration = concentration_calc(None, 1, temp)  # new vapor in air at RH 1
        condensation[rh_mask] = concentration[rh_mask] - saturated_vapor_concentration[
            rh_mask]  # condensation (simulated - saturation)
        concentration[rh_mask] = saturated_vapor_concentration[rh_mask]  # updated air to saturation point
        rh_array[rh_mask] = 1  # update RH to saturation point


# Vectorization of functions
vectorized_regain = np.vectorize(regain_function,
                                 otypes=[float])  # specifies output is float, works when given empty set
vp_equation_greater_than_freezing = np.vectorize(sat_vapor_pressure_eq_greater_0, otypes=[float])
vp_equation_less_than_freezing = np.vectorize(sat_vapor_pressure_eq_less_0, otypes=[float])
vectorized_h_sorp_cal = np.vectorize(h_sorp_calc, otypes=[float])
# rh_equilibrium = np.vectorize(rh_equilibrium, otypes=[float])
# relative_humidity_calc = np.vectorize(relative_humidity_calc)

if __name__ == '__main__':
    temps = np.array([20, 30, 50])
    temps_kelvin = temps + 273.15
    rh_previous = np.array([0.0, 0, 0])
    rh = np.array([0.0, 0.5, 1.2])
    c = concentration_calc(temps, rh)

    print(condensation_checker(rh, c, temps_kelvin))

    # print(vectorized_h_sorp_cal(rh_previous, rh))

    # print(h_vap_calc.__annotations__)
