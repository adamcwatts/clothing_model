import numpy as np
import pandas as pd
import scipy
import MODEL_BC_IC
import math
from time import time


def absorption(rh: 'array fraction') -> float:
    # rh is relative humidity array
    # if relative humidity is bound between 0 and 1, can vectorize function

    number_of_elements = rh.shape[0]
    gain = np.ones(number_of_elements)

    # vect_func = np.vectorize(regain_function, otypes=[float])  # specifies output is float, works when given empty set

    # when rh <= 0
    gain[rh <= 0] = gain[rh <= 0] * 0

    # when 0 <= rh <= 1
    gain[rh <= 1] = gain[rh <= 1] * vectorized_regain(rh[rh <= 1])

    # when rh > 1
    gain[rh > 1] = gain[rh > 1] * vectorized_regain(1)

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


def evap_res_to_diffusion_res(t_celsius: 'celsius', R_ef: 'm ^ 2 Pa / W') -> 's/m':
    # INPUTS:
    # 1) Temperature of Hot Plate when test occurred to get R_ef
    # 2) Evaporative Resistance [m ^ 2 Pa / W ] which reduces to [s / m]
    #
    # Output: Intrinsic evaporate Diffusion Resistance [s / m] Units: Second per meter

    molecular_weight_h20 = 18.01528  # g / mol
    R_gas_constant = 8.3144598  # J / mol / K
    temp_k = t_celsius + 273.15  # Celsius to Kelvin
    enthalpy_vapor = h_vap_calc(t_celsius, temp_k)  # [J / g]

    diffusion_resistance = R_ef * ((molecular_weight_h20 * enthalpy_vapor) / (R_gas_constant * temp_k))

    return diffusion_resistance


def fabric_parameters(fabric_dictionary):
    AIR_SPECIFIC_HEAT_CAPACITY = 1.007  # [J/ kg K]

    fabric_porosity = 1 - fabric_dictionary['% porosity of air in fabric']  # estimated porosity of fabric
    p_fab_dry = fabric_dictionary['dry fiber density'] * fabric_porosity  # [kg/m^3]
    fabric_specific_heat_capacity = (fabric_dictionary['fiber specific heat'] * fabric_porosity) + \
                                    (fabric_dictionary['% porosity of air in fabric'] * AIR_SPECIFIC_HEAT_CAPACITY)
    diffusion_resistance = evap_res_to_diffusion_res(35, fabric_dictionary['R_ef'])  # TODO CHECK WHY 35?
    diffusivity_water_though_fabric = fabric_dictionary['fabric thickness'] / diffusion_resistance

    fabric_dictionary['fabric_porosity'] = fabric_porosity
    fabric_dictionary['dry fiber density'] = p_fab_dry
    fabric_dictionary['fabric specific heat capacity'] = fabric_specific_heat_capacity
    fabric_dictionary['diffusivity of water though fabric'] = diffusivity_water_though_fabric

    return fabric_dictionary


def fabric_1D_meshing(fabric_dictionary, number_of_nodes, fraction_spacing_of_elements=None):
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

    node_names = [f'fabric element {x}' for x in range(number_of_nodes)]  # list comprehension of node names
    data = {'Element': node_names, 'Length': node_length}  # combines into dictionary
    fabric_df = pd.DataFrame(data)  # dictionary to pandas data frame structure
    return fabric_df


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
        c_mol = (saturated_vapor_pressure(t_kelvin) * rh) / (R * t_kelvin)  # [mol/cm^3]
        c_mol *= molecular_weight_H2O  # returns [g/cm^3]
        c_mol *= 10 ** 6  # [g/m^3]

    return c_mol


def saturated_vapor_pressure(t_celsius: 'celsius', number_of_nodes=None) -> 'kPa':
    if number_of_nodes is None:
        number_of_nodes = t_celsius.shape[0]

    pressure_saturated = np.ones(number_of_nodes)
    vp_equation_greater_than_freezing = np.vectorize(sat_vapor_pressure_eq_greater_0, otypes=[float])
    vp_equation_less_than_freezing = np.vectorize(sat_vapor_pressure_eq_less_0, otypes=[float])

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


if __name__ == '__main__':
    vectorized_regain = np.vectorize(regain_function,
                                     otypes=[float])  # specifies output is float, works when given empty set

    # print(fractional_spacing_generator(4, 2))

    a = np.linspace(0, 1, 11)
    temp = np.linspace(20, 40, 11)

    iterations = 1000000

    start_time = time()
    for i in range(iterations):
        absorption(a)
    print('abs_function vectorized Elapsed time:' + str(time() - start_time))

    start_time = time()
    for i in range(iterations):
        absorption_nv(a)
    print('abs_function non-vectorized  Elapsed time:' + str(time() - start_time))

    # print(saturated_vapor_pressure(temp))

    # print(concentration_calc(temp, a))
    # print(concentration_calc.__annotations__)

    # print(absorption(a * 2))
    # print(h_vap_calc.__annotations__)
    # element_fraction = fractional_spacing_generator(MODEL_BC_IC.NUMBER_OF_NODES)
    # print(fabric_1D_meshing(MODEL_BC_IC.FABRIC_INPUT_PARAMETERS, MODEL_BC_IC.NUMBER_OF_NODES, element_fraction))
