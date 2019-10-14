import functions
import MODEL_BC_IC as BC_IC  # BC_IC - Boundary and Initial Conditions
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from time import time

OPI = BC_IC.ODE_PHYSICS_INPUT


def convert_float_to_array(dictionary) -> 'dictionary with numpy array values':
    for key, value in dictionary.items():
        dictionary[key] = np.array([value])
    return dictionary


def fabric_initial_conditions_to_df(fabric_ic: 'temp and relative humidity',
                                    fabric_data_dict: 'fabric data dictionary',
                                    node_thickness: 'array of thickness of nodes') -> 'Pandas array':
    # assumes fabric begins in an isothermal and iso-humid state
    temp_dict = {}
    for key, value in fabric_ic.items():
        # value = value[0]  # turns 1D array into scalar :TODO might not need this
        multiplier = np.ones(node_count)  # multiplies values by node count

        if key == 'initial clothing temp':
            temp_dict['initial clothing temp [C]'] = value * multiplier
            temp_dict['initial clothing temp [K]'] = (value + 273.15) * multiplier
        else:
            temp_dict[key] = value * multiplier

    for key, value in fabric_data_dict.items():
        if ('fiber' or 'thickness') not in key:
            temp_dict[key] = value * multiplier

    temp_dict['R_ef [m^2 Pa / W]'] = temp_dict.pop('R_ef')  # adds units
    temp_dict['total fabric thickness [m]'] = temp_dict.pop('fabric thickness')

    fabric_water_concentration = fabric_data_dict['dry fiber density'] * fabric_data_dict[
        'regain'] * functions.absorption(np.array(fabric_ic['initial clothing rh']))

    temp_dict['water concentration absorbed by fabric [kg/m^3]'] = fabric_water_concentration * multiplier

    water_concentration_in_trapped_air_and_fabric = functions.concentration_calc(fabric_ic['initial clothing temp'],
                                                                                 fabric_ic['initial clothing rh'])
    temp_dict['water concentration in air [g/m^3]'] = water_concentration_in_trapped_air_and_fabric * multiplier

    temp_dict['thickness [m]'] = node_thickness

    df = pd.DataFrame.from_dict(temp_dict)
    names = [f'Fabric Node {x}' for x in range(node_count)]
    df.index = names
    # df.rename(columns={'initial clothing rh': 'Initial clothing rh'}, inplace=True)
    return df


def boundary_conditions_with_kelvin(boundary_condition):
    updated_boundary_conditions = boundary_condition.copy()
    for key, value in boundary_condition.items():
        if 'temp' in key:
            new_key_celsius = key + ' [C]'
            new_key_kelvin = key + ' [K]'
            updated_boundary_conditions[new_key_celsius] = updated_boundary_conditions.pop(key)  # update key
            updated_boundary_conditions[new_key_kelvin] = value + 273.15  # kelvin

    # df = pd.DataFrame.from_dict(boundary_condition)
    return updated_boundary_conditions


def odefun(t, y):
    #  y[0] refers to heat flux needed by plate
    #  y[1] refers to temperature at middle of clothing node

    ode_params = {}
    # 6.5mm effective air layer thickness between outer textile and ambient air
    # n numbers of nodes + 1 node for boundary at plate

    # elements = z.shape[0]
    #
    # temp_at_plate = z[0]

    # temp_at_clothing_boundary = y[1]
    #
    # h_radiation = OPI['eps_clothing'] * OPI['sigma'] * (
    #         temp_at_clothing_boundary + boundary_conditions['air temp [K]']) * (
    #                       temp_at_clothing_boundary ** 2 + boundary_conditions['air temp [K]'] ** 2)

    delta_x = df_fabric_initial_conditions['thickness [m]'].to_numpy()

    K = df_wet_fabric['wet fabric thermal conductivity [W/mk]'].to_numpy()
    C = df_wet_fabric['wet fabric specific heat [J/kg K]'].to_numpy()
    P = df_wet_fabric['wet fabric density [kg/m^3]'].to_numpy()
    U = K / delta_x

    material = 1 / (delta_x * np.multiply(C, P))

    q_1_2 = (y[1] - boundary_conditions['air temp [K]']) * U[0]

    # dydt = np.zeros(2)

    q_gen = (boundary_conditions['plate temp [K]'] - y[1]) * (K[0] / (delta_x[0] * 0.5))  # flux at plate ODE

    dydt_1 = material[0] * (y[0] - q_1_2)  # temp at node ODE

    dydt = [q_gen, dydt_1]
    return dydt


if __name__ == '__main__':
    start = BC_IC.TIME_INPUT_PARAMETERS['start']
    finish_seconds = BC_IC.TIME_INPUT_PARAMETERS['finish min'] * 60
    time_step = BC_IC.TIME_INPUT_PARAMETERS['discrete time step']

    tspan = np.arange(start, finish_seconds + time_step, time_step)
    n = tspan.shape[0]
    node_count = BC_IC.NUMBER_OF_NODES

    node_tuple = (n, node_count)
    zeros_array = np.zeros(node_tuple)

    temps = np.zeros(node_tuple)  # [kelvin]
    raw_water_concentration = zeros_array  # [(g/m^3]  NEED TO VERIFY
    condensation_concentration_history = zeros_array  # [(g/m^3]  NEED TO VERIFY
    fiber_water_concentration_history = zeros_array  # [(g/m^3]  NEED TO VERIFY

    enthalpy_sorption_history = zeros_array

    q_sorption = zeros_array
    q_condensation = zeros_array

    q_evaporation = zeros_array
    q_evaporation_history = zeros_array

    q_1_history = zeros_array
    q_L_history = zeros_array

    relative_humidity_post_absorption_history = zeros_array

    final_temps_and_concentration = zeros_array

    heat_flows_history = np.zeros((n, 1 + node_count))

    fabric_data = functions.fabric_parameters(BC_IC.FABRIC_INPUT_PARAMETERS)
    fabric_node_dimensions = functions.fabric_1D_meshing(fabric_data, node_count)

    boundary_conditions = boundary_conditions_with_kelvin(BC_IC.BOUNDARY_INPUT_PARAMETERS)
    fix_floats_to_numpy_array = convert_float_to_array(BC_IC.FABRIC_IC_INPUT)

    df_fabric_initial_conditions = fabric_initial_conditions_to_df(fix_floats_to_numpy_array, fabric_data,
                                                                   fabric_node_dimensions)

    df_wet_fabric = functions.wet_fabric_calc(df_fabric_initial_conditions,
                                              df_fabric_initial_conditions['initial clothing rh'])

    tspan = np.linspace(0, 10, 1100)
    yinit = [0, 306.85]

    sol = solve_ivp(odefun, [tspan[0], tspan[-1]], yinit, t_eval=tspan)
    print(sol.y)
    sol_df = pd.DataFrame.from_dict({'Time': sol.t, 'Plate Heat Flux': sol.y[0], 'Temp at Fabric': sol.y[1]})
