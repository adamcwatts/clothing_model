import functions
import MODEL_BC_IC as BC_IC  # BC_IC - Boundary and Initial Conditions
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
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
    #  y[n] refers to temperature at nth node : 1 <= n <= number of nodes

    # ode_params = {}
    # n numbers of nodes + 1 node for boundary at plate

    temp_at_clothing_boundary = y[node_count]  # last element of vector associated with temperatures

    h_radiation = OPI['eps_clothing'] * OPI['sigma'] * (
            temp_at_clothing_boundary + boundary_conditions['air temp [K]']) * (
                          temp_at_clothing_boundary ** 2 + boundary_conditions['air temp [K]'] ** 2)

    # r_out = 1 / (h_radiation + OPI['h_convection'])    U_out = 1 / r_out

    r_rad_convect = 1 / (h_radiation + OPI['h_convection'])

    delta_x = df_fabric_initial_conditions['thickness [m]'].to_numpy()

    K = df_wet_fabric['wet fabric thermal conductivity [W/mk]'].to_numpy()
    C = df_wet_fabric['wet fabric specific heat [J/kg K]'].to_numpy()
    P = df_wet_fabric['wet fabric density [kg/m^3]'].to_numpy()
    U = K / delta_x
    U[-1] = 1 / (r_rad_convect + 0.5 * (
            delta_x[-1] / K[-1]))  # corrects last element for radiation, convection heat loss

    material = 1 / (delta_x * np.multiply(C, P))

    temperature_array = np.append(y[1:], boundary_conditions['air temp [K]'])

    q_dry = np.diff(temperature_array) * -1 * U  # -1 due to diff goes a[i+1] - a[i], I need a[i] - a[i+1]
    q_dry = np.concatenate(([y[0]], q_dry), axis=0)  # adds heat flow from plate

    q_gen = (boundary_conditions['plate temp [K]'] - y[1]) * (K[0] / (delta_x[0] * 0.5))  # flux at plate ODE
    clothing_dydt = material * np.diff(q_dry) * -1  # -1 due to diff goes a[i+1] - a[i], I need a[i] - a[i+1]

    dydt = np.concatenate(([q_gen], clothing_dydt), axis=0)

    # q_1_2 = (y[1] - y[2]) * U[0]
    # q_2_3 = (y[2] - y[3]) * U[1]
    # q_3_out = (y[3] - boundary_conditions['air temp [K]']) * U[2]  #
    #
    # dydt_1 = material[0] * (y[0] - q_1_2)  # delta temp at node ODE 1
    # dydt_2 = material[1] * (q_1_2 - q_2_3) # delta temp at node ODE 2
    # dydt3 = material[2] * (q_2_3 - q_3_out) # delta temp at node ODE 3

    # dydt = [q_gen, dydt_1, dydt_2, dydt3]

    return dydt


def solution_to_df(ode_sol):
    # takes ODE solution and creates data frame that's easy to read solutions off of

    fabric_solutions = {f'Temp at Node {i + 1} [C]': ode_sol.y[i + 1] - 273.15 for i in range(node_count)}
    fabric_solutions['Time (Minutes)'] = ode_sol.t / 60
    fabric_solutions['Hot Plate Flux'] = ode_sol.y[0]

    index_names = ode_sol.t
    index_names = np.ndarray.tolist(ode_sol.t)
    index_names = [str(name) + " [s]" for name in index_names]

    # reorder so time is at the beginning of the Data Frame
    sol_df = pd.DataFrame.from_dict(fabric_solutions)
    sol_df = sol_df[sol_df.columns.tolist()[-2:] + sol_df.columns.tolist()[:-2]]  # reorder last 2 as first 2
    sol_df.index = index_names  # add index names

    return sol_df


def solution_plots(solution_df):
    middle_node = (node_count + 1) / 2

    y_values = solution_df[f'Temp at Node {int(middle_node)} [C]'].to_numpy()
    x_values = solution_df['Time (Minutes)'].to_numpy()

    y_max = y_values.max()
    y_min = y_values.min()

    fig, ax = plt.subplots(ncols=1, figsize=(14, 7))

    ax.set_ylim((y_min - 1, y_max + 1))
    ax.set_xlim((0, x_values.max()))

    ax.plot(x_values, y_values, label=f'Temperature at Node {middle_node} [C]')
    ax.set_xlabel('Time [Minutes]')
    ax.set_ylabel(r'Temperature $[\degreeC]$')
    plt.savefig('Middle_Fabric_Node_Temp_C.pdf', bbox_inches="tight")
    plt.show()


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

    # could condense this into a 3D array time X nodes X type of heat flow
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

    # TEST RH Equilibrium Function

    # temps_c = np.array([20, 30, 35, 40])
    # temps_k = temps_c + 273.15
    #
    # incoming_air_rh = np.array([1, 1, 1, 1]) * 0.75
    # prior_rh = np.array([1, 1, 1, 1]) * 0.25
    #
    # concent_grams = functions.concentration_calc(temps_c, incoming_air_rh)
    #
    # sol_temp = functions.rh_equilibrium(df_fabric_initial_conditions, concent_grams, temps_k, prior_rh)
    # print(sol_temp)
    # print(functions.rh_equilibrium.__annotations__)

    # ODE STUFF BELOW

    tspan = np.linspace(start, finish_seconds, finish_seconds + 1)

    # attach 0 at the 0 position to add initial plate heat flux
    yinit = np.insert(df_fabric_initial_conditions['initial clothing temp [K]'].to_numpy(), 0, 0)

    solution = solve_ivp(odefun, [tspan[0], tspan[-1]], yinit, t_eval=tspan)
    print(solution.y)

    sol_df = solution_to_df(solution)
    export_csv = sol_df.to_csv('solution_summary.csv', sep='\t', encoding='utf-8')

    # solution_plots(sol_df)
