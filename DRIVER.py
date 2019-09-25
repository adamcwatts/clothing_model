import functions
import MODEL_BC_IC as MODEL
import numpy as np
import pandas as pd


def fabric_initial_conditions_to_df(fabric_initial_condition: 'temp and relative humidity') -> 'Pandas array':
    # assumes fabric begins in an isothermal and iso-humid state
    initial_condition = {}
    for key, value in fabric_initial_condition.items():

        if key == 'initial clothing temp':
            initial_condition['initial clothing temp [C]'] = value * np.ones(node_count)
            initial_condition['initial clothing temp [K]'] = (value + 273.15) * np.ones(node_count)
        else:
            initial_condition[key] = value * np.ones(node_count)

    df = pd.DataFrame.from_dict(initial_condition)
    names = [f'Node {x}' for x in range(node_count)]
    df.index = names
    return df


if __name__ == '__main__':
    start = MODEL.TIME_INPUT_PARAMETERS['start']
    finish_seconds = MODEL.TIME_INPUT_PARAMETERS['finish min'] * 60
    time_step = MODEL.TIME_INPUT_PARAMETERS['discrete time step']

    tspan = np.arange(start, finish_seconds + time_step, time_step)
    n = tspan.shape[0]
    node_count = MODEL.NUMBER_OF_NODES

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

    fabric_data = functions.fabric_parameters(MODEL.FABRIC_INPUT_PARAMETERS)
    fabric_dimensions = functions.fabric_1D_meshing(fabric_data, node_count)

    boundary_conditions = MODEL.BOUNDARY_INPUT_PARAMETERS
    initial_conditions = fabric_initial_conditions_to_df(MODEL.IC_INPUT_PARAMETERS)
    print(initial_conditions)
