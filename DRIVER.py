import functions
import MODEL_BC_IC as MODEL
import numpy as np
import pandas as pd
from time import time


def convert_float_to_array(dictionary) -> 'dictionary with numpy array values':
    for key, value in dictionary.items():
        dictionary[key] = np.array([value])
    return dictionary


def fabric_initial_conditions_to_df(fabric_ic: 'temp and relative humidity',
                                    fabric_data_dict: 'fabric data dictionary') -> 'Pandas array':
    # assumes fabric begins in an isothermal and iso-humid state
    temp_dict = {}
    for key, value in fabric_ic.items():
        value = value[0]  # turns 1D array into scalar
        multiplier = np.ones(node_count)  # multiplies values by node count

        if key == 'initial clothing temp':
            temp_dict['initial clothing temp [C]'] = value * multiplier
            temp_dict['initial clothing temp [K]'] = (value + 273.15) * multiplier
        else:
            temp_dict[key] = value * multiplier

    fabric_water_concentration = fabric_data_dict['dry fiber density'] * fabric_data_dict[
        'regain'] * functions.absorption(np.array(fabric_ic['initial clothing rh']))

    temp_dict['water concentration absorbed by fabric [kg/m^3]'] = fabric_water_concentration * multiplier

    water_concentration_in_trapped_air_and_fabric = functions.concentration_calc(fabric_ic['initial clothing temp'],
                                                                                 fabric_ic['initial clothing rh'])
    temp_dict['water concentration in air [g/m^3]'] = water_concentration_in_trapped_air_and_fabric * multiplier

    df = pd.DataFrame.from_dict(temp_dict)
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
    model_fabric_initial_conditions = convert_float_to_array(MODEL.FABRIC_IC_INPUT)
    fabric_df = fabric_initial_conditions_to_df(model_fabric_initial_conditions, fabric_data)
    print(fabric_df)
