import functions
import MODEL_BC_IC
import numpy as np


def initialize_arrays():
    start = MODEL_BC_IC.TIME_INPUT_PARAMETERS['start']
    finish_seconds = MODEL_BC_IC.TIME_INPUT_PARAMETERS['finish min'] * 60
    time_step = MODEL_BC_IC.TIME_INPUT_PARAMETERS['discrete time step']

    tspan = np.arange(start, finish_seconds + time_step, time_step)
    n = tspan.shape[0]
    node_count = MODEL_BC_IC.NUMBER_OF_NODES
    temps = np.zeros((n, node_count))

    if __name__ == '__main__':
        input()
