import functions
import INPUT
import numpy as np


def initialize_arrays():
    start = INPUT.TIME_INPUT_PARAMETERS['start']
    finish_seconds = INPUT.TIME_INPUT_PARAMETERS['finish_min'] * 60
    time_step = INPUT.TIME_INPUT_PARAMETERS['discrete time step']

    tspan = np.arange(start, finish_seconds + time_step, time_step)
    n = tspan.shape[0]
    node_count = INPUT.NUMBER_OF_NODES
    temps = np.zeros((n, node_count))

    if __name__ == '__main__':
        input()
