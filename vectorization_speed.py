import numpy as np
from time import time
from functions import absorption
from functions import absorption_nv

if __name__ == '__main__':
    array_sizes = np.logspace(1, 5, 5)
    time_data_no_vector = []
    time_data_vector = []

    for array in array_sizes:
        a = np.linspace(0, 1, array)
        # temp = np.linspace(20, 40, 1100)

        iterations = 1000

        start_time = time()
        for i in range(iterations):
            absorption(a)
        time_end = (time() - start_time)
        time_data_vector.append(time_end)
        print('abs_function vectorized Elapsed time:' + str(time() - start_time))

        start_time = time()
        for i in range(iterations):
            absorption_nv(a)
        time_end = (time() - start_time)
        time_data_no_vector.append(time_end)
        print('abs_function non-vectorized  Elapsed time:' + str(time() - start_time))

    print(array_sizes)
    print(time_data_vector, time_data_no_vector)
