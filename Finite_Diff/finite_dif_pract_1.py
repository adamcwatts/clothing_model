import numpy as np
import math
import pandas as pd
import scipy


def matrix_a_gen(grid_size, fourier_number):
    n = grid_size
    a = np.zeros([n, n])
    fo = fourier_number
    fo_star = (1 - 2 * fo)

    a[0, 0] = 1
    a[n - 1, n - 1] = 1

    for row in range(grid_size):
        # for column in range(grid_size):
        # if row != 0 and row != grid_size-1:
        if (row != 0 and row != grid_size - 1):
            # A[row, row-1:row+2] = fo
            a[row, row - 1], a[row, row + 1] = fo, fo
            a[row, row] = fo_star
    # print(a)
    return a


def dt_generator(min_dt):
    base = math.ceil(math.fabs(math.log10(min_dt)))
    # new_dt = math.floor(min_dt*10**base) / 10**base
    return 10 ** (-base)


def fd_forward_scheme(current_values, stiffness_matrix):
    return np.matmul(stiffness_matrix, current_values)


length = 10  # mm
cuts = int(100)
dx = length / cuts * (1 / 1000)  # meters

p = (0.07) * (1300) + (1 - 0.07) * (1.225)  # density  averaged [kg/m^3]
c = (0.07) * (1360) + (1 - 0.07) * (1.000)  # specific heat capacity averaged [J/kg k]
k = (0.07) * (0.20) + (1 - 0.07) * (0.024)  # thermal conductivity averaged [w/m K]

alpha = k / (p * c)

min_dt = (dx) ** 2 / (2 * alpha)  # less than this number

dt = dt_generator(min_dt)  # seconds
# dt = 0.1  # seconds

temp_0 = np.zeros(cuts)
temp_0.fill(35)
temp_0[0], temp_0[-1] = 0, 0
# print(np.dot(np.array([1, 2]), np.array([1, 2])))
fo = (dt * alpha) / dx ** 2

A = matrix_a_gen(cuts, fo)

time_minutes = 1
time_seconds = time_minutes * 60
iter_n = int(time_seconds / dt)

temps = np.zeros((iter_n + 1, cuts))
temps[0, :] = temp_0

if __name__ == '__main__':
    for i in range(iter_n):
        temps[i + 1, :] = fd_forward_scheme(temps[i, :], A)

    time_array = np.round(dt * np.arange(iter_n + 1), decimals=3)
    solutions = pd.DataFrame(data=temps, index=time_array)
    solutions.to_csv('finite_solution.csv')
