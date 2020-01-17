import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import finite_dif_pract_1

pi = math.pi


def a_n(n):
    try:
        a = 0
    except ZeroDivisionError:
        a = 0
    return a


def b_n(n):
    b = -70 / (n * pi) * ((-1) ** n - 1)
    return b


def fourier_solution(b_coeff, x_input, interval_length, time=0):
    total = 0
    y_output = np.zeros(x_input.shape[0])
    coefficient_size = len(b_coeff)

    def e_function(n, t, a):
        return math.exp(-a * (n * pi / interval_length) ** 2 * t)

    for counter, x in enumerate(x_input):
        total = 0  # rest total
        # print('counter is:', counter, '\nx is: ', x, )

        for i in range(0, coefficient_size, 2):
            total += b_coeff[i] * math.sin(((i + 1) * x * pi) / interval_length) * e_function(i + 1, time, alpha)

        y_output[counter] = total

    return y_output


vector_a = np.vectorize(a_n)
vector_b = np.vectorize(b_n)

# Material Properties
p = finite_dif_pract_1.p  # density  averaged [kg/m^3]
c = finite_dif_pract_1.c  # specific heat capacity averaged [J/kg k]
k = finite_dif_pract_1.k  # thermal conductivity averaged [w/m K]

alpha = finite_dif_pract_1.alpha

interval_size = finite_dif_pract_1.length / 1000  # meters
a_0 = a_n(0)

x_start = 0
x_end = interval_size
x_domain = np.linspace(x_start, x_end, finite_dif_pract_1.cuts)
coefficients_to_plot = 500
n = np.arange(1, coefficients_to_plot + 1, 1, dtype=np.int32)
a = vector_a(n)  # terms are generated
b = vector_b(n)

y_estimates = fourier_solution(b, x_domain, interval_size, time=0)

fig, ax = plt.subplots(ncols=1, figsize=(14, 9))
ax.set_ylim((0, 40))
ax.set_xlim((x_start, x_end))
ax.plot(x_domain, y_estimates)
plt.xlabel('Length [m]')
plt.ylabel('Temperature [C]')
plt.show()

# TODO write function to add fourier and exponential terms for each n
# compare solutions between Analytical and finite scheme

dt = finite_dif_pract_1.dt
temps_analytic = finite_dif_pract_1.temps
temps_analytic[0, :] = y_estimates

for i in range(1, 10000+1):
    print(i*100/10000)
    temps_analytic[i, :] = fourier_solution(b, x_domain, interval_size, time=dt * i)

# time_array = np.round(dt * np.arange(1000 + 1), decimals=3)
solutions = pd.DataFrame(data=temps_analytic)
solutions.to_csv('analytical_solution.csv')
