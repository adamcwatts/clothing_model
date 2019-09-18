import numpy as np
import scipy


def absorption(rh):
    # rh is relative humidity array
    # if relative humidity is bound between 0 and 1, can vector-ize function

    number_of_elements = rh.shape[0]
    rh_tuple = (number_of_elements)

    gain = np.ones(rh_tuple)
    vect_func = np.vectorize(regain_function)
    # gain = vect_func(rh)

    for k in range(number_of_elements):

        if rh[k] <= 0:
            gain[k] = 0

        elif rh[k] <= 1:
            gain[k] = regain_function(rh[k])

        elif rh[k] > 1:
            gain[k] = regain_function(1)
    # gain = gain[rh <= 0] * 0
    #
    # gain = gain[rh < 1] * vect_func(rh[rh < 1])
    #
    # gain = gain[rh > 1] * vect_func(rh[rh > 1] == 1)

    return gain


def regain_function(rh):
    return 0.55 * rh * ((1.0 / (0.25 + rh)) + (1.0 / (1.25 - rh)))

# VERIFY ABSORPTION FUNCTION

# relative_h = np.linspace(0, 1, 10)
# print(absorption(relative_h))
