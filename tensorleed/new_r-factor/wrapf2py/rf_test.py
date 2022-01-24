import numpy as np
from rfactor import r_factor_new as rf
from numpy import asfortranarray
import time
import timeit  # dependency for testing! Don't use in Viperleed build
from matplotlib import pyplot as plt
import pickle
import random

"""
For compiling:
bash test.sh

make sure file .f2py_f2cmap exists!
"""

"""
Signature:
subroutine prepare_beams(n_beams, n_E_in, E_grid_in, intensities, E_start_beams, n_E_beams, &
                         skip_stages, n_averaging_types, averaging_scheme, n_beams_final,&
                         deg, &
                         n_E_out, E_grid_out, &
                         E_start_beams_out, n_E_beams_out, &
                         intpol_intensity, intpol_derivative, y_func)
"""

# Read Hematite test data from file EXPBEAMS.csv
fname = "EXPBEAMS.csv"

data = np.asfortranarray(np.genfromtxt(fname, delimiter=",", skip_header=1))
while np.all(np.isnan(data[:, -1])):
    data = data[:, :-2]  # remove tailing colums with only NaNs

# remove last col for testing...
data = data[:, :-1]

E_grid_in = data[:, 0]
intensities = data[:, 1:]

n_E_in, n_beams = intensities.shape

check_not_nan = np.invert(np.isnan(intensities))

lower_id = np.int32(np.zeros([n_beams]))
upper_id = np.int32(np.zeros([n_beams]))
for i in range(n_beams):
    lower_id[i] = np.argwhere(check_not_nan[:, i])[0]
    upper_id[i] = np.argwhere(check_not_nan[:, i])[-1]

E_start_beams = lower_id + 1  # Fortran index
n_E_beams = upper_id - lower_id

skip_stages = np.int32([0, 0, 0, 0, 0])

n_beams_out = n_beams
averaging_scheme = np.int32(np.arange(n_beams) + 1)  # Fortran index

"""
n_beams_out -= 5
averaging_scheme[1] = 1
averaging_scheme[2] = 1
averaging_scheme[3] = 2
averaging_scheme[4] = 2
averaging_scheme[5] = 2
averaging_scheme[6:] -= 5
"""

deg = 5

n_E_out = n_E_in * 5
E_grid_out = np.linspace(75, 750, num=n_E_out)

(
    e_start_beams_out,
    n_e_beams_out,
    intpol_intensity,
    intpol_derivative,
    y_func,
) = rf.prepare_beams(
    E_grid_in,
    intensities,
    E_start_beams,
    n_E_beams,
    skip_stages,
    n_beams_out,
    averaging_scheme,
    deg,
    E_grid_out,
)


res = (
    n_beams_out,
    e_start_beams_out,
    n_e_beams_out,
    intpol_intensity,
    intpol_derivative,
    y_func,
)

# e_start_beams_out -= 1

print("n_beams_out", n_beams_out)
print("e_start_beams_out", e_start_beams_out)
print("n_e_beams_out", n_e_beams_out)


# plt.clf()
for i in range(n_beams_out):
    x = E_grid_in[E_start_beams[i] : E_start_beams[i] + n_E_beams[i]]
    y = intensities[:, i][E_start_beams[i] : E_start_beams[i] + n_E_beams[i]]
    # plt.plot(x, y, ls= "", marker = ".", markersize= 5, label = "raw")

    x = E_grid_out[e_start_beams_out[i] : e_start_beams_out[i] + n_e_beams_out[i]]
    # plt.plot(x, intpol_intensity[:,i][e_start_beams_out[i]:e_start_beams_out[i]+n_e_beams_out[i]], label = "intpol_intensity")
    # plt.plot(x, intpol_derivative[:,i][e_start_beams_out[i]:e_start_beams_out[i]+n_e_beams_out[i]], label = "intpol_derivative")
    # plt.plot(
    #    x,
    #    y_func[:, i][e_start_beams_out[i] : e_start_beams_out[i] + n_e_beams_out[i]],
    #    label="y_func",
    # )
# plt.legend()
# plt.show()

print(intpol_intensity.shape)

np.savetxt("debug_inten.csv", intpol_intensity)
np.savetxt("debug_deriv.csv", intpol_derivative)
np.savetxt("debug_y_func.csv", y_func)


e_step = E_grid_out[1] - E_grid_out[0]

a = rf.r_pendry_beam_y(
    e_step,
    y_func[:, 0],
    y_func[:, 0],
    e_start_beams_out[0],
    e_start_beams_out[0],
    n_e_beams_out[0],
    n_e_beams_out[0],
    5.0,
)
b = rf.r_pendry_beam_y(
    e_step,
    y_func[:, 0],
    y_func[:, 1],
    e_start_beams_out[0],
    e_start_beams_out[1],
    n_e_beams_out[0],
    n_e_beams_out[1],
    5.0,
)

c = rf.r_pendry_beamset_y(
    e_step,
    y_func,
    y_func,
    e_start_beams_out,
    e_start_beams_out,
    n_e_beams_out,
    n_e_beams_out,
    5.0,
)

d = rf.r_pendry_beamset_y(
    e_step,
    y_func[:,1:],
    y_func[:,:-1],
    e_start_beams_out[1:],
    e_start_beams_out[:-1],
    n_e_beams_out[:-1],
    n_e_beams_out[:-1],
    5.0,
)