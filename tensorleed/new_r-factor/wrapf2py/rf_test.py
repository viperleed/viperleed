import numpy as np
from rfactor import r_factor_new as rf
from numpy import asfortranarray
import time
import timeit # dependency for testing! Don't use in Viperleed build
from matplotlib import pyplot as plt

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
while np.all(np.isnan(data[:,-1])):
    data = data[:,:-2] # remove tailing colums with only NaNs



E_grid_in = data[:,0]
intensities = data[:,1:]

n_E_in, n_beams = intensities.shape

check_not_nan = np.invert(np.isnan(intensities))

lower_id = np.int32(np.zeros([n_beams]))
upper_id = np.int32(np.zeros([n_beams]))
for i in range(n_beams):
    lower_id[i] = np.argwhere(check_not_nan[:,i])[0]
    upper_id[i] = np.argwhere(check_not_nan[:,i])[-1]

E_start_beams = lower_id +1 #Fortran index
n_E_beams = upper_id - lower_id

skip_stages = np.int32([0,0,0,0,0])

n_averaging_types = n_beams
averaging_scheme = np.int32(np.arange(n_beams)+1) # Fortran index

n_beams_final= n_beams

deg = 5

n_E_out = n_E_in*5
E_grid_out = np.linspace(75, 750, num=n_E_out)


n_beams_out, e_start_beams_out, n_e_beams_out, intpol_intensity, intpol_derivative, y_func = rf.prepare_beams(E_grid_in, intensities,
	E_start_beams, n_E_beams, skip_stages, n_averaging_types,
    averaging_scheme,
	deg, E_grid_out)


res = (n_beams_out, e_start_beams_out, n_e_beams_out, intpol_intensity, intpol_derivative, y_func)

print("n_beams_out", n_beams_out)
print("e_start_beams_out", e_start_beams_out)
print("n_e_beams_out", n_e_beams_out)
print("intpol_intensity", intpol_intensity)


plt.clf()
for i in range(max(n_beams_out,2)):
    x = E_grid_out[e_start_beams_out[i]:e_start_beams_out[i]+n_e_beams_out[i]]
    plt.plot(x, intpol_intensity[:,i][e_start_beams_out[i]:e_start_beams_out[i]+n_e_beams_out[i]], label = "intpol_intensity")
    plt.plot(x, intpol_derivative[:,i][e_start_beams_out[i]:e_start_beams_out[i]+n_e_beams_out[i]], label = "intpol_derivative")
    plt.plot(x, y_func[:,i][e_start_beams_out[i]:e_start_beams_out[i]+n_e_beams_out[i]], label = "y_func")
plt.legend()
plt.show()