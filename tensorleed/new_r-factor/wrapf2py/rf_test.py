from cmath import exp
from unicodedata import name
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


def read_beams_csv(fname, delim = ";"):
    data = np.asfortranarray(np.genfromtxt(fname, delimiter=delim, skip_header=1))
    while np.all(np.isnan(data[:, -1])):
        data = data[:, :-1]  # remove tailing colums with only NaNs
    
    E_grid_in = data[:, 0]
    intensities = data[:, 1:]
    return E_grid_in, intensities


def prepare(data, in_grid, out_grid):
    n_E_in, n_beams = data.shape

    check_not_nan = np.invert(np.isnan(data))

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

    deg = 5

    (
        e_start_beams_out,
        n_e_beams_out,
        intpol_intensity,
        intpol_derivative,
        y_func,
    ) = rf.prepare_beams(
        in_grid,
        data,
        E_start_beams,
        n_E_beams,
        skip_stages,
        n_beams_out,
        averaging_scheme,
        deg,
        out_grid,
    )
    return (e_start_beams_out,
        n_e_beams_out,
        intpol_intensity,
        intpol_derivative,
        y_func)


if __name__ == "__main__":

    grid_exp, exp  = read_beams_csv("testdata/EXPBEAMS.csv", delim=",")
    grid_theo, theo = read_beams_csv("testdata/THEOBEAMS.csv")

    min_energy = 71
    max_energy = 799.5
    n_E = 1458
    out_grid = np.linspace(min_energy, max_energy, n_E)



    prep_exp = prepare(exp, grid_exp, out_grid)
    prep_theo = prepare(theo, grid_theo, out_grid)


    # # plt.clf()
    # for i in range(n_beams_out):
    #     x = E_grid_in[E_start_beams[i] : E_start_beams[i] + n_E_beams[i]]
    #     y = intensities[:, i][E_start_beams[i] : E_start_beams[i] + n_E_beams[i]]
    #     # plt.plot(x, y, ls= "", marker = ".", markersize= 5, label = "raw")

    #     x = E_grid_out[e_start_beams_out[i] : e_start_beams_out[i] + n_e_beams_out[i]]
    #     # plt.plot(x, intpol_intensity[:,i][e_start_beams_out[i]:e_start_beams_out[i]+n_e_beams_out[i]], label = "intpol_intensity")
    #     # plt.plot(x, intpol_derivative[:,i][e_start_beams_out[i]:e_start_beams_out[i]+n_e_beams_out[i]], label = "intpol_derivative")
    #     # plt.plot(
    #     #    x,
    #     #    y_func[:, i][e_start_beams_out[i] : e_start_beams_out[i] + n_e_beams_out[i]],
    #     #    label="y_func",
    #     # )
    # # plt.legend()
    # # plt.show()


    e_step = out_grid[1] - out_grid[0]
    with open("testdata/THEOBEAMS.csv") as file:
        line = file.readline()
        theolist  = [item.replace(" ", "") for item in line.split(";")]
        theolist[-1] = theolist[-1][:-1]
        theolist.pop(0)

    with open("testdata/EXPBEAMS.csv") as file:
        line = file.readline()
        explist = [ item.split() for item in line.split(",")]
        explist = explist[1:-2]
        explist = [item[0] for item in explist]

    exp_corr = []
    theo_corr = []
    for i, exp in enumerate(explist):
        try:
            theo_corr.append(theolist.index(exp))
            exp_corr.append(i)
        except ValueError:
            continue

    a = rf.r_pendry_beamset_y(
        e_step,
        prep_exp[-1][:, exp_corr],
        prep_theo[-1][:, theo_corr],
        prep_exp[0][exp_corr],
        prep_theo[0][theo_corr],
        prep_exp[1][exp_corr],
        prep_theo[1][theo_corr],
        0.0
    )

    # plot first curves:

    fig = plt.figure()
    i_ax = fig.add_subplot(2,1,1)
    y_ax = fig.add_subplot(2,1,2)

    i_ax.plot(out_grid, prep_theo[2][:,0]/(np.max(prep_theo[2])), label = "Theoretical")
    i_ax.plot(out_grid, prep_exp[2][:,0]/(np.max(prep_exp[2])), label = "Experimental")
    i_ax.legend()

    y_ax.plot(out_grid, prep_theo[4][:,0], label = "Theoretical")
    y_ax.plot(out_grid, prep_exp[4][:,0], label = "Experimental")
    y_ax.legend()

    fig.show()