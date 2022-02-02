from cmath import exp
from re import S
from unicodedata import name
import numpy as np
from rfactor import r_factor_new as rf
from interpolation import interpolation as intpol
import time
import timeit  # dependency for testing! Don't use in Viperleed build
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import random
import os



# Read Hematite test data from file EXPBEAMS.csv


def read_beams_csv(fname, delim = ";"):
    data = np.asfortranarray(np.genfromtxt(fname, delimiter=delim, skip_header=1))
    while np.all(np.isnan(data[:, -1])):
        data = data[:, :-1]  # remove tailing colums with only NaNs
    
    E_grid_in = data[:, 0]
    intensities = data[:, 1:]
    return E_grid_in, intensities

# Stolen from ivplot - reading of .column data
def read_rfactor_columns(cols_dir=''):
    """
    Reads data from the theo.column and exp.column files in a given directory.

    Parameters
    ----------
    cols_dir : str, optional
        The directory to read from. The default is ''.

    Returns
    -------
    list [theo_beams, exp_beams]
        Both the theoretical and the experimental beams are formatted as lists
        of 2D numpy arrays, each of which has the form [[en1, intens1], ...]
    """

    fnames = ['theo.column', 'exp.column']
    xxyy = []
    for fname in fnames:
        try:
            f = open(os.path.join(cols_dir, fname), 'r')
        except FileNotFoundError:
            logger.error("read_rfactor_columns: File {} not found. Aborting."
                         .format(fname))
            return []
        except PermissionError:
            logger.error("read_rfactor_columns: Cannot open file {}. Aborting."
                         .format(fname))
            return []

        cols = np.array([[float(col) for col in line.split()] for line in f])
        xy = np.split(cols, np.shape(cols)[1]/2, axis=1)
        # xy is now a list of 2D arrays.
        # Each array has the form [[en1, intens1], ...]
        #
        # for each beam, get rid of the points that have (en, intens) = (0, 0)
        # so that they don't screw up the plots later
        xy = [coords[~np.all(coords < 1e-3, axis=1)] for coords in xy]
        if xy:
            xxyy.append(xy)
        else:
            logger.warning("File " + fname + " contains no usable data.")
    # xxyy now contains first the theoretical, then the experimental beams
    return xxyy

def plot_csv(prep_exp, prep_theo, exp_corr, theo_corr, beam):
    fig = plt.figure()
    i_ax = fig.add_subplot(2,1,1)
    y_ax = fig.add_subplot(2,1,2)

    i_ax.plot(out_grid, prep_theo[2][:,beam]/(np.max(prep_theo[2])), label = "Theoretical")
    i_ax.plot(out_grid, prep_exp[2][:,beam]/(np.max(prep_exp[2])), label = "Experimental")
    i_ax.legend()

    y_ax.plot(out_grid, prep_theo[4][:,beam], label = "Theoretical")
    y_ax.plot(out_grid, prep_exp[4][:,beam], label = "Experimental")
    y_ax.legend()

    fig.show()

def plot_xxyy(xxyy, beam):
    theo_col = xxyy[0][beam]
    exp_col = xxyy[1][beam]
    plt.plot(theo_col[:,0], theo_col[:,1], label = "Theoretical")
    plt.plot(exp_col[:,0], exp_col[:,1], label = "Theoretical")
    plt.legend()
    plt.show()

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
        ierr,
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
        y_func, ierr)

def test_y_intpol(super_s_y, beam):

    

    grid_exp, exp  = read_beams_csv("testdata/EXPBEAMS.csv", delim=",")

    min_energy = 71
    max_energy = 799.5
    n_orig = 1458
    n_E = n_orig*100
    out_grid = np.linspace(min_energy, max_energy, n_E)


    prep_exp = prepare(exp, grid_exp, out_grid)

    x = out_grid[prep_exp[0][beam]:prep_exp[0][beam]+prep_exp[1][beam]]
    y_func_test = prep_exp[-1][prep_exp[0][beam]:prep_exp[0][beam]+prep_exp[1][beam], beam]
    plt.plot(x, y_func_test, label="x100", color="blue")

    ######

    min_energy = 71
    max_energy = 799.5
    n_orig = 1458
    n_E = n_orig
    out_grid = np.linspace(min_energy, max_energy, n_E)

    prep_exp = prepare(exp, grid_exp, out_grid)

    x = out_grid[prep_exp[0][beam]:prep_exp[0][beam]+prep_exp[1][beam]]
    y_func_test = prep_exp[-1][prep_exp[0][beam]:prep_exp[0][beam]+prep_exp[1][beam], beam]
    plt.plot(x, y_func_test, label=f"Y func original data grid", color="orange", ls = "", marker = "x", ms =5)

    ######

    for i in super_s_y:
        n_E = n_orig * i
        out_grid = np.linspace(min_energy, max_energy, n_E)


        prep_exp = prepare(exp, grid_exp, out_grid)

        x = out_grid[prep_exp[0][beam]:prep_exp[0][beam]+prep_exp[1][beam]]
        y_func_test = prep_exp[-1][prep_exp[0][beam]:prep_exp[0][beam]+prep_exp[1][beam], beam]

        for deg in (3, 5,):

            n_knots, nt, LHS_rows = intpol.get_array_sizes(y_func_test.shape[0], deg)

            knots, coeffs, ierr = intpol.single_calc_spline(x, y_func_test, deg, n_knots, nt, LHS_rows)

            n_new = 100*(y_func_test.shape[0]-1)

            
            # With this, the abstract representation of the Y function is achived!

            x_new = np.linspace(min_energy, max_energy, n_new)

            y_func_out = intpol.single_interpolate_coeffs_to_grid(knots, coeffs, deg, x_new, 0)


            plt.plot(x_new, y_func_out, label=f"Re-interpolated deg={deg}, Y_supersample_rate={i}", ls  = "--")

    plt.xlim(out_grid[prep_exp[0][beam]],out_grid[prep_exp[0][beam]+prep_exp[1][beam]])
    plt.ylim([-1,1])
    plt.legend()
    plt.show()


def parabola(x,a, b, c):
    return a*(x + b)**2 + c

def para2(x, a, b, c):
    return a*x**2 + b*x +c 

def poly4(x, a, b, c, d):
    return a*(x + b)**4 + c*(x + b)**2 + d

def test_parafit():
    x = np.linspace(-15, 15, 100)
    y = (x-1)**2 + 10*np.random.rand(100)
    weights = np.ones(len(x))
    plt.plot(x, y, label = "data")
    para_opt, _ = curve_fit(parabola, x, y)
    plt.plot(x, parabola(x, *para_opt), label = "scipy fit")
    fit,  ierr = rf.parabola_lsq_fit(x, y, weights)
    plt.plot(x, para2(x, fit[0], fit[1], fit[2]), label = "manual fit")
    plt.legend()
    plt.show()
    return fit

grid_exp, exp  = read_beams_csv("testdata/EXPBEAMS.csv", delim=",")
grid_theo, theo = read_beams_csv("testdata/THEOBEAMS.csv")

min_energy = 71
max_energy = 799.5
n_E = 1458 
out_grid = np.linspace(min_energy, max_energy, n_E)



prep_exp = prepare(exp, grid_exp, out_grid)
prep_theo = prepare(theo, grid_theo, out_grid)


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

xxyy = read_rfactor_columns('testdata')
theo_col = xxyy[0]
exp_col = xxyy[1]
exp_array = np.zeros([max([len(col) for col in exp_col]), len(exp_col)])
theo_array = np.zeros([max([len(col) for col in theo_col]), len(exp_col)])


if __name__ == "__main__":


    V0r_range = np.array([-20, 20], dtype="int32")
    start_guess = np.array([-3, 0, 3])
    fast_search = False

    rf.r_pendry_beamset_v0r_opt_on_grid(
        V0r_range,
        start_guess,
        fast_search,
        e_step,
        prep_exp[-2][:, exp_corr],
        prep_theo[-2][:, theo_corr],
        prep_exp[0][exp_corr],
        prep_theo[0][theo_corr],
        prep_exp[1][exp_corr],
        prep_theo[1][theo_corr])





