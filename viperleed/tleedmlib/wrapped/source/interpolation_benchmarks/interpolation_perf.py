from scipy import interpolate
import numpy as np
import time
import matplotlib
import matplotlib.pyplot as plt

# Performance test
n_repeats = 3
supersample = 5
max_exponent = 21
min_exponent = 3


def prof_cubic_s(x,y, x_supersampled):
    cubic_s = interpolate.CubicSpline(x, y, bc_type="not-a-knot")
    y_supersampled = cubic_s(x_supersampled)
    return

def prof_univ5(x,y, x_supersampled):
    univ5 = interpolate.UnivariateSpline(x,y, w = None, k = 5, s = 0)
    y_supersampled = univ5(x_supersampled)
    return

def prof_univ3(x,y, x_supersampled):
    univ3 = interpolate.UnivariateSpline(x,y, w = None, k = 3, s = 0)
    y_supersampled = univ3(x_supersampled)
    return

def prof_bspline3(x,y, x_supersampled):
    bspline3 = interpolate.make_interp_spline(x,y , k=3, t = None, bc_type="not-a-knot", check_finite=False)
    y_supersampled = bspline3(x_supersampled)
    return

def prof_bspline5(x,y, x_supersampled):
    bspline5 = interpolate.make_interp_spline(x,y , k=5, t = None, bc_type="not-a-knot", check_finite=False)
    y_supersampled = bspline5(x_supersampled)
    return

def prof_cubic_hermite(x,y, x_supersampled):
    dx = x[1] - x[0]
    dydx = np.gradient(y, dx)
    chermite = interpolate.CubicHermiteSpline(x, y, dydx)
    y_supersampled = chermite(x_supersampled)
    return

def prof_akima1D(x,y, x_supersampled):
    akima1D = interpolate.Akima1DInterpolator(x,y)
    y_supersampled = akima1D(x_supersampled)
    return

# Performance
splines = {
    "Piecewise Cubic standard":(prof_cubic_s, []),
    "Univariate O3 - free Bspline":(prof_univ3,[]),
    "Univariate O5 - free Bspline":(prof_univ5,[]),
    "Bspline O3 forced":(prof_bspline3,[]),
    "Bspline 5th O":(prof_bspline5,[]),
    "Cubic Hermite":(prof_cubic_hermite,[]),
    "Akima PP O3":(prof_akima1D,[])
}



for name, spline in splines.items():
    func = spline[0]
    exec_times = spline[1]

    for i in range(min_exponent, max_exponent + 1):
        # initialize spline with smaller object
        nr_points = 2 ** i

        x = np.linspace(-3, 3, nr_points)
        y = np.cos(np.pi * x) * np.exp(-np.abs(x)) + np.random.randn(nr_points) / 10
        x_supersampled = np.linspace(-3, 3, nr_points * supersample)
        # time with large object


        start_time = time.time()
        for i in range(n_repeats):
            func(x, y, x_supersampled)
        end_time = time.time()
        elapsed = (end_time - start_time) / n_repeats
        splines[name][1].append(elapsed)

used_points = [2 ** i for i in range(min_exponent, max_exponent + 1)]

for name, spline in splines.items():
    plt.plot(used_points, spline[1], '-', marker = "+", label = name)


plt.legend()
plt.xscale('log')
plt.xlabel("Number of base points")
plt.ylabel("Time [ms]")
#plt.show()
plt.savefig("interpolation_profiler_log.pdf")
plt.yscale('log')
plt.savefig("interpolation_profiler_loglog.pdf")