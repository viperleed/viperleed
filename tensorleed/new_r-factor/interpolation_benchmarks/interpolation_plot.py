from scipy import interpolate
import numpy as np
import time
import matplotlib
matplotlib.use("MacOSX")
import matplotlib.pyplot as plt

nr_points = 60

a = 1.4*np.pi

x_ss2 = np.linspace(-3,3,nr_points*5)
x_ss3 = x_ss2 + 0.1
extra_random = np.random.randn(len(x_ss2))/500
x = x_ss2[::10]


def f(x, a):
    return np.cos(a*x)*np.exp(-x**2) +1

def f_1stder(x, a):
    return -1*(np.cos(a*x)*2*x_ss2+ a*np.sin(a*x))*np.exp(-x**2)

def f_2ndder(x, a):
    return np.exp(-x**2)*(
        2*x*(np.cos(a*x)*2*x+ a*np.sin(a*x))
        - ((a**2+2)*np.cos(a*x)-2*a*x*np.sin(a*x))
        )

def Y(i_prime, i):
    l = i_prime/i
    return l/(1+(5*l)**2)

shift = 0.2
ratio = 0.7
y_real = f(x_ss2, a) + ratio*f(x_ss2+shift, 2*a)
y_1st_derv = f_1stder(x_ss2, a) + ratio*f_1stder(x_ss2+shift, 2*a)
y_2nd_derv = f_2ndder(x_ss2, a) + ratio*f_2ndder(x_ss2+shift, 2*a)
y_real[0] = 100
y = y_real[::10]

simple_intp_1d = interpolate.interp1d(x, y,kind = "cubic", assume_sorted=True)
cubic_s = interpolate.CubicSpline(x, y, bc_type="not-a-knot")
bspline3 = interpolate.make_interp_spline(x,y , k=3, t = None, bc_type="natural", check_finite=False)
bspline5 = interpolate.make_interp_spline(x,y , k=5, t = None, bc_type="natural", check_finite=False)
univ3 = interpolate.UnivariateSpline(x,y, w = None, k = 3, s = None)
univ5 = interpolate.UnivariateSpline(x,y, w = None, k = 5, s = None)
dx = x[1]-x[0]
dydx = np.gradient(y, dx)
chermite = interpolate.CubicHermiteSpline(x,y, dydx)
akima1D = interpolate.Akima1DInterpolator(x,y)

sci_splines = {
    "interp1d cubic": (simple_intp_1d, "."),
    "CubicSpline": (cubic_s, "+"),
    "Univariate O3 - free Bspline":(univ3, "-."),
    "Univariate O5 - free Bspline": (univ5, "*"),
    "Bspline O3 forced": (bspline3,"--"),
    "Bspline O5 forced": (bspline5,"-.."),
    "Cubic Hermite": (chermite, "x"),
    "Akima PP O3": (akima1D, "^")
}

# Plot performance



x_supersampled = np.linspace(-2.5,2.5,nr_points*5)



plt.plot(x_ss2,y_real,"-", label = "real")

for name, spline in sci_splines.items():
    y_supersampled = spline[0](x_supersampled)
    plt.plot(x_supersampled, y_supersampled, spline[1], label = name)

plt.scatter(x,y, marker= "s", s=15, color = "black", label ="sampled points", zorder = 100)
plt.title("Interpolations")
plt.legend()
plt.figure()
plt.plot(x_ss2,y_1st_derv,"-", label = "real 1st derivative")
for name, spline in sci_splines.items():
    y_supersampled = spline[0](x_supersampled)
    try:
        deriv = spline[0].derivative()
    except:
        print("No drivative for " + name)
        continue
    y_derv = deriv(x_supersampled)
    plt.plot(x_supersampled, y_derv, spline[1], label = name)
plt.title("1st derivative")
plt.legend()


plt.figure()
plt.plot(x_ss2,y_2nd_derv,"-", label = "real 2nd derivative")
for name, spline in sci_splines.items():
    y_supersampled = spline[0](x_supersampled)
    try:
        deriv = spline[0].derivative().derivative()
    except:
        print("No 2nd drivative for " + name)
        continue
    y_derv = deriv(x_supersampled)
    plt.plot(x_supersampled, y_derv, spline[1], label = name)
plt.title("2nd derivative")
plt.legend()

plt.figure()
plt.plot(x_ss2,Y(y_1st_derv, y_real),"-", label = "real 2nd derivative")
for name, spline in sci_splines.items():
    y_supersampled = spline[0](x_supersampled)
    try:
        deriv = spline[0].derivative()
    except:
        print("No 2nd drivative for " + name)
        continue
    y_derv = deriv(x_supersampled)
    plt.plot(x_supersampled, Y(y_derv,y_supersampled), spline[1], label = name)
plt.title("Y function")
plt.legend()

plt.show()




