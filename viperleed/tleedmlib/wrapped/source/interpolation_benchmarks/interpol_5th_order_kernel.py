from scipy import interpolate
import numpy as np
import time
import matplotlib
matplotlib.use("MacOSX")
import matplotlib.pyplot as plt

n_points = 30

x_ss = np.linspace(-6,6, n_points*10+1)
x= np.linspace(-6,6, n_points+1)

y = np.zeros([n_points+1])
y[(n_points)//2] = 1



cubic_s = interpolate.CubicSpline(x, y, bc_type="not-a-knot")
bspline5 = interpolate.make_interp_spline(x,y , k=5, t = None, bc_type="not-a-knot", check_finite=False)
univ5 = interpolate.UnivariateSpline(x,y, w = None, k = 5, s = 0)


y_ss_cubic =  cubic_s(x_ss)
y_bsp_5 = bspline5(x_ss)
y_univ_5 =univ5(x_ss)


np.savetxt("x.csv", x, delimiter= ",", header="x input")
np.savetxt("x_supersampled.csv", x_ss, delimiter= ",", header="x supersampled")
np.savetxt("y.csv", x, delimiter= ",", header="y input")
np.savetxt("kernel_cubic.csv", y_ss_cubic, delimiter= ",", header="y Cubic")
np.savetxt("kernel_5th_O.csv", y_bsp_5, delimiter= ",", header="y 5th Order Bspline")


plt.plot(x,y, "+", label = "input")
plt.plot(x_ss, y_ss_cubic, "-", label = "Cubic")
plt.plot(x_ss, y_bsp_5, "-.", label = "Bspline 5th O.")

plt.legend()
#plt.show()
plt.savefig("kernels_3_5.pdf", dpi = 1000)

