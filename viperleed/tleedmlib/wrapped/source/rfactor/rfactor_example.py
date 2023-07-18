import numpy as np
from rfactor import r_factor_new as rf
import time
import timeit # dependency for testing! Don't use in Viperleed build

"""
For compiling:
f2py -m rfactor rfactor.f90 -h rfactor.pyf --overwrite-signature
f2py -c rfactor.pyf rfactor.f90 
"""

# Prepare some example arrays with curve data

y1 = np.array([-1, -0.5 ,0, 0.5, 1, 1.5, 2.5, 1.9,1,0])
y2_match = np.array([-1, -0.5 ,0, 0.5, 1, 1.5, 2.5, 1.9,1,0])
y2_bad_match = np.array([-1, -1.5 ,-1, -0.5, 0.5, 1.5, 3.5, 1.9,-2,0])

E_start1 = 20.0
E_start2 = 20.0
E_step = 1.0

V0r_shift = 0

R_pendry = None
numerator = None
denominator = None
N_overlapping_points = None

# Test calls
print("Test exec:")
print(rf.test_exec(1.0))
print("\n")

R_pendry, numerator, denominator, N_overlapping_points = rf.r_factor_beam(y1, y2_match, E_start1, E_start2, E_step, V0r_shift)

print("Matching: R, num, dem, N_ovl")
print(R_pendry)
print(numerator)
print(denominator)
print(N_overlapping_points)
print("\n")

R_pendry, numerator, denominator, N_overlapping_points = rf.r_factor_beam(y1, y2_bad_match, E_start1, E_start2, E_step, V0r_shift)

print("Non-Matching: R, num, dem, N_ovl")
print(R_pendry)
print(numerator)
print(denominator)
print(N_overlapping_points)
print("\n")

########
print("\n")
print("Trapez Integration Test:")
dx = 1
array = np.arange(-67, 185, step=dx)
print("Expected result: " + str(np.trapz(array)))
result = rf.trapez_integration_const_dx(array, dx)
print("Result: " + str(result))

############################################
#  BEAMSET
############################################

y1 = np.array([y1,y1,y1])
y2 = np.array([y2_match,y2_bad_match, y2_match])
sizes_y1 = [10,10,10]
sizes_y2 = [10,10,10]
E_start1 = np.array([20,20,20])
E_start2 = np.array([20,18,20])
E_step = 1
V0r_shift = 0
beamtypes = np.array([0,0,0])

R_Pe_weighted, R_Pe_beams, N_overlapping_points = rf.rfactor_beamset(y1, sizes_y1, y2, sizes_y2, E_start1, E_start2, E_step, V0r_shift, beamtypes)

print("Beamset Test:")
print(R_Pe_weighted)
print(R_Pe_beams)
print(N_overlapping_points)