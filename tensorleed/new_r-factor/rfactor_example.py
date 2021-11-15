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

y1 = []
size_y1 = 10
y2 = []
size_y2 = 10

E_start1 = 20.0
E_start2 = 20.0
E_step = 1.0

V0r_shift = 0

R_pendry = None
numerator = None
denominator = None
N_overlapping_points = None

# Test calls
print(rf.test_exec(1.0))

rf.r_factor_beam(y1, size_y1, y2, size_y2, E_start1, E_start2, E_step, V0r_shift,
                      R_pendry, numerator, denominator, N_overlapping_points)

print(R_pendry)
print(N_overlapping_points)