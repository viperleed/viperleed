from typing_extensions import assert_type
import numpy as np
import os
import sys

vpr_path = "../"

if os.path.abspath(vpr_path) not in sys.path:
    sys.path.append(os.path.abspath(vpr_path))
        
import viperleed
from viperleed.tleedmlib.wrapped.rfactor import r_factor_new as rf
from viperleed.tleedmlib.wrapped.error_codes import error_codes

class Test_rfac:
    def test_r_compile(self):
        ierr = rf.test_compilation()
        assert ierr == 0
        
        cos_func = np.cos(np.linspace(1,15,100))
        sin_func = np.cos(np.linspace(1,15,100) + np.pi)

        step = np.linspace(1,15,100)[1] - np.linspace(1,10,100)[0]
        rfac, *_ =  rf.r_pendry_beam_y(step, cos_func, cos_func, 1, 1, 100, 100, 0)
        assert abs(rfac) < 1e-10
        
    def test_r_prop2(self):
        cos_func = np.cos(np.linspace(1,15,100))
        sin_func = np.cos(np.linspace(1,15,100) + np.pi)

        step = np.linspace(1,15,100)[1] - np.linspace(1,10,100)[0]
        rfac, *_ = rf.r_pendry_beam_y(step, cos_func, sin_func, 1, 1, 100, 100, 0)
        assert abs(rfac) -1 < 1e-6
        
