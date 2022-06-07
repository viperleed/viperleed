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