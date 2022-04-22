# -*- coding: utf-8 -*-

"""
Created on Apr 19 2022

@author: Alexander M. Imre

Interface between generalized search algorithm and atomic paramters
"""

from viperleed.tleedmlib.files.iorefcalc import readFdOut
from viperleed.tleedmlib.leedbase import (fortran_compile_batch, getTLEEDdir,
                                          getTensors)
import viperleed.tleedmlib.files.iorfactor as io

from viperleed.tleedmlib.wrapped.rfactor import r_factor_new as rf
from viperleed.tleedmlib.wrapped.error_codes import error_codes


# Collapse the parameter space by taking the symmetries into account
sl.atlist # atoms

# the main question now is: how big is the parameter space? i.e. how many parameters can be varies independently?
