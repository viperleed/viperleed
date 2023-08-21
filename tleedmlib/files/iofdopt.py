# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 16:37:23 2021

@author: fkrau

Functions for writing output from full-dynamic optimization
"""
import csv
import numpy as np
import logging
import copy
from numpy.polynomial import Polynomial

try:
    import matplotlib
    matplotlib.rcParams.update({'figure.max_open_warning': 0})
    matplotlib.use('Agg')  # !!! check with Michele if this causes conflicts
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    plt.style.use('viperleed.tleedm')
    from matplotlib import cm
except Exception:
    _CAN_PLOT = False
else:
    _CAN_PLOT = True

from viperleed.tleedmlib.files.iorfactor import read_rfactor_columns
from viperleed.tleedmlib.files.ivplot import plot_iv

logger = logging.getLogger("tleedm.files.iofdout")
logger.setLevel(logging.INFO)


