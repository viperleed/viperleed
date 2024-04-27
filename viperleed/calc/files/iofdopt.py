"""Functions for writing output from full-dynamic optimization."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2021-10-25'
__license__ = 'GPLv3+'

import copy
import logging

import numpy as np
from numpy.polynomial import Polynomial

try:
    import matplotlib
    matplotlib.rcParams.update({'figure.max_open_warning': 0})
    matplotlib.use('Agg')  # !!! check with Michele if this causes conflicts
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    plt.style.use('viperleed.calc')
    from matplotlib import cm
except Exception:
    _CAN_PLOT = False
else:
    _CAN_PLOT = True

from viperleed.calc.files.iorfactor import read_rfactor_columns
from viperleed.calc.files.ivplot import plot_iv


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


