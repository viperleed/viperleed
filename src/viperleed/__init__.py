"""=================
    ViPErLEED
=================
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-08-03'
__license__ = 'GPLv3+'
__version__ = '0.14.0'

import numpy as np

# Name of environment variable specifying the path to the tensor-LEED
# source code. Typically the viperleed-tensorleed repository.
VIPERLEED_TENSORLEED_ENV = 'VIPERLEED_TENSORLEED'

NUMPY2_OR_LATER = np.__version__ >= '2'
