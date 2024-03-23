"""=================
    ViPErLEED
=================
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-01'
__license__ = 'GPLv3+'

import sys

from viperleed.cli import main


EXIT_CODE = main()
if EXIT_CODE is None:
    EXIT_CODE = 0
sys.exit(EXIT_CODE)
