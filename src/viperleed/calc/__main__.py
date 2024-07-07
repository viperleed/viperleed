"""Module __main__ of ViPErLEED (viperleed) calc.

See viperleed.calc.__init__.py for more information.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

import sys

from viperleed.calc.cli import ViPErLEEDCalcCLI

ViPErLEEDCalcCLI.run_as_script()
