"""Module __main__ of viperleed.calc.bookkeeper.

See viperleed.calc.bookkeeper.__init__.py for more information.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-08-03'
__license__ = 'GPLv3+'

from viperleed.calc.bookkeeper.cli import BookkeeperCLI

BookkeeperCLI.run_as_script()
