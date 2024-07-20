"""Module constants of viperleed.calc.bookkeeper.

Collects the definitions of a few constants used
in multiple modules of the bookkeeper package.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2020-01-30'
__license__ = 'GPLv3+'

import logging

from viperleed.calc.lib.base import parent_name


LOGGER = logging.getLogger(parent_name(__name__))
HISTORY_INFO_NAME = 'history.info'
