"""Module constants of viperleed.gui.measure."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-10-03'
__license__ = 'GPLv3+'

from pathlib import Path

DEFAULTS = Path(__file__).parent / '_defaults'
SRC_ALIASES_PATH = DEFAULTS / '_aliases.ini'
SRC_SYS_SETTINGS_PATH = DEFAULTS / '_system_settings.ini'
