"""Module __main__ of viperleed.gui.

Entry point for starting the ViPErLEED Graphical User Interface.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-05-25'
__license__ = 'GPLv3+'

from viperleed.gui.cli import ViPErLEEDGUICLI

ViPErLEEDGUICLI.run_as_script()
