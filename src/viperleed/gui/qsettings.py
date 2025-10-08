"""Module qsettings of viperleed.gui.

Defines functions used to work with QSettings.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-09-26'
__license__ = 'GPLv3+'

from PyQt5 import QtCore as qtc


def get_qsettings(file_name):
    """Return a QSettings instance of the given settings.

    Parameters
    ----------
    file_name : str
        The target settings for which the QSettings instance is created.

    Returns
    -------
    qsettings : QSettings
        A QSettings instance in IniFormat for the given settings.
    """
    return qtc.QSettings(qtc.QSettings.IniFormat, qtc.QSettings.UserScope,
                         'viperleed', file_name.lower())
