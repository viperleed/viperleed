"""Module gui.py of viperleed.

======================================
  ViPErLEED Graphical User Interface
======================================

Author: Michele Riva
Created: 2020-01-11

This is the module that invokes the execution of the
Graphical User Interface
"""

import sys
import os

from viperleed import GLOBALS

try:
    import PyQt5.QtCore as qtc
except ImportError:
    GLOBALS['USE_GUI'] = False
else:
    import PyQt5.QtGui as qtg
    import PyQt5.QtWidgets as qtw
    GLOBALS['USE_GUI'] = True

from viperleed import guilib as gl
from viperleed.cli_base import ViPErLEEDCLI


class ViPErLEEDGUICLI(ViPErLEEDCLI, cli_name='gui'):
    """Main entry point for GUI. Both command-line and graphical."""

    def __call__(self, args=None):
        """Call either the CLI or graphical versions of the GUI."""
        if is_commandline_mode():
            return commandline_main()
            _ = super().__call__(args)                                          # TODO: This line is unreachable because we don't handle the gui arguments correctly yet
        return gui_main()


def is_commandline_mode():
    """Return whether the system requires to run in command line mode.

    This is the case if, e.g., the correct modules are not installed)
    or if the user asked to run in command line mode via the --nogui
    command line argument.
    """
    needs_commandline = not GLOBALS['USE_GUI']
    needs_commandline |= gl.BACKEND is None
    needs_commandline |= '--nogui' in sys.argv
    return needs_commandline


def commandline_main():
    """Body of the command-line version of the GUI.

    In this case, there is very little graphics involved.
    """
    print('Running in command line version...', flush=True, end='')
    print('Not implemented yet.', flush=True)


def resources_path(dir_name):
    """Return the correct path to dir_name.

    This is useful when building an executable from pyinstaller.
    When built from pyinstaller, it takes the path relative to the
    temporary path in which the "exe" is extracted during execution.

    Parameters
    ----------
    dir_name : str
        Path relative to the one from which resources_path is called.

    Returns
    -------
    str
        Modified path, only if the GUI is running from a pyinstaller
        'executable', otherwise returns dir_name unchanged.
    """
    # EVENTUALLY IT IS PROBABLY BETTER TO INCLUDE THE WHOLE /fonts FOLDER from
    # '/guilib' IN THE CORRECT PLACE, AND HAVE resources_path RETURN ITS BASE
    # PATH (i.e., the top-level folder in which the exe is) IF PYINSTALLER IS
    # USED
    if hasattr(sys, '_MEIPASS'):
        return os.path.join(sys._MEIPASS, dir_name)
    return dir_name


def gui_main():
    """Body of the GUI when running with graphics capability.

    Body of the functionality that invokes the ViPErLEED
    Graphical User Interface.
    """
    gl.catch_gui_crash()

    print('Loading GUI...', flush=True, end='')
    qtg.QGuiApplication.setAttribute(qtc.Qt.AA_UseHighDpiPixmaps)
    app = qtw.QApplication(sys.argv)

    # Import some fonts from ./fonts folder
    font_path = resources_path('guilib/fonts')
    # * Text: family =  'DejaVu Sans'
    qtg.QFontDatabase.addApplicationFont(os.path.join(font_path,
                                                      'DejaVuSans.ttf'))
    # * Math: family =  'CMU Serif'
    qtg.QFontDatabase.addApplicationFont(os.path.join(font_path, 'cmunrm.otf'))

    # leed_gui = gl.LEED_GUI()
    # leed_gui.show()
    gl.show_use_betatest_version_popup()

    print('Done', flush=True)

    # sys.exit(app.exec_())
    sys.exit()


if __name__ == '__main__':
    ViPErLEEDGUICLI.run_as_script()
