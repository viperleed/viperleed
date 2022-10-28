"""Module gui.py of viperleed.

======================================
  ViPErLEED Graphical User Interface
======================================

Author: Michele Riva
Created: 2020-01-11

This is the module that invokes the execution of the
Graphical User Interface
"""

import logging
from pathlib import Path
import signal
import sys

vpr_path = Path(__file__).resolve().parents[1]
if vpr_path not in sys.path:
    sys.path.append(str(vpr_path))

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
    # USED. This is similarly the case for other resources the user may be
    # allowed to edit (e.g., config files). For this to work, one can use
    # Path(sys.executable).resolve().parent as the base path of the .exe
    # file when running from a bundled --onefile pyinstaller.
    if hasattr(sys, '_MEIPASS'):
        return str(Path(sys._MEIPASS, dir_name).resolve())
    return dir_name


def gui_main():
    """Body of the GUI when running with graphics capability.

    Body of the functionality that invokes the ViPErLEED
    Graphical User Interface.
    """
    log_path = Path(__file__).resolve().parent.parent / "_logs"
    if not log_path.exists():
        log_path.mkdir()

    print('Loading GUI...', flush=True, end='')
    qtg.QGuiApplication.setAttribute(qtc.Qt.AA_EnableHighDpiScaling)
    qtg.QGuiApplication.setAttribute(qtc.Qt.AA_UseHighDpiPixmaps)
    app = qtw.QApplication(sys.argv)
    app.setWindowIcon(qtg.QIcon(gl.pluginsbase.LOGO))

    gl.widgetslib.catch_gui_crash(log_path)
    gl.widgetslib.raise_on_qt_messages()

    # Import some fonts from ./fonts folder
    font_path = resources_path("guilib/fonts")
    # * Text: family =  'DejaVu Sans'
    qtg.QFontDatabase.addApplicationFont(
        str(Path(font_path, "DejaVuSans.ttf").resolve())
        )
    # * Math: family =  'CMU Serif'
    qtg.QFontDatabase.addApplicationFont(
        str(Path(font_path, "cmunrm.otf").resolve())
        )

    plugin_selector_window = gl.ViPErLEEDSelectPlugin()
    plugin_selector_window.show()

    print('Done', flush=True)
    
    # An awful hack to allow keyboard interrupts to be accepted also
    # while the next app.exec_ runs and the python interpreter is not
    # running. See https://stackoverflow.com/questions/4938723
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    sys.exit(app.exec_())


if __name__ == '__main__':
    if is_commandline_mode():
        commandline_main()
    else:
        gui_main()
