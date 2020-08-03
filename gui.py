"""
======================================
  ViPErLEED Graphical User Interface
======================================

Author: Michele Riva
Created: 2020-01-11

This is the module that invokes the execution of the Graphical User Interface

"""

import sys
import os

cd = os.path.realpath(os.path.dirname(__file__))
vpr_path = os.path.realpath(os.path.join(cd, '..'))
for import_path in (cd, vpr_path):
    if import_path not in sys.path:
        sys.path.append(import_path)

from vprglobals import GLOBALS

try:
    import PyQt5.QtCore as qtc
except ImportError:
    GLOBALS['USE_GUI'] = False
else:
    import PyQt5.QtGui as qtg
    import PyQt5.QtWidgets as qtw
    GLOBALS['USE_GUI'] = True

import guilib as gl


def is_commandline_mode():
    """
    Returns whether the system requires to run in command line mode (e.g.,
    the correct modules are not installed) or if the user asked to run in
    command line mode via the --nogui command line argument.
    """
    needs_commandline = not GLOBALS['USE_GUI']
    needs_commandline |= gl.BACKEND is None
    needs_commandline |= '--nogui' in sys.argv
    return needs_commandline


def commandline_main():
    """
    Body of the command-line version of the ViPErLEED Graphical User Interface.
    In this case, there is very little graphics involved.
    """
    print('Running in command line version...', flush=True, end='')
    print('Not implemented yet.', flush=True)


def gui_main():
    """
    Body of the functionality that invokes the ViPErLEED Graphical User
    Interface
    """
    gl.catch_gui_crash()

    print('Loading GUI...', flush=True, end='')
    qtg.QGuiApplication.setAttribute(qtc.Qt.AA_UseHighDpiPixmaps)
    app = qtw.QApplication(sys.argv)

    # Import some fonts from ./fonts folder
    # * Text: family =  'DejaVu Sans'
    qtg.QFontDatabase.addApplicationFont("guilib/fonts/DejaVuSans.ttf")
    # * Math: family =  'CMU Serif'
    qtg.QFontDatabase.addApplicationFont("guilib/fonts/cmunrm.otf")

    leedGUI = gl.LEED_GUI()
    leedGUI.show()
    
    print('Done', flush=True)
    
    sys.exit(app.exec_())


if __name__ == '__main__':    
    if is_commandline_mode():
        commandline_main()
    else:
        gui_main()