"""Module gui.py of viperleed.

======================================
  ViPErLEED Graphical User Interface
======================================

Author: Michele Riva
Created: 2020-01-11

This is the module that invokes the execution of the
Graphical User Interface
"""

from pathlib import Path
import sys

try:
    import PyQt5.QtCore as qtc
except ImportError:
    pass
else:
    import PyQt5.QtGui as qtg
    import PyQt5.QtWidgets as qtw

from viperleed.cli_base import ViPErLEEDCLI
from viperleed.guilib.base import catch_gui_crash
from viperleed.guilib.constants import LOGO
from viperleed.guilib.detect_graphics import has_graphics
from viperleed.guilib.detect_graphics import has_pyqt
from viperleed.guilib.helpers import resources_path

if has_pyqt():
    from viperleed.guilib.selectplugin import ViPErLEEDSelectPlugin


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

    This is the case if, e.g., the correct modules are not installed,
    or if the user asked to run in command line mode via the --nogui
    command line argument.
    """
    return (
        not has_pyqt()
        or not has_graphics()
        or '--nogui' in sys.argv
        )


def commandline_main():
    """Body of the command-line version of the GUI."""
    print('Running in command line version', flush=True, end='')
    if not has_pyqt():
        print(' because PyQt5 was not found', flush=True, end='')
    print('...', flush=True, end='')
    print('Not implemented yet.', flush=True)


def gui_main():
    """Body of the GUI when running with graphics capability.

    Body of the functionality that invokes the ViPErLEED
    Graphical User Interface.
    """
    catch_gui_crash()

    print('Loading GUI...', flush=True, end='')
    qtg.QGuiApplication.setAttribute(qtc.Qt.AA_EnableHighDpiScaling)
    qtg.QGuiApplication.setAttribute(qtc.Qt.AA_UseHighDpiPixmaps)
    app = qtw.QApplication(sys.argv)
    app.setWindowIcon(qtg.QIcon(LOGO))

    # Import some fonts from ./fonts folder
    font_path = resources_path('guilib/fonts')
    # * Text: family =  'DejaVu Sans'
    qtg.QFontDatabase.addApplicationFont(
        str(Path(font_path, 'DejaVuSans.ttf').resolve())
        )
    # * Math: family =  'CMU Serif'
    qtg.QFontDatabase.addApplicationFont(
        str(Path(font_path, 'cmunrm.otf').resolve())
        )

    plugin_selector_window = ViPErLEEDSelectPlugin()
    plugin_selector_window.show()

    ########## TODO: stuff from master to deactivate GUI
    # leed_gui = gl.LEED_GUI()
    # leed_gui.show()
    # gl.show_use_betatest_version_popup()

    print('Done', flush=True)

    sys.exit(app.exec_())
    # sys.exit()   ######## TODO: Also from master


if __name__ == '__main__':
    ViPErLEEDGUICLI.run_as_script()
