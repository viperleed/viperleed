"""Module cli of viperleed.gui.

Defines the ViPErLEEDGUICLI class, the main entry point for starting
the ViPErLEED Graphical User Interface.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2020-01-11'
__license__ = 'GPLv3+'

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
from viperleed.gui.base import catch_gui_crash
from viperleed.gui.constants import LOGO
from viperleed.gui.detect_graphics import has_graphics
from viperleed.gui.detect_graphics import has_pyqt
from viperleed.gui.helpers import resources_path

if has_pyqt():
    from viperleed.gui.selectplugin import ViPErLEEDSelectPlugin


class ViPErLEEDGUICLI(ViPErLEEDCLI, cli_name='gui'):
    """Main entry point for GUI. Both command-line and graphical."""

    def __call__(self, args=None):
        """Call either the CLI or graphical versions of the GUI."""
        args = self.parse_cli_args(args)
        if is_commandline_mode(args):
            return commandline_main()
        return gui_main()

    def add_parser_arguments(self, parser):
        """Add CLI arguments for viperleed.gui to `parser`."""
        super().add_parser_arguments(parser)
        parser.add_argument(
            '--nogui',
            help=('run the ViPErLEED graphical user interface in '
                  'command-line mode, i.e., without any windows'),
            action='store_true',
            )


def is_commandline_mode(args):
    """Return whether the GUI should run in command line mode.

    This is the case if, e.g., the correct modules are not installed,
    or if the user asked to run in command line mode via the --nogui
    command line argument.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    bool
        Whether the GUI should run in command-line mode, i.e.,
        without any windows.
    """
    return (
        not has_pyqt()
        or not has_graphics()
        or args.nogui
        )


def commandline_main():
    """Body of the command-line version of the GUI."""
    print('Running in command line version', flush=True, end='')
    if not has_pyqt():
        print(' because PyQt5 was not found', flush=True, end='')
    print('...', flush=True, end='')
    print('Not implemented yet.', flush=True)
    return 0


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
    font_path = resources_path('gui/fonts')
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
    print('Done', flush=True)
    return app.exec_()
