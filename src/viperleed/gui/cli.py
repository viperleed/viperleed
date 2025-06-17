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

from importlib import import_module
from itertools import chain
from pathlib import Path
import sys

from viperleed.cli_base import ViPErLEEDCLI
from viperleed.gui.base import catch_gui_crash
from viperleed.gui.constants import LOGO
from viperleed.gui.detect_graphics import Qt5DependencyFinder
from viperleed.gui.detect_graphics import find_missing_qt_dependencies
from viperleed.gui.detect_graphics import has_graphics
from viperleed.gui.detect_graphics import has_pyqt
from viperleed.gui.detect_graphics import suppress_file_permission_warnings
from viperleed.gui.helpers import resources_path


class ViPErLEEDGUICLI(ViPErLEEDCLI, cli_name='gui'):
    """Main entry point for GUI. Both command-line and graphical."""

    def __call__(self, args=None):
        """Call either the CLI or graphical versions of the GUI."""
        args = self.parse_cli_args(args)
        if args.nogui:
            return commandline_main()
        self.check_can_run_gui()
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

    def check_can_run_gui(self):
        """Raise SystemExit if the GUI cannot execute.

        Raises
        ------
        SystemExit
            If the GUI cannot execute because
            - PyQt5 is not installed,
            - PyQt5 misses any system dependencies,
            - the current machine has no graphics capability.
        """
        err_msg = (
            'Cannot execute the ViPErLEED graphical user interface because {}'
            )
        if not has_pyqt():
            self.parser.error(err_msg.format('PyQt5 is not installed.'))
        missing = tuple(chain.from_iterable(
            find_missing_qt_dependencies().values()
            ))
        if missing:
            sep = '\n    '
            fmt_missing = sep.join(missing)
            fmt_missing = sep + fmt_missing
            how_to_install = Qt5DependencyFinder.find_install_for_libs(missing)
            install_msg = 'Try again after installing them'
            install_msg += ('.' if not how_to_install
                            else f' via{sep}{how_to_install}')
            self.parser.error(err_msg.format(
                'the following PyQt5 dependencies are missing on your system:'
                f'{fmt_missing}\n{install_msg}'
                ))
        if not has_graphics():
            self.parser.error(err_msg.format(
                'the system appears to have no graphics capability (i.e., no '
                'monitor was detected). Try once more if this is the first '
                'time you execute the GUI, or have just updated viperleed.'
                ))


def commandline_main():
    """Body of the command-line version of the GUI."""
    print('Running in command line version...', flush=True, end='')
    print('Not implemented yet.', flush=True)
    return 0


def gui_main():
    """Body of the GUI when running with graphics capability.

    Body of the functionality that invokes the ViPErLEED
    Graphical User Interface.
    """
    (qtc,
     qtg,
     qtw,
     plugin_selector) = import_graphics_modules()

    catch_gui_crash()

    print('Loading GUI...', flush=True, end='')
    qtg.QGuiApplication.setAttribute(qtc.Qt.AA_EnableHighDpiScaling)
    qtg.QGuiApplication.setAttribute(qtc.Qt.AA_UseHighDpiPixmaps)
    suppress_file_permission_warnings()  # Next line emits the warnings
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

    plugin_selector_window = plugin_selector()
    plugin_selector_window.show()
    print('Done', flush=True)
    return app.exec_()


def import_graphics_modules():
    """Dynamically import modules for graphical capability."""
    if not has_pyqt():
        raise ImportError('PyQt5 is not available.')
    gui_select = import_module('viperleed.gui.selectplugin')
    return (
        import_module('PyQt5.QtCore'),  # Not graphics per-se
        import_module('PyQt5.QtGui'),
        import_module('PyQt5.QtWidgets'),
        gui_select.ViPErLEEDSelectPlugin,
        )
