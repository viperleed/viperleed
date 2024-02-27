"""Module selectmodule of viperleed.guilib.

======================================
  ViPErLEED Graphical User Interface
======================================

Defines the ViPErLEEDSelectModule class, i.e., the primary window of
the ViPErLEED Graphical User Interface that allows to pick what the
user wants to do.

Author: Michele Riva
Created: 2021-06-28
"""

import os

from PyQt5 import (QtCore as qtc,
                   QtWidgets as qtw,
                   QtGui as qtg)

from viperleed import guilib as gl
from viperleed.gui import resources_path
from viperleed.guilib.pluginsbase import logo_one_line, ViPErLEEDPluginBase
from viperleed.guilib.measure.uimeasurement import Measure
from viperleed.guilib.widgetslib import move_to_front

from viperleed.guilib import decorators as dev_

PRE_RELEASE = True


def show_pre_release_popup():
    """Show a pop-up informing about pre-release status.

    The current ViPErLEED version is a not meant to be distributed
    without permission, as it is a pre-release.

    Returns
    -------
    None.
    """
    txt = (f"ViPErLEED v{gl.GLOBALS['version']} is a pre-release not meant "
           "for public distribution!<p> Contact us at <a href="
           u''"'mailto:riva@iap.tuwien.ac.at'"'>riva@iap.tuwien.ac.at</a> '
           "</p>")
    msg_box = qtw.QMessageBox(qtw.QMessageBox.Information,
                             "Pre-release notice", txt)
    msg_box.exec()


class ViPErLEEDSelectPlugin(ViPErLEEDPluginBase):
    """Window to pick what to do."""

    # modules is a dictionary of the available modules.
    # Keys are tuples, with key[0] == filename of icon
    # to use, and key[1] == class that opens the main
    # window of the module
    modules = {'pattern_simulator': ('pattern_simulator.png',
                                     gl.LEEDPatternSimulator),
               # 'measure': ('measure.png', Measure),}
               'measure': ('measure.png', None),}

    def __init__(self, parent=None):
        """Initialize window."""
        if PRE_RELEASE:
            show_pre_release_popup()

        super().__init__(parent)

        self._btns = {k: qtw.QPushButton('') for k in self.modules}
        self._open_modules = dict.fromkeys(self.modules)
        self._compose()

    def closeEvent(self, event):
        """Reimplement closeEvent to also close open modules."""
        open_modules = [module
                        for module in self._open_modules.values()
                        if module]
        if any(open_modules):
            reply = qtw.QMessageBox.question(
                self, 'Message',
                'Are you sure to quit?\nThis will close all modules!',
                qtw.QMessageBox.Yes | qtw.QMessageBox.No,
                qtw.QMessageBox.No
                )
            if reply == qtw.QMessageBox.No:
                event.ignore()
                return
            for module in open_modules:
                module.close()
        super().closeEvent(event)

    def keyPressEvent(self, event):  # pylint: disable=invalid-name
        """Reimpement to catch Ctrl+Q and Ctrl+W for quitting."""
        if (event.key() in (qtc.Qt.Key_Q, qtc.Qt.Key_W)
                and event.modifiers() == qtc.Qt.ControlModifier):
            self.close()
        super().keyPressEvent(event)

    def _compose(self):
        """Place children widgets."""
        layout = qtw.QVBoxLayout()

        layout.addWidget(logo_one_line(), alignment=qtc.Qt.AlignHCenter)

        for module, (icon_fname, cls) in self.modules.items():
            icon = qtg.QIcon(os.path.join(resources_path('guilib/icons'),
                                          icon_fname))
            button = self._btns[module]
            button.setIcon(icon)
            button.setIconSize(qtc.QSize(200, 122))
            button.clicked.connect(self._on_module_open_requested)
            button.setSizePolicy(qtw.QSizePolicy.Fixed, qtw.QSizePolicy.Fixed)
            layout.addWidget(button)

            if cls is None:
                button.setEnabled(False)

        central_widget = qtw.QWidget()
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)

        self.setWindowTitle('ViPErLEED')
        self.adjustSize()
        self.setFixedSize(self.size())

        flags = self.windowFlags()
        flags &= ~qtc.Qt.WindowMinimizeButtonHint
        self.setWindowFlags(flags)

        self.move(5, 5)

    @dev_.print_call
    def _on_module_closed(self, module):
        """React to a module being closed."""
        for name, open_module in self._open_modules.items():
            if open_module is module:
                self._open_modules[name] = None
                break

    @dev_.print_call
    def _on_module_open_requested(self, __checked):
        """Open the requested module.

        If the module is already open, bring it to the front.

        Returns
        -------
        None.

        Raises
        ------
        RuntimeError
            If the module being opened is not a valid ViPErLEED module
        """
        # Disable pylint warning, as _on_module_open_requested can
        # only run if there are buttons, one of which will always be
        # the sender. Thus name always has a well-defined value.

        # pylint: disable=undefined-loop-variable

        for name, button in self._btns.items():
            if self.sender() is button:
                break

        module = self._open_modules[name]

        # Module is already open
        if module:
            move_to_front(module)
            return

        _, cls = self.modules[name]

        module = cls()
        self._open_modules[name] = module
        module.show()
        try:
            module.module_closed.connect(self._on_module_closed)
        except AttributeError as err:
            raise RuntimeError(f"Module {module.__class__.__name__} is "
                               "not a valid ViPErLEED module. All modules "
                               "should be a concrete implementation of the "
                               "ViPErLEEDPluginBase class") from err
