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


PRE_RELEASE = True


def show_pre_release_popup():
    """
    Show a pop-up dialog informing that the current version is a not meant
    to be distributed without permission, as it is a pre-release
    """
    txt = (f"ViPErLEED v{gl.GLOBALS['version']} is a pre-release not meant "
           "for public distribution!<p> Contact us at <a href="
           u''"'mailto:riva@iap.tuwien.ac.at'"'>riva@iap.tuwien.ac.at</a> '
           "</p>")
    msgBox = qtw.QMessageBox(qtw.QMessageBox.Information,
                             "Pre-release notice", txt)
    msgBox.exec()


class ViPErLEEDSelectModule(qtw.QMainWindow):
    """Window to pick what to do."""
    # modules is a dictionary of the available modules.
    # Keys are tuples, with key[0] == filename of icon
    # to use, and key[1] == class that opens the main
    # window of the module
    modules = {'pattern_simulator': ('pattern_simulator.png', gl.LEED_GUI),
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
        if any(self._open_modules.values()):
            reply = qtw.QMessageBox.question(
                self, 'Message',
                'Are you sure to quit?\nThis will close all modules!',
                qtw.QMessageBox.Yes | qtw.QMessageBox.No,
                qtw.QMessageBox.No
                )
            if reply == qtw.QMessageBox.No:
                event.ignore()
                return
        event.accept()
        qtw.qApp.quit()

    def keyPressEvent(self, event):
        """Reimpement to catch Ctrl+Q and Ctrl+W for quitting."""
        if (event.key() in (qtc.Qt.Key_Q, qtc.Qt.Key_W)
                and event.modifiers() == qtc.Qt.ControlModifier):
            self.close()

    def _compose(self):
        """Place children widgets."""
        layout = qtw.QVBoxLayout()
        logo_oneline = qtg.QPixmap(resources_path(
            'guilib/icons/viperleed_logo_oneline_200x38.png'
            ))
        logo = qtw.QLabel()
        logo.setPixmap(logo_oneline)
        layout.addWidget(logo, alignment=qtc.Qt.AlignHCenter)

        for module, (icon_fname, cls) in self.modules.items():
            icon = qtg.QIcon(os.path.join(resources_path('guilib/icons'),
                                          icon_fname))
            button = self._btns[module]
            button.setIcon(icon)
            button.setIconSize(qtc.QSize(200, 122))
            button.clicked.connect(self._on_module_open_requested)
            button.setSizePolicy(qtw.QSizePolicy.Fixed, qtw.QSizePolicy.Fixed)
            # button.setCheckable(True)
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

    @gl.print_call
    def _on_module_closed(self, module):
        """React to a module being closed."""
        for name, open_module in self._open_modules.items():
            if open_module is module:
                self._open_modules[name] = None
                break

    @gl.print_call
    def _on_module_open_requested(self, __checked):
        """Open the requested module.

        If the module is already open, bring it to the front.

        Returns
        -------
        None.
        """
        for name, button in self._btns.items():
            if self.sender() is button:
                break

        module = self._open_modules[name]
        if module:
            module.show()
            module.raise_()
            module.activateWindow()
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
                               "ViPErLEEDModuleBase class") from err
