"""Module badpxfinderdialog of viperleed.guilib.measure.dialogs.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-11-22
Author: Michele Riva

Defines the BadPixelsFinderDialog class that handles user interaction
while finding bad pixels for a camera.
"""

from PyQt5 import (QtWidgets as qtw,
                   QtCore as qtc)

from viperleed.guilib.measure.hardwarebase import get_devices


class BadPixelsFinderDialog(qtw.QDialog):
    """Dialog to handle user interaction when finding bad pixels."""

    def __init__(self, parent=None):
        """Initialize dialog."""
        super().__init__(parent=parent)

        self.__ctrls = {'camera': qtw.QComboBox(),
                        'update_list': qtw.QPushButton("Update devices"),
                        'total_progress': qtw.QProgressBar(),
                        'section_progress': qtw.QProgressBar(),}
        self.__buttons = {self.__ctrls['update_list'],}

        self.__available_cameras = {}

        self.setWindowTitle("Find bad pixels")

        self.__compose()
        self.__connect()
        
        self.update_available_camera_list()

    def update_available_camera_list(self, *_):
        """Update the list of available cameras."""
        camera_combo = self.__ctrls['camera']
        old_selection = camera_combo.currentText()
        self.__available_cameras = get_devices('camera')
        camera_combo.clear()
        camera_combo.addItems(self.__available_cameras.keys())
        if old_selection:
            camera_combo.setCurrentText(old_selection)

    def __compose(self):
        """Place children widgets."""
        layout = qtw.QVBoxLayout()
        
        camera_layout = qtw.QHBoxLayout()
        camera_layout.addWidget(qtw.QLabel("Camera:"))
        camera_layout.addWidget(self.__ctrls['camera'], stretch=1)
        camera_layout.addWidget(self.__ctrls['update_list'])
        
        for btn in self.__buttons:
            btn.setAutoDefault(False)
        
        layout.addLayout(camera_layout)
        self.setLayout(layout)

    def __connect(self):
        """Connect children signals."""
        self.__ctrls['update_list'].clicked.connect(
            self.update_available_camera_list
            )
        self.__ctrls['camera'].currentTextChanged.connect(
            self.__on_camera_selected
            )
    
    def __on_camera_selected(self, camera_name):
        """React to selection of a new camera."""
        print(f"{camera_name=}")
