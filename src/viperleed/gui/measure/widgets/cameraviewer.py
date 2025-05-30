"""Module cameraviewer of viperleed.gui.measure.widgets.

Defines the CameraViewer class, a QScrollArea that allows
visualizing frames from a concrete subclass of CameraABC.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-10-29'
__license__ = 'GPLv3+'

from copy import deepcopy
from functools import partial
from math import ceil
import weakref

import numpy as np
from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg
from PyQt5 import QtWidgets as qtw

from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.camera.imageprocess import ImageProcessor
from viperleed.gui.measure.dialogs.settingsdialog import SettingsDialog
from viperleed.gui.measure.widgets.imageviewer import ImageViewer
from viperleed.gui.measure.widgets.roi import RegionOfInterest
from viperleed.gui.widgetslib import move_to_front
from viperleed.gui.widgetslib import screen_fraction


# TODO: ImageViewer.optimum_size is not updated when screen is changed
# TODO: ROI show size in image coordinates as it is resized (tooltip?)
# TODO: CameraViewer add controls for same actions as context menu,
#       plus some basic properties (probably only exposure & gain).
#       May be done with an expandable toolbar with a "..." button?
#       One could also have an "auto-contrast strength" field that
#       does not make much sense in a context menu. F5 could be used
#       for the purpose when "auto-contrast" is active.
# TODO: WEIRD: auto-contrast + saturation overlay is very slow,
#       while the two separately are not too bad (saturation mask
#       seems worse) <-- optimize a bit?


_RED_BORDER_WIDTH = 3  # pixels
_UNIQUE = qtc.Qt.UniqueConnection

# Each entry is ("short_name", default_value, "long name", always_active)
# Those with an empty long name will not be added as context
# menu entries, and will thus not be editable by the user
_DEFAULT_FLAGS = (
    ("roi_visible", True, "Allow setting ROI", False),
    ("show_auto", True, "Show on new frames", False),
    ("stop_on_close", True, "Stop camera when closed", False),
    ("auto_contrast", False, "Auto-adjust contrast", True),
    ("interactions_enabled", True, "", None),
    )

# pylint: disable=too-many-instance-attributes
# Disabled because pylint counts also class attributes, but
# pyqtSignals cannot be grouped into a container like other
# attributes and returned as @property.
class CameraViewer(qtw.QScrollArea):
    """Class for displaying camera frames with scroll bars.

    Zooming can be performed when pressing '+' and '-' or by
    holding down the Control key (Command on Mac) while scrolling.
    A region of interest can be selected and moved/resized.

    Attributes
    ----------
    image_size : QtCore.QSize
        Size in pixel (without scaling) of the image shown.
    roi : RegionOfInterest
        Region of interest object.
    roi_visible : bool
        Whether the ROI is made visible upon user interaction.
    screen : QtGui.QScreen or None
        The screen on which this widget is currently shown.
    show_auto : bool
        Whether the widget is made visible when new frames arrive.
    stop_on_close : bool
        Whether the camera is stopped when the window is closed.

    Public methods
    --------------
    scale_image(by_factor)
        Resize self and the image shown by the given factor.

    Reimplemented methods
    ---------------------
    changeEvent(event)
        Handle window maximization.
    closeEvent(event)
        Stop camera when user closes the window.
    keyPressEvent(event)
        Handle zooming and ROI visibility/scaling/movements.
    mouseDoubleClickEvent(event)
        Prevent propagation to parent if ROI is visible.
    mouseMoveEvent(event)
        Draw a ROI.
    mousePressEvent(event)
        Initiate ROI drawing.
    mouseReleaseEvent(event)
        Prevent propagation to parent if ROI is visible.
    showEvent(event)
        Resize self and image when shown.
    sizeHint()
        Optimal size is bound to image size.
    wheelEvent(event)
        Handle zooming.

    Signals
    -------
    shown : no argument
        Emitted each time this object is made visible.
    image_size_changed : no argument
        Emitted each time the number of pixels acquired by the
        camera changes. Access the new size via self.image_size.
    zoom_changed : float
        Emitted each time the scaling factor used to display
        camera frames changes. Carries the new zoom factor.
    """

    # _viewers is a dict of {camera: viewer} used to prevent
    # creating multiple viewers for the same camera. It holds
    # weak references to the viewer objects, so that we don't
    # risk creating cyclic references.
    _viewers = {}

    shown = qtc.pyqtSignal()
    image_size_changed = qtc.pyqtSignal()
    zoom_changed = qtc.pyqtSignal(float)  # new zoom factor

    def __new__(cls, camera, *args, **kwargs):
        """Return a CameraViewer instance.

        This method is extended such that only one viewer
        can exists for a given camera.

        Parameters
        ----------
        camera : CameraABC
            The camera object whose frames are displayed.
        *args : object
            Other positional arguments not used within__new__.
        **kwargs : object
            Optional arguments not used within __new__.

        Returns
        -------
        viewer : CameraViewer
            The viewer object.
        """
        if camera in cls._viewers:
            viewer = cls._viewers[camera]()
            if viewer is not None:
                # Viewer is alive
                return viewer
            # Viewer has been garbage-collected. Remove it from cache.
            cls._viewers.pop(camera)
        new_viewer = super().__new__(cls, camera, *args, **kwargs)
        new_viewer._initialized = False   # Call __init__ only once
        cls._viewers[camera] = weakref.ref(new_viewer)
        return new_viewer

    def __del__(self):
        """Remove self from the dictionary of viewers."""
        try:
            self._viewers.pop(self.camera)
        except KeyError:
            # Probably cache already cleared
            pass
        try:
            super().__del__()
        except AttributeError:
            pass

    def __init__(self, camera, *args, **kwargs):
        """Initialize widget.

        Parameters
        ----------
        camera : CameraABC
            The camera object whose frames are displayed
        *args : object
            Other unused positional arguments, passed to
            QScrollArea.__init__
        parent : QWidget, optional
            Parent widget of self. Default is None.
        stop_on_close : bool, optional
            Stop the camera when this widget is closed.
            Default is True.
        show_auto : bool, optional
            If True this widget is made visible automatically
            whenever a new frame arrives from the camera (if
            it is not visible already). Otherwise, it should
            be made explicitly visible by the caller. This
            behaviour can be changed at any time by setting
            the .show_auto attribute. Default is True.
        roi_visible : bool, optional
            If True, a region of interest can be set and
            manipulated. This behaviour can be changed at any
            time by setting the .roi_visible attribute.
            Default is True.
        interactions_enabled : bool, optional
            Whether actions in the right-click menu are enabled.
            This attribute can be set to False in case the user
            is not supposed to interact with the camera. Default
            is True.
        auto_contrast : bool, optional
            Whether the visual representation of images should
            be automatically adjusted to fit the 5-95% range of
            pixel intensities. This is only a visual effect. The
            underlying data is not modified. Default is True.
        **kwargs : object
            Other unused optional arguments, passed to
            QScrollArea.__init__

        Raises
        ------
        TypeError
            If camera is not a subclass of CameraABC
        """
        _initialized = getattr(self, "_initialized", False)
        if not _initialized:
            self.__flags = {"needs_roi_limits": False,
                            "image_saturated": False}

        # Extract acceptable keyword argument flags. The second item
        # in the flags will be the corresponding QAction, constructed
        # in __compose_context_menu
        for key, default, *_ in _DEFAULT_FLAGS:
            enabled = bool(kwargs.pop(key, default))
            if not _initialized:
                self.__flags[key] = [enabled, None]
                continue
            setattr(self, key, enabled)

        if _initialized:
            return

        super().__init__(*args, **kwargs)
        self._initialized = True
        self.__glob = {
            "image_size": qtc.QSize(),
            "camera": camera,
            "mouse_button": None,
            "img_array": None,  # Used for snapping images
            "intensity_limits": (0, 2**15 - 1),
            # img_square_array is used for auto-contrast to compute
            # the standard deviation of the current image. See
            # self.__get_image()
            "img_square_array": np.zeros(1, dtype=np.uint32),
            "auto_contrast_info": (0, 1),  # intensity offset & scaling
            "actions_with_connected_camera": set(),
            }
        self.__children = {
            "viewer": ImageViewer(),
            "context_menu": {"self": qtw.QMenu(parent=self),
                             "roi":  qtw.QMenu(parent=self)},
            "settings_dialog": None,  # Created upon user request
            }
        self.__children["roi"] = RegionOfInterest(parent=self.__img_view)

        try:
            self.roi.limits = camera.get_roi_size_limits()
        except camera.exceptions:
            # Most likely camera is not open
            self.__flags["needs_roi_limits"] = True

        self.__compose()
        self.__connect()

    @classmethod
    def clear_cache(cls):
        """Clear the cache of known viewers."""
        cls._viewers = {}

    def _flag_getter(self, flag_name=''):
        """Return a flag value by its name."""
        return self.__flags[flag_name][0]

    def _flag_setter(self, enabled, flag_name=''):
        """Enable or disable a flag given its name."""
        # Usable with functools.partial(_flag_setter, flag_name=...)
        # as a property setter default in case no peculiar action has
        # to be done. Otherwise, it can be called directly in a setter
        enabled = bool(enabled)
        self.__flags[flag_name][0] = enabled
        action = self.__flags[flag_name][1]
        if action:
            action.setChecked(enabled)

    # Set up simple checkable properties. Those that require
    # more complex actions on set/get are defined below with
    # the @property decorator syntax
    auto_contrast = property(partial(_flag_getter, flag_name='auto_contrast'),
                             partial(_flag_setter, flag_name='auto_contrast'))
    show_auto = property(partial(_flag_getter, flag_name='show_auto'),
                         partial(_flag_setter, flag_name='show_auto'))
    stop_on_close = property(partial(_flag_getter, flag_name='stop_on_close'),
                             partial(_flag_setter, flag_name='stop_on_close'))

    @property
    def camera(self):
        """Return the camera whose frames are displayed."""
        return self.__glob["camera"]

    @property
    def interactions_enabled(self):
        """Return whether the user can interact with the camera."""
        return self._flag_getter("interactions_enabled")

    @interactions_enabled.setter
    def interactions_enabled(self, enabled):
        """Set whether the user can interact with the camera."""
        self._flag_setter(enabled, "interactions_enabled")
        self.__update_context_menu_action_state()

    @property
    def image_size(self):
        """Return the true size in pixel of the image currently shown."""
        return self.__glob["image_size"]

    @image_size.setter
    def image_size(self, new_size):
        """Store the true size in pixel of the image currently shown."""
        if isinstance(new_size, tuple):
            new_size = qtc.QSize(*new_size)
        if (not new_size
                or any(s <=0 for s in (new_size.width(), new_size.height()))):
            return  # Invalid size

        if new_size != self.image_size:
            self.__glob["image_size"] = new_size
            self.__img_view.updateGeometry()
            self.updateGeometry()
            self.__img_view.optimum_size = qtc.QSize()
            self.__img_view.adjustSize()
            self.adjustSize()
            self.image_size_changed.emit()

    @property
    def roi(self):
        """Return the RegionOfInterest object of self."""
        return self.__children["roi"]

    @property
    def roi_visible(self):
        """Return whether the ROI is made visible upon user interaction."""
        return self._flag_getter("roi_visible")

    @roi_visible.setter
    def roi_visible(self, visible):
        """Set the ROI to be made visible upon user interaction."""
        self._flag_setter(visible, "roi_visible")
        if not self.roi_visible:
            self.roi.hide()

    @property
    def screen(self):
        """Return a the QScreen on which self is shown."""
        try:
            screen = self.window().windowHandle().screen()
        except AttributeError:
            screen = None
        return screen

    def changeEvent(self, event):  # pylint: disable=invalid-name
        """Extend changeEvent to react to window maximization."""
        # Keep track of quantities for later adjusting the ROI
        old_scaling = self.__img_view.image_scaling
        old_roi_size = self.roi.size()
        old_roi_pos = self.roi.pos()

        if event.type() == event.WindowStateChange:
            was_maximized = event.oldState() & qtc.Qt.WindowMaximized
            if not was_maximized and self.isMaximized():
                # If window got maximized, adjust image to fit
                self.__img_view.get_scaling_to_fit(self.size()
                                                   - self.__extra_size)
                self.__img_view.scale_to_image()

            elif was_maximized and not self.isMaximized():
                # If it is turned from maximized to normal,
                # scale it to optimum (i.e., same as after shown)
                self.__img_view.scale_to_optimum()
                self.adjustSize()

            by_factor = self.__img_view.image_scaling / old_scaling
            if abs(by_factor - 1) > 1e-3:
                # Image has been scaled. Adjust also the ROI.
                new_roi = qtc.QRect(by_factor * old_roi_pos,
                                    by_factor * old_roi_size)
                self.roi.setGeometry(new_roi)

        super().changeEvent(event)

    def closeEvent(self, event):  # pylint: disable=invalid-name
        """Extend to stop camera when window is closed by the user."""
        camera = self.camera
        if event.spontaneous() and camera.is_running and self.stop_on_close:
            # Disconnect signals that may pop up the window
            # again should a frame arrive in the meantime
            if self.show_auto:
                base.safe_disconnect(camera.image_processed, self.__show_image)
                base.safe_disconnect(camera.frame_ready, self.__show_image)
            camera.stop()
            qtw.qApp.processEvents()
        try:
            self.__children['settings_dialog'].reject()
        except AttributeError:
            # Probably never requested
            pass
        super().closeEvent(event)

    def keyPressEvent(self, event):  # pylint: disable=invalid-name
        """Extend keyPressEvent for zooming and ROI movements.

        Pressing '+' or '-' resizes the image (and self),
        pressing 'Esc' hides the ROI, pressing arrow keys
        while holding down the Control key (Command on Mac)
        affects the ROI: If Alt is also pressed, the ROI is
        resized, otherwise it is moved. If Shift is pressed,
        movements/scalings are 10 times larger.

        Parameters
        ----------
        event : QtGui.QKeyEvent
            The event corresponding to a key press

        Returns
        -------
        None.
        """
        # pylint: disable=no-member
        #    I could not find a solution to the no-member
        #    errors raised by pylint for Key.Key_*
        #    Should be fixed by white-listing PyQt5, but it
        #    is possible that Key_* attributes are added at
        #    runtime depending on the system.
        if event.key() == qtc.Qt.Key.Key_Plus:
            self.__zoom('in')
        elif event.key() == qtc.Qt.Key.Key_Minus:
            self.__zoom('out')
        elif event.key() == qtc.Qt.Key.Key_Escape:
            self.roi.hide()
        elif (self.roi.isVisible()
              and event.modifiers() & qtc.Qt.ControlModifier
              and event.key() in (qtc.Qt.Key_Left, qtc.Qt.Key_Right,
                                  qtc.Qt.Key_Up, qtc.Qt.Key_Down)):
            # Move/resize ROI if visible
            move_resize = qtc.QPoint(0, 0)
            delta = 10 if event.modifiers() & qtc.Qt.ShiftModifier else 1
            if event.key() == qtc.Qt.Key_Left:
                move_resize += qtc.QPoint(-delta, 0)
            if event.key() == qtc.Qt.Key_Right:
                move_resize += qtc.QPoint(delta, 0)
            if event.key() == qtc.Qt.Key_Up:
                move_resize += qtc.QPoint(0, -delta)
            if event.key() == qtc.Qt.Key_Down:
                move_resize += qtc.QPoint(0, delta)
            if event.modifiers() & qtc.Qt.AltModifier:
                self.roi.scale(move_resize)
            else:
                self.roi.translate(move_resize)
            return
        # pylint: enable=no-member
        super().keyPressEvent(event)

    # pylint: disable=invalid-name
    def mouseDoubleClickEvent(self, event):
        """Prevent parent propagation if ROI is visible."""
        if not self.roi.isVisible():
            super().mouseDoubleClickEvent(event)
    # pylint: enable=invalid-name

    def mouseMoveEvent(self, event):     # pylint: disable=invalid-name
        """Extend mouseMoveEvent to draw a ROI (if visible).

        If Control (Command on Mac) is held down, the ROI will
        be square.

        Parameters
        ----------
        event : QtGui.QMouseEvent
            The mouse move event that triggered this call.

        Returns
        -------
        None.
        """
        if not self.roi_visible or self.__mouse_button != qtc.Qt.LeftButton:
            super().mouseMoveEvent(event)
            return
        mouse_pos = self.roi.parent().mapFromGlobal(event.globalPos())
        self.roi.setGeometry(qtc.QRect(self.roi.origin, mouse_pos))
        if not self.roi.isVisible():
            self.roi.show()

    def mousePressEvent(self, event):    # pylint: disable=invalid-name
        """Extend mousePressEvent to begin drawing a ROI."""
        self.__mouse_button = event.button()
        if not self.roi_visible or event.button() != qtc.Qt.LeftButton:
            super().mousePressEvent(event)
            return
        origin = self.roi.parent().mapFromGlobal(event.globalPos())
        self.roi.initiate_draw(origin)

    def mouseReleaseEvent(self, event):  # pylint: disable=invalid-name
        """Extend mouseReleaseEvent to handle ROI drawing."""
        self.__mouse_button = None
        self.roi.finish_draw()
        if not self.roi.isVisible():
            # Allow propagating event to parent only if
            # not handled already by the ROI movement
            super().mouseReleaseEvent(event)

    def scale_image(self, by_factor):
        """Scale self and image by the given factor."""
        # Keep track of the (scaled) offset and size of the region
        # of interest. NB: this should be done before actually scaling
        # the underlying image, as this also changes ROI size limits.
        scaled_roi_offset = (by_factor - 1) * self.roi.pos()
        scaled_roi_size = by_factor * self.roi.size()

        # Rescale the image viewer
        self.__img_view.image_scaling = (by_factor
                                         * self.__img_view.image_scaling)
        self.__img_view.scale_to_image()

        # Rescale and reposition the region of interest, if visible
        if self.roi.isVisible():
            new_roi = qtc.QRect(self.roi.pos() + scaled_roi_offset,
                                scaled_roi_size)
            self.roi.setGeometry(new_roi)

        # Now decide whether it makes sense to scale also self
        new_self_size = self.__img_view.size() + self.__extra_size
        if self.isMaximized():
            # Do not rescale if maximized
            scr_fraction = -1
        else:
            # See if rescaling would make self go out of the screen
            bottom_right_delta = (
                qtc.QSize(self.pos().x(), self.pos().y())
                + new_self_size
                + qtc.QSize(self.frameWidth(), self.frameWidth())
                - self.screen.availableSize())
            if any(s > 0 for s in (bottom_right_delta.width(),
                                   bottom_right_delta.height())):
                scr_fraction = -1
            else:
                scr_fraction = screen_fraction(self, new_self_size)

        # Use 10% and 80% of the screen area as limits.
        if 0.1 < scr_fraction < 0.8:
            self.resize(new_self_size)

        self.__adjust_scroll_bars(by_factor)

    def showEvent(self, event):          # pylint: disable=invalid-name
        """Extend showEvent to rescale self and the image to optimum."""
        super().showEvent(event)
        # Distinguish whether the showEvent was triggered because
        # the window is shown by the application (.spontaneous()
        # == False), or if the user restored a minimized window.
        # In the latter case, we do not want to rescale the window
        # but merely restore it to what it was last.
        if not event.spontaneous():
            self.__img_view.scale_to_optimum()
            self.adjustSize()
        self.shown.emit()

    def sizeHint(self):                  # pylint: disable=invalid-name
        """Override sizeHint to bind to the underlying image."""
        return self.__img_view.sizeHint() + self.__extra_size

    def wheelEvent(self, event):         # pylint: disable=invalid-name
        """Extend wheelEvent for zooming while Control is pressed."""
        if not event.modifiers() & qtc.Qt.ControlModifier:
            super().wheelEvent(event)
            return
        # Mouse wheel turned with Ctrl down --> zoom
        direction = 'in' if event.angleDelta().y() > 0 else 'out'
        self.__zoom(direction=direction)

    @property
    def __checkable_flags(self):
        """Return the flags that are user-editable."""
        # pylint: disable=unsubscriptable-object
        # Appears to be a bug related to control-flow and isinstance
        return {k: v
                for k, v in self.__flags.items()
                if (isinstance(v, list) and v[1])}

    @property
    def __extra_size(self):
        """Return a bit larger size for self."""
        # Useful so that (i) scroll bars are only visible
        # if the image does not fit, and (ii) scroll bars
        # do not appear also in case we have a red border
        # drawn for saturated images.
        extra = qtc.QSize(2, 2)
        if self.__flags["image_saturated"]:
            # Make a little room for the extra red border
            _delta = ceil(_RED_BORDER_WIDTH / 2) * 2
            extra += qtc.QSize(_delta, _delta)
        return extra

    @property
    def __img_view(self):
        """Return the ImageViewer children."""
        return self.__children["viewer"]

    @property
    def __mouse_button(self):
        """Return the mouse button currently pressed, or None."""
        return self.__glob["mouse_button"]

    @__mouse_button.setter
    def __mouse_button(self, button):
        """Set the mouse button currently pressed."""
        self.__glob["mouse_button"] = button

    @property
    def __settings_roi(self):
        """Return the RoiEditor in the settings dialog or None."""
        _dlg = self.__children["settings_dialog"]
        if _dlg is None:
            return None
        return _dlg.handler['camera_settings']['roi'].handler_widget

    def __adjust_scroll_bars(self, by_factor):
        """Adjust the position of scrollbars when zooming."""
        if self.__img_view.underMouse():
            # Zoom in at the current mouse position
            mouse_offset = (self.mapFromGlobal(qtg.QCursor.pos())
                            - self.widget().pos()) * (by_factor - 1)
            for scroll_bar, mouse_offs in zip(
                    (self.horizontalScrollBar(), self.verticalScrollBar()),
                    (mouse_offset.x(), mouse_offset.y())
                    ):
                scroll_bar.setValue(scroll_bar.value() + mouse_offs)
            return

        # Zoom centered when the mouse cursor is somewhere else
        for scroll_bar in (self.horizontalScrollBar(),
                           self.verticalScrollBar()):
            scroll_bar.setValue(
                int(by_factor * scroll_bar.value()
                    + ((by_factor - 1) * scroll_bar.pageStep()/2))
                )

    def __apply_roi(self):
        """Set a new ROI in the camera."""
        # Give the camera new, updated settings. This is the easier way
        # to go, as it will update all the information in the camera.
        settings = deepcopy(self.camera.settings)
        settings["camera_settings"]["roi"] = str(self.__get_roi_abs_coords())
        was_running = self.camera.is_running
        self.camera.settings = settings
        self.camera.settings.update_file()
        if self.__settings_roi:
            self.__settings_roi.original_roi = self.camera.roi

        # Remove the rubber-band, and restart the camera.
        self.roi.hide()
        if was_running:
            self.camera.start()

    def __compose(self):
        """Place children widgets."""
        self.setSizePolicy(self.sizePolicy().Preferred,
                           self.sizePolicy().Preferred)
        self.setBackgroundRole(qtg.QPalette.Shadow)
        self.setWidget(self.__img_view)
        self.setAlignment(qtc.Qt.AlignCenter)
        self.adjustSize()
        self.setWindowTitle(self.camera.name)

        # Set up to receive context-menu requests (i.e.,
        # right click), and prepare the context menus
        self.setContextMenuPolicy(qtc.Qt.CustomContextMenu)
        self.__compose_context_menu_self()
        self.__compose_context_menu_roi()

    def __compose_context_menu_self(self):
        """Set up the context menu for self."""
        menu = self.__children["context_menu"]["self"]
        menu.addAction("Reset ROI")
        menu.addAction("Snap image")

        # Flags
        menu.addSeparator()
        for key, _, act_text, always_active in _DEFAULT_FLAGS:
            if not act_text:
                continue
            act = self.__flags[key][1] = menu.addAction(act_text)
            act.setCheckable(True)
            act.setData(always_active)
            enabled = getattr(self, key)
            act.setChecked(qtc.Qt.Checked if enabled else qtc.Qt.Unchecked)

        # Device settings
        menu.addSeparator()
        act = menu.addAction("Properties")
        self.__glob['actions_with_connected_camera'].add(act)
        act.setData(self.camera.connected)

        self.__update_context_menu_action_state()

    def __compose_context_menu_roi(self):
        """Set up the context menu for the ROI rubber-band."""
        menu = self.__children["context_menu"]["roi"]

        menu.addAction("Apply ROI")
        act = menu.addAction("Set ROI coordinates")
        act.setEnabled(False)                                                   # TODO: open up a dialog with a ROIEditor

    def __connect(self):
        """Connect signals."""
        self.camera.connection_changed.connect(
            self._on_camera_connection_changed
            )
        self.image_size_changed.connect(self.__update_title)
        self.__img_view.image_scaling_changed.connect(self.__on_image_scaled)
        self.camera.started.connect(self.__on_camera_started)
        self.camera.error_occurred.connect(self.__on_camera_error)
        self.customContextMenuRequested.connect(self.__show_context_menu)
        for menu in self.__children["context_menu"].values():
            menu.triggered.connect(self.__on_context_menu_triggered)
        self.roi.roi_changed.connect(self.__on_roi_rubberband_changed)

    def __get_image(self):
        """Return a QImage ready for display."""
        img_array = self.__glob['img_array']

        # Adjust contrast if necessary. Computing 5-95% percentiles
        # as one would normally do (using numpy.percentile) is too
        # slow. Let's use mean +- 2.5 * st_dev. The standard deviation
        # is computed manually, as numpy.std is quite slow. Notice that
        # we use uint32 for the array that contains the square of the
        # image, as this prevents overflowing.
        _min, _scale = 0, 1
        if self.auto_contrast:
            _dmax = 65535  # 2**16 - 1
            _mean = img_array.mean()
            _square = self.__glob['img_square_array']
            if _square.shape != img_array.shape:
                _square = np.zeros_like(img_array, dtype=np.uint32)
                self.__glob['img_square_array'] = _square
            _square[:, :] = img_array
            _square *= _square
            _std = np.sqrt(_square.mean() - _mean**2)
            _min = np.uint16(max(_mean - 2.5 * _std, 0))
            _max = np.uint16(min(_mean + 2.5 *_std, _dmax))
            _scale = np.uint16(_dmax / (_max - _min))                          # TODO: Michael suggested having a _clip_max
                                                                               #= _dmax/scale + _min
            img_array = img_array.clip(_min, _max, dtype=img_array.dtype)
            # Adjust range to 0 -- (2**16 - 1)
            img_array -= _min
            img_array *= _scale
        self.__glob["auto_contrast_info"] = (_min, _scale)

        # Notice the swapping of width and height!
        height, width = img_array.shape
        return qtg.QImage(img_array, width, height, img_array.strides[0],
                          qtg.QImage.Format_Grayscale16)

    def __get_roi_abs_coords(self):
        """Return the absolute coordinates of the current roi."""
        # This method is necessary because the ROI rubber-band
        # offset is always relative to the ROI currently applied
        # in the camera. This method returns, instead the ROI
        # in the sensor's reference. The currently applied ROI
        # is stored in self.__settings_roi.original_roi (if a
        # settings dialog was already created) or in camera.roi
        # (before a dialog was created)
        if self.__settings_roi:
            original_roi = self.__settings_roi.roi
        else:
            original_roi = self.camera.roi
        offs_x, offs_y, *_ = original_roi
        rel_x, rel_y, width, height = self.roi.image_coordinates
        return rel_x + offs_x, rel_y + offs_y, width, height

    def __get_saturation_mask(self, img_array):
        """Return a QRegion saturation mask for image."""
        _min, _max = self.__glob['intensity_limits']
        height, width = img_array.shape

        # Pick pixels that essentially saturate,
        # taking into account auto-contrast scaling
        threshold = _min + .99*(_max - _min)
        offs, scale = self.__glob["auto_contrast_info"]

        threshold = scale*(threshold - offs)                                   # TODO: Make Threshold independent of scale.
                                                                               # Just take high uint16 value.
        # TODO: here we may get it wrong if we have                            # Currently the threshold is multiplied by
        # many bad pixels! Consider applying a bad                             # the scale which may lower it to
        # pixels correction to images before showing.                          # unreasonable values.
        # Alternatively: use self.camera.bad_pixels
        # to decide depending on how many there are
        saturation_arr = img_array > threshold
        if np.sum(saturation_arr) <= 5:
            return None

        # Now construct a QRegion mask. We will use it as
        # a clipping region for filling with solid red
        saturation_arr = np.packbits(saturation_arr, axis=1, bitorder='little')
        mask = qtg.QImage(saturation_arr, width, height,
                          # pylint: disable-next=E0606,E1136  # Bug?
                          saturation_arr.strides[0],
                          qtg.QImage.Format_MonoLSB)
        return qtg.QRegion(qtg.QBitmap.fromImage(mask))

    def __make_settings_dialog(self):
        """Make and return the camera settings dialog."""
        dlg = self.__children['settings_dialog']
        if dlg:
            return dlg

        dlg = SettingsDialog(handled_obj=self.camera)
        self.__children['settings_dialog'] = dlg

        dlg.settings_changed.connect(self.__on_settings_changed)
        dlg.settings_saved.connect(self.__on_settings_saved)
        dlg.finished.connect(self.__restore_roi_visibility)
        self.__settings_roi.roi_changed.connect(self.__on_roi_settings_changed)
        return dlg

    @qtc.pyqtSlot(bool)
    def _on_camera_connection_changed(self, connected):
        """React to a connection status change of the camera."""
        # Enable/disable stuff that should only be accessed when a
        # camera is currently connected. Also, store information
        # in the action data, so it can be correctly restored
        # when interactions are enabled/disables externally.
        enabled = self.interactions_enabled and connected
        for action in self.__glob['actions_with_connected_camera']:
            action.setEnabled(enabled)
            action.setData(connected)

    @qtc.pyqtSlot(tuple)
    def __on_camera_error(self, _):
        """Close viewer when errors occur."""
        if self.isVisible():                                                    # TODO: a better way -- close only on "fatal" errors [once we have the error levels implemented]
            self.close()

    @qtc.pyqtSlot()
    def __on_camera_started(self):
        """React to a start of the camera."""
        mode = self.camera.mode
        if mode == 'triggered':
            connect_to = self.camera.image_processed
            disconnect_from = self.camera.frame_ready
        else:
            connect_to = self.camera.frame_ready
            disconnect_from = self.camera.image_processed

        base.safe_disconnect(disconnect_from, self.__show_image)
        base.safe_connect(connect_to, self.__show_image, type=_UNIQUE)
        self.__glob["intensity_limits"] = self.camera.intensity_limits
        if self.__flags["needs_roi_limits"]:
            self.__flags["needs_roi_limits"] = False
            self.roi.limits = self.camera.get_roi_size_limits()

    def __on_context_menu_triggered(self, action):
        """React to a selection in the context menu."""
        # See if action is one of the checkable flags
        checkable = {a: k for k, (_, a) in self.__checkable_flags.items()}
        key = checkable.get(action, None)
        if key:
            # Get the property object, and set its value
            setattr(self, key, action.isChecked())
            return

        # Resort to looking at the text, as we
        # are not storing the other actions
        text = action.text().lower()
        if "reset" in text:  # reset ROI
            self.__reset_roi()
            return
        if "snap" in text:
            self.__on_snap_image()
            return
        if "properties" in text:
            self.__on_settings_dialog_open()
            return
        if "apply" in text:
            self.__apply_roi()
            return
        print(action.text(), ": not implemented yet")

    def __on_image_scaled(self):
        """React to a change of zoom factor."""
        new_scaling = self.__img_view.image_scaling
        self.roi.image_scaling = new_scaling

        title, *_ = self.windowTitle().split(' - ')
        title += f' - {100*new_scaling:.1f}%'
        self.setWindowTitle(title)
        self.zoom_changed.emit(new_scaling)

    @qtc.pyqtSlot()
    def __on_roi_rubberband_changed(self):
        """React to a change in the ROI rubber-band."""
        if not self.roi.is_being_edited or not self.__settings_roi:
            # Not a user-induced edit but just an update,
            # or the settings dialog was not created yet
            return

        new_roi = self.__get_roi_abs_coords()
        if new_roi != self.__settings_roi.roi:
            self.__settings_roi.roi = new_roi

    @qtc.pyqtSlot(tuple)
    def __on_roi_settings_changed(self, new_roi):
        """React to a user change of the ROI coordinates."""
        if self.roi.is_being_edited:
            return

        if not self.roi_visible:
            self.roi_visible = True
        # The image coordinates are relative to the ROI currently applied
        offs_x, offs_y, *_ = self.__settings_roi.original_roi
        roi_x, roi_y, roi_w, roi_h = new_roi
        new_roi = (roi_x - offs_x, roi_y - offs_y, roi_w, roi_h)
        if new_roi != self.roi.image_coordinates:
            self.roi.image_coordinates = new_roi
        if not self.roi.isVisible():
            self.roi.show()

    def __on_settings_changed(self):
        """React to a user press of "Apply"/"Ok" in the settings dialog."""
        _dialog = self.__children["settings_dialog"]
        self.__restore_roi_visibility()
        was_running = self.camera.is_running
        self.camera.settings = _dialog.settings
        if was_running:
            self.camera.start()

    def __restore_roi_visibility(self):
        """Hide rubberband if it was hidden when settings_dialog opened."""
        dlg = self.__children["settings_dialog"]
        try:
            self.roi_visible = dlg.roi_was_editable
        except AttributeError:
            # Dialog was never opened. No point
            # doing anything about the ROI
            return

        if not dlg.roi_was_visible:
            self.roi.hide()

    def __on_settings_dialog_open(self):
        """React to a user requesting to view the settings dialog."""
        dlg = self.__children['settings_dialog']
        if dlg is None:
            dlg = self.__make_settings_dialog()

        if dlg.isVisible():
            move_to_front(dlg)
            return

        # Remember the visibility state of the ROI, so
        # we can set it back when we close the dialog
        dlg.roi_was_visible = self.roi.isVisible()
        dlg.roi_was_editable = self.roi_visible
        current_roi = ()
        if not self.roi.isVisible():
            # Temporarily disconnecting prevents that changes
            # to the values in the ROIEditor upon showEvent
            # trigger a settings_changed that in turn makes a
            # hidden ROI visible.
            base.safe_disconnect(self.__settings_roi.roi_changed,
                                 self.__on_roi_settings_changed)
        else:
            # If the rubberband was already visible, keep track of
            # its current coordinates, as we will use them to update
            # the ROI editor in the dialog.
            current_roi = self.__get_roi_abs_coords()
        dlg.open()

        if current_roi:
            self.__settings_roi.roi = current_roi
        base.safe_connect(self.__settings_roi.roi_changed,
                          self.__on_roi_settings_changed, type=_UNIQUE)

    def __on_settings_saved(self, saved):
        """Restore original settings in case of unsaved changes."""
        if saved:
            return
        was_running = self.camera.is_running
        old_settings = self.camera.settings
        old_settings.read_again()
        self.camera.settings = old_settings
        if was_running:
            self.camera.start()

    def __on_snap_image(self):
        """Save current frame to file."""
        image = self.__glob['img_array'].copy()
        fname, _ = qtw.QFileDialog.getSaveFileName(                             # TODO: add default filename with date/time and directory?
            parent=self,
            filter="TIFF Image (*.tiff *.tif)"
            )
        if not fname:  # No file selected
            return
        if not (fname.endswith('.tif') or fname.endswith('.tiff')):
            fname += '.tiff'
        process_info = self.camera.process_info.copy()
        process_info.filename = fname
        process_info.n_frames = 1
        processor = ImageProcessor()
        processor.prepare_to_process(process_info, image)
        processor.process_frame(image)

    def __reset_roi(self):
        """Reset camera ROI to full sensor and update settings."""
        was_running = self.camera.is_running
        if was_running:
            self.camera.stop()
        self.camera.set_roi(no_roi=True)
        self.camera.settings.set("camera_settings", "roi",
                                 str(self.camera.get_roi()))
        self.camera.settings.update_file()
        if self.__settings_roi:
            self.__settings_roi.original_roi = self.camera.roi
            self.__settings_roi.roi = self.camera.roi
        self.roi.hide()
        try:
            self.camera.bad_pixels.apply_roi(no_roi=True)
        except RuntimeError:
            # There were no bad pixels to begin with.
            # Notice that we cannot use if .bad_pixels, as
            # this checks only those in the (now old) ROI
            pass
        if was_running:
            self.camera.start()

    def __show_context_menu(self, position):
        """Show a context menu when right-clicking at position."""
        if self.roi.isVisible() and self.roi.geometry().contains(position):
            menu = self.__children["context_menu"]["roi"]
        else:
            menu = self.__children["context_menu"]["self"]
        menu.popup(self.mapToGlobal(position))

    @qtc.pyqtSlot(np.ndarray)
    def __show_image(self, img_array):
        """Show the gray-scale image in the img_array numpy.ndarray."""
        if img_array.dtype != np.uint16:
            # images coming from triggering are ready for TIFF,
            # and have an improper data type for viewing as a
            # QImage with Grayscale16 format.
            img_array = np.uint16(img_array)

        # If the camera does not support ROI at the hardware
        # level, show only the portion of the image inside
        # the ROI.
        supports_roi = ('roi' in self.camera.hardware_supported_features
                        or self.camera.get_roi())
        if not supports_roi:
            roi_x, roi_y, roi_w, roi_h = self.camera.roi
            img_array = img_array[roi_y:roi_y+roi_h, roi_x:roi_x+roi_w]
        self.__glob['img_array'] = img_array

        # Get image, including contrast adjustment if necessary
        image = self.__get_image()
        mask = self.__get_saturation_mask(img_array)
        self.__flags["image_saturated"] = mask is not None
        self.__update_frame_style()

        self.__img_view.set_image(image, mask)
        self.image_size = image.size()
        if self.show_auto and not self.isVisible():
            self.show()

    def __update_context_menu_action_state(self):
        """Update enabled state of action in the context menu."""
        camera_connected_only = self.__glob['actions_with_connected_camera']
        for action in self.__children["context_menu"]["self"].actions():
            active = action.data()
            if action in camera_connected_only:
                enable = active and self.interactions_enabled
            else:
                enable = active or self.interactions_enabled
            action.setEnabled(enable)

    def __update_frame_style(self):
        """Pick frame style depending on whether the image is saturated."""
        stylesheet = ""
        if self.__flags["image_saturated"]:
            stylesheet = (f"{self.__class__.__name__} "
                         f"{{border: {_RED_BORDER_WIDTH}px solid red}}")
        self.setStyleSheet(stylesheet)

    def __update_title(self):
        """Update the window title with camera information."""
        img_w, img_h = self.image_size.width(), self.image_size.height()
        full_w, full_h = self.camera.sensor_size
        bin_factor = self.camera.binning

        img_size = f"{img_w}x{img_h}"
        if bin_factor == 1 and (img_w, img_h) != (full_w, full_h):
            img_size += f" of {full_w}x{full_h}"
        try:
            scale = self.windowTitle().split(' - ')[1]
        except IndexError:
            scale = ""
        title = f"{self.camera.name} ({img_size})"
        if scale:
            title += f" - {scale}"
        self.setWindowTitle(title)

    def __zoom(self, direction='in'):
        """Make image larger(smaller) by 25%(20%).

        Zooming is limited to:
            direction=='in':
                At most a factor of 8 of the original image size
            direction=='out':
                At most such that the image occupies >=10% of the screen

        Parameters
        ----------
        direction : {'in', 'out'}
            Whether zooming should enlarge the view of
            the image ('in') or make it smaller ('out')

        Returns
        -------
        None.
        """
        if direction == 'in' and self.__img_view.image_scaling < 8:
            self.scale_image(1.25)
        elif direction == 'out':
            # Limit zoom-out to 10% of the screen size
            img_size = 0.8*self.__img_view.scaled_image_size
            if screen_fraction(self, img_size) < 0.1:
                return
            self.scale_image(0.8)
# pylint: enable=too-many-instance-attributes
