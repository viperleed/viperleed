"""Module calibrationtask of viperleed.gui.measure.classes.

Defines the abstract CalibrationTask class that is to be used
as the base class of device-calibration tasks.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2022-10-12'
__license__ = 'GPLv3+'

from abc import abstractmethod
from copy import deepcopy
from enum import IntEnum

from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw

from viperleed.gui.measure import hardwarebase as base
from viperleed.gui.measure.classes.abc import QMetaABC
from viperleed.gui.measure.classes.abc import QObjectWithError


_INVOKE = qtc.QMetaObject.invokeMethod
_UNIQUE = qtc.Qt.UniqueConnection
_QUEUED = qtc.Qt.QueuedConnection


class CalibrationTaskError(base.ViPErLEEDErrorEnum):
    """Class for device-calibration-task errors."""

    DEVICE_BUSY = (500, "Cannot {} while {} is busy.")


class _InfoBox(qtw.QMessageBox):
    """Modal, information-level subclass of QMessageBox."""

    def __init__(self, *args, **kwargs):
        """Initialize message box."""
        _title = kwargs.pop('title', '')
        self.slots_ = kwargs.pop('slots', {})

        super().__init__(*args, **kwargs)
        self.setModal(True)
        self.setIcon(self.Information)
        btn = self.addButton(self.Ok)
        btn.setText("I'm ready")
        self.addButton(self.Abort)

        if _title:
            self.setWindowTitle(_title)

        for signal_name, slot in self.slots_.items():
            try:
                signal = getattr(self, signal_name)
            except AttributeError:
                continue
            base.safe_connect(signal, slot, type=_UNIQUE)

    @classmethod
    def from_info_box(cls, other):
        """Return an _InfoBox identical to other.

        The new _InfoBox instance may have affinity to a thread
        different than other.

        Parameters
        ----------
        other : _InfoBox
            The instance to copy from.

        Returns
        -------
        new_info_box : _InfoBox
            The newly created instance.
        """
        new_info_box = cls(parent=other.parent(),
                           title=other.windowTitle(),
                           slots=other.slots_)
        new_info_box.setText(other.text())
        new_info_box.setInformativeText(other.informativeText())
        return new_info_box

    @qtc.pyqtSlot(str)
    def set_info_text(self, text):
        """Set informative text in the correct thread."""
        self.setInformativeText(text)


class CalibrationTaskOperation(IntEnum):
    """Abstract class for calibration task operations.

    Subclasses of CalibrationTaskOperation are typically used
    to report progress.

    Elements should be 1-based integers, in the order in which
    operations are to be performed.

    An "operation" could be, e.g., acquisition of a set of camera
    images with given settings. A step in this operation is the
    acquisition of single images. Other examples are complex
    calculations with multiple intermediate steps.
    """

    @property
    def long_name(self):
        """Return a descriptive name for this operation.

        Must be overridden. Subclasses should return
        self._format_long_name(long_name) if the name
        has runtime-formatted portions.

        Returns
        -------
        long_name : str
            Descriptive name for this operation. Typically
            used when displaying a progress bar.

        Raises
        ------
        NotImplementedError
            Always, as this property needs to be
            overridden in subclasses.
        """
        raise NotImplementedError

    @property
    def n_steps(self):
        """Return the number of steps in this operation. Must be overridden."""
        raise NotImplementedError

    def __iter__(self):
        """Return an iterable version of self."""
        return iter((self.long_name, self.value, self.n_operations()))

    @classmethod
    def first(cls):
        """Return the first CalibrationTaskOperation."""
        return min(cls)

    @classmethod
    def n_operations(cls):
        """Return the total number of operations in this Enum."""
        return len(cls)

    def _format_long_name(self, name):
        """Return a formatted version of name."""
        _args, _kwargs = getattr(self, '_fmt_args', ((), {}))
        return name.format(*_args, **_kwargs)

    def is_last(self):
        """Return whether this is the last CalibrationTaskOperation."""
        return self.value >= self.n_operations()

    def next_(self):
        """Return the next section (or the first if this is last)."""
        next_idx = self.value + 1
        _tot = self.n_operations()
        if next_idx > _tot:
            next_idx %= _tot
        return self.__class__(next_idx)

    def set_format(self, *fmt_args, **fmt_kwargs):
        """Set positional and keyword arguments for self.long_name."""
        # pylint: disable=attribute-defined-outside-init
        # Much easier this way than extending __new__/__init__
        # for an Enum. Attribute is used only internally.
        self._fmt_args = (fmt_args, fmt_kwargs)


class CalibrationTask(QObjectWithError, metaclass=QMetaABC):
    """Abstract class for device-calibration tasks.

    Notice that CalibrationTask instances are parent-less. As such
    they can be executed in separate threads. When accessing the
    device or any GUI element, make sure to use signals and slots,
    or QtCore.QMetaObject.invokeMethod to ensure execution in the
    correct thread.
    """

    aborted = qtc.pyqtSignal()
    done = qtc.pyqtSignal()

    # Typical arguments of progress_occurred:
    # - the first one is always self.progress_name
    # - the following three can be "*operation", with
    #   operation a CalibrationTaskOperation element
    # - "current progress of the current operation" depends on
    #   when the signal is emitted. Can be an index, or 100%-based
    # - "total number of steps" can be operation.n_steps if the
    #   previous one is index-based, 100 if 100%-based.
    progress_occurred = qtc.pyqtSignal(
        str,  # task name
        str,  # current operation
        int,  # current operation index (1-based)
        int,  # total no. operations
        int,  # current progress of the current operation
        int   # total number of steps in this operation
        )

    name = ''
    progress_name = ''  # Used for progress bar label
    description = ''

    def __init__(self, device,  # pylint: disable=unused-argument
                 *args, timeout=None, **kwargs):
        """Initialize task.

        Parameters
        ----------
        device : CameraABC, or ControllerABC, or None
            The device to which this calibration task belongs.
        *args : object
            Unused positional arguments that subclasses may need.
        timeout : int, or None, optional
            The maximum duration for this task, in milliseconds.
            If negative or None, this task will not time out.
            Default is None.
        **kwargs : object, optional
            Other unused keyword arguments that subclasses may need.
        """
        super().__init__(parent=None)

        self.__device = device
        self.__flags = {'running': False,
                        'to_be_done': True,
                        'does_time_out': True,
                        'aborted': False,}

        timer = qtc.QTimer(self)
        timer.setSingleShot(True)
        timer.timeout.connect(self.on_timeout)
        if timeout is None or timeout < 0:
            self.__flags['does_time_out'] = False
        else:
            timer.setInterval(round(timeout))

        self._task_settings = None  # Device settings while running
        self._original = {'settings': None}  # Restored at the end

        # Prepare an _InfoBox instance. This instance is only
        # temporary, though. It will be made anew very soon in
        # the main GUI thread using the __dispatcher below (a
        # reference is kept to avoid the garbage collector)
        info_title = self.name
        if self.device:
            info_title += f" for {self.device.name}"
        info_slots = {'accepted': self.continue_,  # "I'm ready" button
                      'rejected': self.abort}      # "Abort"/"X" button
        self._info = _InfoBox(parent=None, title=info_title, slots=info_slots)
        self.__glob = {
            'timer': timer,
            'dispatch': base.QMainThreadDispatcher(self.__remake_info_box),
            }

        # Keep a list of (signal, slot_or_signal) that will be
        # connected by calling self._connect_device_signals()
        # Notice that the base implementation of CalibrationTask
        # DOES NOT CALL _connect_device_signals()
        self._to_be_connected = []
        base.safe_connect(self.error_occurred, self.abort,
                          type=_QUEUED|_UNIQUE)

        if not self.device:
            return

        self._original['settings'] = deepcopy(self.device.settings)
        self._task_settings = deepcopy(self.device.settings)
        self._to_be_connected.append(
            (self.device.error_occurred, self._on_device_error)
            )

    @property
    def device(self):
        """Return the device to which this task belongs."""
        return self.__device

    @property
    def to_be_done(self):
        """Return whether this task has been already performed or not."""
        return self.__flags['to_be_done']

    @property
    def is_aborted(self):
        """Return whether this task is in the process of being aborted."""
        return self.__flags['aborted']

    @property
    def is_running(self):
        """Return whether this task is currently being carried out."""
        return self.__flags['running']

    @property
    def original_settings(self):
        """Return the settings the device had before self was instantiated."""
        return self._original['settings']

    @qtc.pyqtSlot()
    def abort(self):
        """Abort this task.

        This slot can be extended in subclasses. This means that
        subclasses should call super().abort() at the end of their
        code. Subclasses should also make sure to return immediately
        if self.is_aborted, as this limits the risk of infinite
        recursion, should aborting occur as a result of an error that
        may repeat itself when executing abort().

        Emits
        -----
        aborted
            Right before returning.
        """
        if self.is_aborted:
            return
        self.restore_device()
        self.__flags['aborted'] = True
        self.__flags['running'] = False
        self.mark_as_done(False)
        self.aborted.emit()

    def _connect_device_signals(self):
        """Connect signals to slots. Use only when starting self."""
        for signal, slot in self._to_be_connected:
            base.safe_connect(signal, slot, type=_UNIQUE)

    @qtc.pyqtSlot()
    @abstractmethod
    def continue_(self, *_):
        """Continue this task after user confirmation."""

    def mark_as_done(self, is_done):
        """Mark this task as to be done or not."""
        if not is_done:
            self.__flags['to_be_done'] = True
            return
        self.__flags['running'] = False
        self.__flags['to_be_done'] = False
        for signal, slot in self._to_be_connected:
            base.safe_disconnect(signal, slot)
        self.done.emit()

    @qtc.pyqtSlot(tuple)
    def _on_device_error(self, error):
        """Respond to a device error.

        This slot can be extended in subclasses that may want to
        filter some of the device errors. The implementation in
        the base-class immediately re-emits self.error_occurred.

        Parameters
        ----------
        error : tuple or ViPErLEEDErrorEnum
            Error information. If a tuple, the first element
            is the error code, the second the error message.

        Returns
        -------
        silenced : bool
            Whether the error was silenced. Silenced errors do not
            cause emission of self.error_occurred

        Emits
        -----
        self.error_occurred(error)
            Unless it returns True.
        """
        self.error_occurred.emit(error)
        return False

    @qtc.pyqtSlot()
    def on_timeout(self):
        """React to a timeout event.

        This method must be overridden in all
        subclasses that can time out.

        Raises
        ------
        NotImplementedError
            If the method is not overridden
            in subclasses that can time out
        """
        if self.__flags['does_time_out']:
            raise NotImplementedError(
                f"{self.__class__.__name__} does time out, "
                "but on_timeout() was not overridden."
                )

    def restore_device(self):
        """Restore settings and other device attributes.

        This method is called each time the task is aborted
        or has finished. The original device settings can be
        accessed via self.original_settings.

        Subclasses can extend this method, and may modify the
        _original dictionary before calling the base-class
        implementation. The super().restore_device() call
        should be at the end of reimplementations.

        Returns
        -------
        None.
        """
        if self.original_settings is None:
            return
        _INVOKE(self.device, 'set_settings',
                qtc.Q_ARG(object, self.original_settings))

    def set_info_text(self, text):
        """Set informative text in the info message box."""
        _INVOKE(self._info, 'set_info_text', qtc.Q_ARG(str, text))

    @qtc.pyqtSlot()
    def show_info(self):
        """Show the information message box."""
        _INVOKE(self._info, 'exec')

    @qtc.pyqtSlot()
    @abstractmethod
    def start(self):
        """Start this task.

        This method should be extended in subclasses, i.e.,
        subclasses should check that super().start() returns
        True. Only in this case they can start the task.

        A simple implementation for a task that DOES NOT
        require any user confirmation can simply do:
        >>> ok_to_start = super().start()
        >>> if ok_to_start:
        >>>     self.continue_()
        >>> return ok_to_start

        A simple implementation for a task that DOES
        require user confirmation can simply do:
        >>> ok_to_start = super().start()
        >>> if ok_to_start:
        >>>     self.set_info_text("<text to be shown to the user>")
        >>>     self.show_info()
        >>> return ok_to_start

        In both cases, .continue_() should be extended to actually
        begin execution of the task at the hardware level.

        Returns
        -------
        ok_to_start : bool
            Whether starting the device is possible.

        Emits
        -----
        self.device.error_occurred(CalibrationTaskError.DEVICE_BUSY)
            If .start() is called while the device is busy.
        """
        if self.device is not None and self.device.busy:
            base.emit_error(self.device, CalibrationTaskError.DEVICE_BUSY,
                            self.name, self.device.name)
            return False

        self.update_device_settings()
        self.__flags['aborted'] = False
        self.__flags['running'] = True
        self.__flags['to_be_done'] = True
        if self.__flags['does_time_out']:
            self.__glob['timer'].start()
        return True

    def update_device_settings(self):
        """Update settings in device from _task_settings."""
        if self._task_settings is None:
            return
        if self._task_settings == self.device.settings:
            return
        _INVOKE(self.device, 'set_settings',
                qtc.Q_ARG(object, deepcopy(self._task_settings)))

    def __remake_info_box(self):
        """Overwrite self._info with an identical _InfoBox instance.

        This method is supposed to be used uniquely with
        QMainThreadDispatcher to recreate self._info with
        affinity to the main GUI thread. This is needed in
        case self was created in a non-GUI thread.

        Returns
        -------
        None.
        """
        self._info = _InfoBox.from_info_box(self._info)
