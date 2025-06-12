"""Module imageprocess of viperleed.gui.measure.camera.

This module defines the ImageProcessor class that is used to process
single camera frames into a final image. The following processing
steps take place:
    (1) Summing frames
    (2) Application of a region of interest
    (3) Removal of bad pixels
    (4) Binning
    (5) Averaging (i.e., sum/counts)
Steps 2 and 4 are skipped if the camera supports them as a native
hardware feature; steps 1 and 5 are also skipped if the camera
supports frame averaging at the hardware level.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    'Florian Dörr (@FlorianDoerr)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-07-16'
__license__ = 'GPLv3+'

import copy
from dataclasses import dataclass, field
import typing
from pathlib import Path

import numpy as np
from numpy.lib.stride_tricks import as_strided
from PyQt5 import QtCore as qtc

from viperleed.gui.measure.camera import tifffile as tiff


@dataclass
class ImageProcessInfo:  # pylint: disable=too-many-instance-attributes
    """Data class storing information needed for processing images.

    Attributes
    ----------
    bad_pixels : BadPixels
        BadPixels object holding information about which of the camera
        pixels are bad and how to correct them.
    filename : str
        Name of the file to be saved to disk after processing is over.
        If not given, no image will be saved after processing. It must
        have a '.tiff' suffix.
    base_path : str, optional
        Path to which filename will be saved. If not given, the current
        directory is used.
    n_frames : int
        The number of frames to be averaged. The truth value of this
        attribute should be False if the camera handles frame averaging
        internally
    roi : Sequence
        Region of interest. Four elements: x_offset, y_offset, width,
        height. Offsets are with respect to the top-left corner. The
        truth value of this attribute should be False if the camera
        handles the ROI internally
    binning : int
        Binning factor. Binning is executed on a (binning x binning)
        mesh. If the size of frames is not a integer multiple of
        binning, the bottom-right pixels will be removed. The truth
        value of this attribute should be False if the camera handles
        binning internally
    frame_times : list
        List of frame-arrival times.
    camera : CameraABC
        The camera to which these image processing infos are related.
    energy : str or float, optional
        The current LEED energy used, in electron-volts. Only used for
        storing information.
    date_time : datetime, optional
        The datetime.now() object for the current image.
    comment : str, optional
        Any user comment that will be stored with the image (e.g.,
        sample surface, preparation details, etc.).
    count : int, optional
        This is used exclusively if filename contains "{__count__".
        It can be used to progressively number frames. It is used
        to produce a filename like follows:
        filename = filename.format(__count__=count)

    Methods
    ------
    copy()
        Return a deep-copy of self
    """

    __filename: str = ''
    base_path: str = ''
    n_frames: int = 0
    roi: typing.Sequence = tuple()
    binning: int = 0
    frame_times: typing.Sequence = field(default_factory=list)
    camera: typing.Any = None     # the camera object
    energy: typing.Any = None     # energy in eV, typically float
    date_time: typing.Any = None  # datetime object
    comment: str = ''
    count: int = 0

    @property
    def bad_pixels(self):
        """Return bad pixel info as a BadPixels object."""
        return self.camera.bad_pixels

    @property
    def filename(self):
        """Return the name of the file to save to."""
        fname = self.__filename
        if "{__count" in fname:
            fname = fname.format(__count__=self.count)
        return fname

    @filename.setter
    def filename(self, new_name):
        """Set the name of the file to save to."""
        self.__filename = new_name

    def copy(self):
        """Return a deep-copy of self."""
        return copy.deepcopy(self)

    def clear_times(self):
        """Clear the frame arrival times."""
        self.frame_times = []

    def restore_from(self, other):
        """Restore this instance from another one."""
        if not isinstance(other, ImageProcessInfo):
            raise ValueError(f"{self.__class__.__name__}: Cannot restore "
                             f"from a {other.__class__.__name__} instance")
        for attr, value in vars(other).items():
            if not hasattr(self, attr):
                continue
            setattr(self, attr, value)


class ImageProcessor(qtc.QObject):
    """Class that processes frames."""

    all_frames_acquired = qtc.pyqtSignal()
    image_processed = qtc.pyqtSignal(np.ndarray)
    image_saved = qtc.pyqtSignal(str)  # path to file

    def __init__(self):
        """Initialize the instance."""
        super().__init__()
        self.process_info = ImageProcessInfo()
        self.processed_image = np.zeros(0)
        self.n_frames_received = 0
        self.frame_bits = 16
        self.busy = False

        # Prevent double-processing of images (Issue #122).
        # _should_process_image is made True in prepare_to_process.
        self._should_process_image = False

    @property
    def missing_frames(self):
        """Return the number of frames that were not yet received."""
        return self.n_frames_to_process - self.n_frames_received

    @property
    def n_frames_to_process(self):
        """Return the number of frames to be processed.

        No image will be saved until n_frames_to_process
        frames have been accumulated.

        Returns
        -------
        n_frames_to_process : int
            The number of frames to process
        """
        return self.process_info.n_frames

    def prepare_to_process(self, process_info, first_frame):
        """Prepare the necessary attributes for processing.

        This method should be called before moving the ImageProcessor
        to a new thread, where the calculations will actually occur.

        Parameters
        ----------
        process_info : ImageProcessInfo
            See documentation of ImageProcessInfo for details.
        first_frame : numpy.ndarray
            The first frame. No processing is done, but it is
            used to properly set up the size and data type of
            the array that will contain the processed image.
        """
        self.busy = True
        self.process_info = process_info
        self.n_frames_received = 0
        self._should_process_image = True

        # Select an appropriate data type to avoid overflowing
        # [kind is either 'u' (unsigned) or '' (signed)]
        kind, frame_bits = first_frame.dtype.name.split('int')
        self.frame_bits = int(frame_bits)

        process_bits = 64 if self.frame_bits > 16 else 32

        self.processed_image = np.zeros_like(first_frame,
                                             dtype=f"{kind}int{process_bits}")

    @qtc.pyqtSlot(np.ndarray)
    def process_frame(self, frame):
        """Process a new frame."""
        if self.missing_frames > 0:
            # Add the new frame
            self.processed_image += frame
            self.n_frames_received += 1

        if self.missing_frames > 0:
            return

        if self._should_process_image:  # All frames arrived
            self._should_process_image = False
            self.all_frames_acquired.emit()
            self.apply_roi()
            self.remove_bad_pixels()
            self.bin_and_average()
            fname = self.save()
            self.image_saved.emit(fname)

    def remove_bad_pixels(self):
        """Remove bad pixels by neighbor averaging."""
        bad = self.process_info.bad_pixels
        if not bad:
            return
        bad_coords = bad.bad_pixel_coordinates.T
        repl_offsets = bad.replacement_offsets.T

        # In the following, the conversion to tuple when indexing
        # is critical, as it changes the indexing behavior to the
        # correct one: it is like doing (x1+x2, y1+y2). Without it,
        # indexing of an array A with another one (ind) would return
        # as many arrays as there are indices in ind[0], with
        # return[i] == A[ind[0][i]] (i.e., the rows of A with indices
        # ind[0][i]).
        replacements = (
            self.processed_image[tuple(bad_coords + repl_offsets)]
            + self.processed_image[tuple(bad_coords - repl_offsets)]
            ) >> 1  # right shift by 1 bit is division by 2
        self.processed_image[tuple(bad_coords)] = replacements

    def apply_roi(self):
        """Crop the image to a region of interest."""
        roi = self.process_info.roi
        if not roi:
            # ROI is taken care of by the camera itself
            return

        roi_x, roi_y, roi_w, roi_h = roi
        self.processed_image = self.processed_image[roi_y:roi_y+roi_h,
                                                    roi_x:roi_x+roi_w]

    def bin_and_average(self):
        """Apply binning and frame averaging to the processed image."""
        bin_fact = self.process_info.binning

        if not bin_fact or bin_fact == 1:
            # Binning is done in the camera, or
            # no actual binning should be done.
            # Do only averaging, if necessary.
            if self.n_frames_to_process == 1:
                return
            averaging_fact = 1 / self.n_frames_to_process
            self.processed_image = self.processed_image * averaging_fact
            return

        # Binning factors should be compatible with ROI settings,
        # except when no ROI is applied at all
        img_w, img_h = self.processed_image.shape
        if img_w % bin_fact or img_h % bin_fact:
            img_w = (img_w // bin_fact)*bin_fact
            img_h = (img_h // bin_fact)*bin_fact
            self.processed_image = self.processed_image[:img_w, :img_h]

        # The following is a fast binning code coming from
        # https://stackoverflow.com/questions/28906298
        strided_shape = (img_w//bin_fact, img_h//bin_fact, bin_fact, bin_fact)
        strided_strides = ((self.processed_image.strides[0]*bin_fact,
                            self.processed_image.strides[1]*bin_fact)
                           + self.processed_image.strides)
        strided_view = as_strided(self.processed_image,
                                  shape=strided_shape,
                                  strides=strided_strides)

        averaging_fact = 1 / (self.n_frames_to_process * bin_fact**2)

        self.processed_image = strided_view.sum(axis=(-2,-1))*averaging_fact

    def save(self):
        """Save image to disc as uncompressed TIFF.

        The image is saved only if there is a filename in the
        ImageProcessInfo stored in self.process_info.

        The TIFF data type is 16 or 32 bits depending on what
        was the original size of the camera frames. The data
        is stored in big-endian byte order.

        Returns
        -------
        full_path : str
            Path to file saved. An empty string is returned
            if no file was saved (because no filename was given).

        Raises
        ------
        ValueError
            If this method is called while self.process_info
            contains a non-empty .filename that des not have
            a .tif* extension.
        """
        dtype = '>u4' if self.frame_bits > 16 else '>u2'

        data = self.processed_image.astype(dtype)
        self.image_processed.emit(data)
        fname = self.process_info.filename
        if not fname:
            self.busy = False
            return ""

        fname = Path(fname)
        if self.process_info.base_path:
            fname = self.process_info.base_path / fname
        if '.tif' not in fname.suffix:
            raise ValueError("Can only save .tiff files. Found invalid "
                             f"extension {fname.suffix}")

        cam = self.process_info.camera
        bin_factor = self.process_info.binning
        comment = (f"Average of {self.process_info.n_frames} frames; "
                   f"Exposure: {cam.exposure} ms; "
                   f"Gain: {10**(cam.gain/20):.1f} ({cam.gain:.1f} dB); "
                   f"Binning: {bin_factor}x{bin_factor} pixels; ")

        info = {'model': cam.name, 'comment': comment,}
        if self.process_info.comment:
            info['comment'] = self.process_info.comment + '\n' + comment
        if self.process_info.energy:
            info['energy'] = self.process_info.energy
        if self.process_info.date_time:
            info['date_time'] = self.process_info.date_time

        img = tiff.TiffFile.from_array(data, **info)
        img.write(fname)
        self.busy = False
        return str(fname.resolve())
