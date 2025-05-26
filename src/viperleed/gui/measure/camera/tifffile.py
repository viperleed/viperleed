"""Module tifffile of viperleed.gui.measure.camera.

This module defines the TiffImage class that is used to read and write
from/to disk uncompressed, gray-scale Tagged Images into/from numpy
arrays.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-11-15'
__license__ = 'GPLv3+'

from datetime import datetime
import getpass  # Used for storing PC user-name
from numbers import Number
from pathlib import Path
import struct

import numpy as np

from viperleed import __version__


class TiffTag:
    """Class that handles a single TIFF tag."""

    __field_types = ('byte', 'string', 'short', 'long', 'rational',
                     'signed_byte', 'undefined', 'signed_short',
                     'signed_long', 'signed_rational', 'float',
                     'double')

    __bytes_per_entry = {
        'byte': 1,
        'string': 1,
        'short': 2,
        'long': 4,
        'rational': 8,
        'signed_byte': 1,
        'undefined': 1,
        'signed_short': 2,
        'signed_long': 4,
        'signed_rational': 8,
        'float': 4,
        'double': 8
        }

    # Each id is 'Name': (number, type, count, default)
    # count and default are -1 if they depend on other tags,
    # they are None if no default exists.
    # See also: https://www.awaresystems.be/imaging/tiff/tifftags/
    __tag_info = {
        'NewSubfileType': (254, 'long', 1, 0),
        'SubfileType': (255, 'short', 1, None),
        'ImageWidth': (256, 'long', 1, None),
        'ImageLength': (257, 'long', 1, None),
        'BitsPerSample': (258, 'short', -1, 1),
        'Compression': (259, 'short', 1, 1),
        'PhotometricInterpretation': (262, 'short', 1, None),
        'Threshholding': (263, 'short', 1, 1),
        'DocumentName': (269, 'string', None, None),
        'ImageDescription': (270, 'string', None, None),
        'Make': (271, 'string', None, None),
        'Model': (272, 'string', None, None),
        'StripOffsets': (273, 'long', -1, None),
        'Orientation': (274, 'short', 1, 1),
        'SamplesPerPixel': (277, 'short', 1, 1),
        'RowsPerStrip': (278, 'long', 1, 2**32 - 1),
        'StripByteCounts': (279, 'long', -1, None),
        'XResolution': (282, 'rational', 1, None),
        'YResolution': (283, 'rational', 1, None),
        'PlanarConfiguration': (284, 'short', 1, 1),
        'XPosition': (286, 'rational', 1, None),
        'YPosition': (287, 'rational', 1, None),
        'GrayResponseUnit': (290, 'short', 1, 2),
        'GrayResponseCurve': (291, 'short', -1, None),
        'Group3Options': (292, 'long', 1, 0),
        'Group4Options': (293, 'long', 1, 0),
        'ResolutionUnit': (296, 'short', 1, 2),
        'PageNumber': (297, 'short', 2, None),
        'TransferFunction': (301, 'short', -1, -1),
        'Software': (305, 'string', None, None),
        'DateTime': (306, 'string', 20, None),
        'Artist': (315, 'string', None, None),
        'HostComputer': (316, 'string', None, None),
        'Predictor': (317, 'short', 1, 1),
        'WhitePoint/ColorImageType': (318, 'rational', 2, None),
        'PrimaryChromaticities/ColorList': (319, 'rational', 6, None),
        'ColorMap': (320, 'short', -1, None),
        'SampleFormat': (339, 'short', -1, 1),
        'ElectronEnergy': (34000, 'string', None, None),  # private!
        }

    __tag_ids = {k: v[0] for k, v in __tag_info.items()}
    __tag_names = {v[0]: k for k, v in __tag_info.items()}
    __tag_types = {v[0]: v[1] for v in __tag_info.values()}
    __tag_counts = {v[0]: v[2] for v in __tag_info.values()}
    __tag_defaults = {v[0]: (v[3],) for v in __tag_info.values()}

    def __init__(self, name='', index=-1, value=None, n_entries=0):
        """Initialize object."""
        if name and index >= 0:
            raise ValueError("Cannot create tag from both name and index.")
        if value is not None and not (name or index >=0):
            raise ValueError("Cannot set value of tag without "
                             "a name or an index.")
        # Attributes. Set further down.
        self.index = -1           # ordinal identifier
        self.__field_type = None  # type of the values
        self.n_entries = 0        # no. of values
        self.__value = None       # the value itself
        self.offset = None        # offset of value in Tiff file

        if index >= 0:
            self.index = index
        elif name:
            self.index = self.__tag_ids.get(name, None)
            if self.index is None:
                self.index = -1
                raise ValueError(f"Could not find tag name {name}.")
        else:
            self.index = -1

        self.__get_default_entries()  # self.n_entries is >= 0 now

        if self.n_entries and n_entries and n_entries != self.n_entries:
            raise ValueError(f"Number of entries given ({n_entries}) "
                             "inconsistent with those expected "
                             f"({self.n_entries})")
        if n_entries:
            self.n_entries = n_entries

        value = self.__check_value_consistency(value)
        if value is None or not value:
            value = self.__tag_defaults.get(self.index, (-1,))
            if value[0] is not None and value[0] < 0:
                value = None
        self.__value = value

    @property
    def has_value_outside(self):
        """Return whether value is stored outside of the tag."""
        return self.size > 4

    @property
    def field_type(self):
        """Return the field type of self."""
        if self.__field_type is None:
            raise RuntimeError(f"{self.__class__.__name__}: "
                               "field type was never set.")
        return self.__field_type

    @field_type.setter
    def field_type(self, new_field_type):
        """Set the field type of this tag."""
        if isinstance(new_field_type, str):
            if new_field_type not in self.__field_types:
                raise ValueError(f"{self.__class__.__name__}: Unknown "
                                 f"TiffTag field type {new_field_type}")
            self.__field_type = new_field_type
            return

        if not isinstance(new_field_type, int):
            raise ValueError(f"{self.__class__.__name__}: field "
                             "type can only be int or str.")

        try:
            self.__field_type = self.__field_types[new_field_type - 1]
        except IndexError as err:
            raise ValueError(
                f"{self.__class__.__name__}: Unknown "
                f"field type ordinal {new_field_type}"
                ) from err

    @property
    def item_size(self):
        """Return the size in bytes of each entry in self."""
        return self.__bytes_per_entry[self.field_type]

    @property
    def name(self):
        """Return the TiffTag name for self as a string."""
        return self.__tag_names.get(self.index, str(self.index))

    @property
    def size(self):
        """Return the total size of self in bytes."""
        return self.n_entries * self.item_size

    @property
    def value(self):
        """Return the value of this TiffTag.

        Returns
        -------
        value : Sequence
            The value of this TiffTag. If only one value,
            it is returned as the first element of a tuple.

        Raises
        ------
        RuntimeError
            If this property is accessed before the value of
            this tag was not set nor read from bytes before.
        """
        if self.__value is None:
            raise RuntimeError(f"{self.__class__.__name__}: value "
                               "was never set nor read from bytes.")
        return self.__value

    @value.setter
    def value(self, new_value):
        """Set the value (i.e., data) of this TiffTag."""
        new_value = self.__check_value_consistency(new_value)
        self.__value = new_value

    def as_bytes(self, byte_order):
        """Return self as bytes with a given byte_order."""
        if self.index < 0:
            raise RuntimeError("Unknown tag cannot be formatted as bytes")
        if self.has_value_outside and self.offset is None:
            raise RuntimeError("Tag has data outside, but no offset was given")

        tag_id = self.index.to_bytes(2, byte_order)
        field_type_ordinal = self.__field_types.index(self.field_type) + 1
        field_type = field_type_ordinal.to_bytes(2, byte_order)
        n_entries = self.n_entries.to_bytes(4, byte_order)

        if self.has_value_outside:
            # The last 4 bytes are the offset
            data_or_offset = self.offset.to_bytes(4, byte_order)
        else:
            # The last 4 bytes are the value itself,
            # left-justified into 4 bytes if the data
            # type is less than 4-bytes long
            data_or_offset = self.value_as_bytes(byte_order).ljust(4, b'\x00')

        return tag_id + field_type + n_entries + data_or_offset

    def value_as_bytes(self, byte_order):
        """Return the value of self as bytes."""
        if any(v is None for v in self.value):
            # print(f"{self.name=}, {self.n_entries=}, {self.value=}")
            raise RuntimeError(f"No valid value for Tag {self.name}")
        byte_fmt = '>' if byte_order == 'big' else '<'
        bytes_value = b''
        if self.field_type == 'string':
            bytes_value = ''.join(self.value).encode('utf-8') + b'\x00'
        elif self.field_type == 'float':
            bytes_value = struct.pack(
                f"{byte_fmt}{self.n_entries}f{self.item_size}",
                *self.value
                )
        elif self.field_type == 'double':
            bytes_value = struct.pack(
                f"{byte_fmt}{self.n_entries}d{self.item_size}",
                *self.value
                )
        elif self.field_type in ('byte', 'signed_byte'):
            bytes_value = b''.join(self.value)
        else:  # Better with struct.pack??
            for value in self.value:
                if self.field_type in ('rational', 'signed_rational'):
                    bytes_value += value.to_bytes(self.item_size // 2,
                                                  byte_order)
                else:
                    bytes_value += value.to_bytes(self.item_size, byte_order)

        return bytes_value

    def value_from_bytes(self, value_as_bytes, byte_order):
        """Calculate the value of self from its bytes."""
        # Decode the value depending on the field type.
        # I'm actually not 100% sure that I'm treating signed
        # and unsigned quantities correctly. Probably would
        # be better to actually use struct.unpack for all?
        byte_fmt = '>' if byte_order == 'big' else '<'
        if self.field_type == 'string':
            # Tag value is a null-terminated string.
            # Skip the b'\0' at the end and decode
            value = value_as_bytes[:-1].decode('utf-8')
        elif self.field_type == 'float':
            value = struct.unpack(
                f"{byte_fmt}{self.n_entries}f{self.item_size}",
                value_as_bytes
                )
        elif self.field_type == 'double':
            value = struct.unpack(
                f"{byte_fmt}{self.n_entries}d{self.item_size}",
                value_as_bytes
                )
        elif self.field_type in ('byte', 'signed_byte'):
            # Tag value is bytes
            value = value_as_bytes
        else:
            # Tag is a series of integer values
            values = (value_as_bytes[i:i+self.item_size]
                      for i in range(0, self.size, self.item_size))
            if self.field_type in ('rational', 'signed_rational'):
                # Tag is a series of numerator/denominator
                # ratios in the most- and least-significant
                # four bytes of each of the values
                value = tuple(int.from_bytes(v[:4], byte_order)
                              / int.from_bytes(v[4:], byte_order)
                              for v in values)
            else:
                # some sort of integer
                value = tuple(int.from_bytes(v, byte_order) for v in values)

        # Finally, make value a single element, if there
        # is only one, unless it is a tag with info on
        # the arrangement of image data
        if (len(value) == 1
                and not self.name in ('StripOffsets', 'StripByteCounts')):
            value = value[0]
        self.__value = value

    def __get_default_entries(self):
        """Set .n_entries and .field_type to the default."""
        if self.index >= 0:
            self.field_type = self.__tag_types[self.index]
            n_entries = self.__tag_counts[self.index]
            if n_entries is not None and n_entries > 0:
                self.n_entries = n_entries

    def __check_value_consistency(self, value):
        """Check that value fits with the expected number of entries.

        Parameters
        ----------
        value : object
            The value to be checked

        Returns
        -------
        adjusted_value : iterable
            Always returns an iterable

        Raises
        ------
        ValueError
            If the number of "elements" in value does not match
            what is expected for this tag.
        """
        if value is None or not value:
            return tuple()

        # Check that the value passed fits
        # with the number of entries expected
        if not hasattr(value, '__len__'):
            value = (value,)
        n_values = len(value)

        n_expected = self.n_entries
        if self.field_type == 'string':
            n_values += 1  # The '\x00' terminator
        elif self.field_type in ('rational', 'signed_rational'):
            n_expected *= 2  # Numerator & denominator

        if n_expected <= 0:
            n_expected = n_values

        if n_values != n_expected:
            raise ValueError(f"Invalid number of values for tag {self.name}. "
                             f"Expected {n_expected}, found {n_values}.")
        self.n_entries = n_expected
        return value


class TiffFile:
    """Class for reading/writing gray-scale TIFF images."""

    def __init__(self):
        """Initialize TiffImage object.

        Use the class methods .read(filename) or
        .from_array(array) to fill in the data.

        Returns
        -------
        None.
        """
        self.__filename = ''
        self.__images = []
        self.__image_info = []
        self.__tiff_file = None
        self.__byte_order = 'big'

    @classmethod
    def from_array(cls, array, **infos):
        """Return a TiffImage object from an array.

        Parameters
        ----------
        array : numpy.ndarray
            Image data. If array is 2D, it will be interpreted
            as an image, if it is 3D, it will be interpreted as
            a sequence of images (i.e., an image stack). In this
            case, array[i] will be taken as the i-th image.
        **infos : dict
            Information to be stored together with the images.
            If an image stack is given, each key should be a
            sequence of information (one per image). Should
            some of the information is missing, None can be used.
            The following keys will be looked for:
            'make' : Sequence of str
                Manufacturer of the camera that acquired the image
            'model' : Sequence of str
                Camera model
            'energy' : Sequence of float or str
                Electron energy (in electron-volts) for this LEED
                image.
            'date_time' : Sequence of datetime or str
                Date and time of image acquisition. If a string,
                it should be formatted as "YYYY:MM:DD HH:MM:SS".
            'comment' : Sequence of str
                Extra comments to be included in the image
                (e.g., sample description)

        Returns
        -------
        tiff_file : TiffFile
            The TiffFile object with pixel data from array

        Raises
        ------
        TypeError
            If array data type is invalid
        ValueError
            If array does not have exactly two dimensions
        """
        array = np.asarray(array)
        dtype = array.dtype
        if dtype.kind not in 'ifub':
            raise TypeError(f"{cls.__name__}: invalid data type {dtype} "
                            "for TIFF data. Expected integer or float.")
        if len(array.shape) < 2 or len(array.shape) > 3:
            raise ValueError(f"{cls.__name__}: array with {len(array.shape)} "
                             "dimensions is not a valid gray-scale image.")
        self = cls()
        if len(array.shape) == 3:
            self.__images = array
        else:
            self.__images = [array]
        self.__set_extra_info(**infos)
        return self

    @classmethod
    def read(cls, filename):
        """Read image from a file.

        Parameters
        ----------
        filename : str or Path
            Path to file to be read. Should have a .tif
            or .tiff extension.

        Returns
        -------
        tiff_image : TiffImage
            TiffImage object with the file contents read in.

        Raises
        ------
        FileNotFoundError
            If filename does not point to an existing file.
        ValueError
            If filename has incorrect extension.
        """
        filename = Path(filename)

        if not filename.is_file():
            raise FileNotFoundError(
                f"{cls.__name__}: No file {filename.resolve()} found"
                )
        if not filename.suffix.startswith('.tif'):
            raise ValueError(
                f"{cls.__name__}: Invalid file extension {filename.suffix!r}"
                )
        self = cls()
        self.__filename = filename

        with open(filename, 'rb') as tiff_file:
            byte_order = tiff_file.read(2)
            self.__byte_order = 'big' if byte_order == b'MM' else 'little'
            tiff_id = int.from_bytes(tiff_file.read(2), self.__byte_order)
            if tiff_id != 42:
                raise ValueError(
                    f"{cls.__name__}: {filename.resolve()} "
                    "is not a valid TIFF file."
                    )
            self.__tiff_file = tiff_file
            # Pointer to first image-file directory (IFD):
            next_ifd = int.from_bytes(self.__tiff_file.read(4),
                                      self.__byte_order)
            while next_ifd:
                tiff_file.seek(next_ifd)
                next_ifd = self.__read_ifd()
            self.__tiff_file = None

        if not self.__images:
            raise ValueError(f"{cls.__name__}: {filename.resolve()} "
                             "does not contain any images.")
        return self

    @property
    def images(self):
        """Return images as an iterable of numpy.ndarrays."""
        return self.__images

    def __set_extra_info(self, **infos):
        """Set optional images information from a dictionary."""
        self.__check_extra_infos(infos)  # Raises if not OK
        self.__image_info = [None]*len(self.images)

        zzip = zip(
            ('make', 'model', 'energy', 'date_time', 'comment'),
            ('Make', 'Model', 'ElectronEnergy', 'DateTime', 'ImageDescription')
            )

        for key, tag_name in zzip:
            # Pick an appropriate string-converter
            # pylint: disable=redefined-variable-type
            _to_str = str
            if key == 'date_time':
                _to_str = self.__datetime_to_str
            elif key == 'energy':
                _to_str = self.__energy_to_str

            # Pack converted values into __image_info
            values = infos.get(key, tuple())
            for i, value in enumerate(_to_str(v) for v in values):
                if not self.__image_info[i]:
                    self.__image_info[i] = {}
                if value != 'None':
                    self.__image_info[i][tag_name] = value

        # And add the user name of this PC
        user = getpass.getuser()
        for i, info in enumerate(self.__image_info):
            if info is None:
                info = {}
                self.__image_info[i] = info
            info['HostComputer'] = user

    def __check_extra_infos(self, infos):
        """Check consistency of infos with number of images."""
        # Check that no. of entries == no. images
        # while allowing also to pass single values
        # if there is ony one image to be saved.
        n_images = len(self.images)
        for key, val in infos.items():
            if isinstance(val, (str, Number, datetime)):
                infos[key] = (val,)
            # Check number of entries == number of images
            n_entries = len(infos[key])
            if n_entries != n_images:
                raise ValueError(f"Not enough entries ({n_entries}) "
                                 f"in info {key}. Expected {n_images}")

    @staticmethod
    def __datetime_to_str(date_time):
        """Convert a date-time-like object to string."""
        if not date_time:
            return 'None'

        # Expect the format "YYYY:MM:DD HH:MM:SS"
        if isinstance(date_time, datetime):
            date_time = date_time.strftime("%Y:%m:%d %H:%M:%S")

        if (len(date_time) != 19
                or not all(date_time[i] == ':' for i in (4, 7, 13, 16))
                or date_time[10] != ' '):
            raise ValueError(
                f"Invalid date_time info {date_time}. "
                "Requires format 'YYYY:MM:DD HH:MM:SS'. "
                "Use strftime('%Y:%m:%d %H:%M:%S')"
                )
        return date_time

    @staticmethod
    def __energy_to_str(energy):
        """Convert an energy value to the expected string format."""
        if energy in (None, 'None', ''):
            return 'None'
        if isinstance(energy, Number):
            return f"{energy:.1f} eV"
        return str(energy)

    def info(self, img_index):
        """Return the info of an image.

        Returns
        -------
        info : dict
            keys are tag names, values are tag values.

        Raises
        ------
        RuntimeError
            If this method is called before the file was saved.
        """
        try:
            info = self.__image_info[img_index]
        except IndexError as err:
            raise RuntimeError("Image info is not available "
                               "until the file is saved.") from err
        if info is None:
            return {}
        return info

    def write(self, filepath):
        """Write self.images to filepath."""
        filepath = Path(filepath).resolve()
        if not filepath.parent.exists():
            filepath.parent.mkdir(parents=True)
        if not self.__image_info:
            self.__image_info = [None]*len(self.images)

        with open(filepath, 'wb') as tiff_file:
            self.__tiff_file = tiff_file

            # Begin with the header
            b_order = b'MM' if self.__byte_order == 'big' else b'II'
            tiff_id = (42).to_bytes(2, self.__byte_order)
            # First IFD offset: b_oder (==2) + tiff_id (==2) +
            # ifd_offset (==4)
            ifd_offset = (8).to_bytes(4, self.__byte_order)

            # Write down the header:
            self.__tiff_file.write(b_order + tiff_id + ifd_offset)

            # Then write one IFD per image.
            for i, _ in enumerate(self.images):
                self.__write_ifd(i)

    def __read_ifd(self):
        """Read one image-file directory.

        Images found are stored in self.__images.

        Returns
        -------
        next_ifd_start : int
            Offset in bytes from the beginning of the file
            at which the next image-file directory starts.
            Will be zero if this was the last IFD.

        Raises
        ------
        ValueError
            If the file from which the IFD is to be read does
            not contain the minimum set of tags that allows
            interpretation as a TIFF grayscale image.
        RuntimeError
            If an unknown "SampleFormat" tag value (i.e., data
            type per pixel) is found.
        """
        n_tags = int.from_bytes(self.__tiff_file.read(2), self.__byte_order)
        tags = dict(self.__read_tag() for _ in range(n_tags))

        next_ifd_start = int.from_bytes(self.__tiff_file.read(4),
                                        self.__byte_order)

        # Now read the image data for this IFD: We
        # expect ImageWidth * ImageLength entries of
        # BitsPerSample/8 length (in bytes), arranged
        # into strips of RowsPerStrip number of rows,
        # each starting at one of the StripOffsets
        if not all(tag in tags for tag in ('ImageWidth', 'ImageLength',
                                           'BitsPerSample', 'StripByteCounts',
                                           'StripOffsets',
                                           'PhotometricInterpretation')):
            raise ValueError(
                f"{self.__class__.__name__}: file {self.__filename.resolve()} "
                "is missing data for interpreting a gray-scale image."
                )

        # Prepare the correct data type
        bytes_per_px = int(tags['BitsPerSample']//8)
        order = '>' if self.__byte_order == 'big' else '<'
        if 'SampleFormat' not in tags or tags['SampleFormat'] == 1:
            px_type = 'u'
        elif tags['SampleFormat'] == 2:
            px_type = 'i'
        elif tags['SampleFormat'] == 3:
            px_type = 'f'
        else:
            raise RuntimeError(
                f"Unknown pixel format id {tags['SampleFormat']}"
                )

        dtype = f"{order}{px_type}{bytes_per_px}"
        pixel_intensity = []
        for b_count, offset in zip(tags['StripByteCounts'],
                                   tags['StripOffsets']):
            self.__tiff_file.seek(offset)
            strip = self.__tiff_file.read(b_count)
            pixel_intensity.append(
                np.frombuffer(strip, dtype=dtype)
                )

        image = np.concatenate(pixel_intensity).reshape(
            tags['ImageLength'],
            tags['ImageWidth']
            )

        # pylint: disable=compare-to-zero
        # Disable as it can only be zero or 1, and
        # explicit comparison feels clearer here.
        if tags['PhotometricInterpretation'] == 0:
            image = 2**tags['BitsPerSample'] - 1 - image

        self.__images.append(image)
        self.__image_info.append(tags)

        return next_ifd_start

    def __read_tag(self):
        """Read and decode the info in a TIFF tag.

        The tag at the current file position is read.

        Returns
        -------
        tag_id : str
            One of the supported TIFF tag ids
        tag_value : object
            The value of the tag.
        """
        # Each tag is:
        # 2 bytes for numeric identifier
        # 2 bytes for field type (int/string/etc...)
        # 4 bytes for number of entries of the type in the tag
        # 4 bytes for the absolute offset (in bytes) of the value

        tag = TiffTag()
        tag.index = int.from_bytes(self.__tiff_file.read(2),
                                   self.__byte_order)
        tag.field_type = int.from_bytes(self.__tiff_file.read(2),
                                        self.__byte_order)
        tag.n_entries = int.from_bytes(self.__tiff_file.read(4),
                                       self.__byte_order)
        offset = self.__tiff_file.read(4)

        if tag.size <= 4:
            # Value is not at offset, but stored
            # directly in the offset bytes
            value_as_bytes = offset[:tag.size]
        else:
            # Look at the value pointed to by the offset,
            # storing the current position in the file
            # beforehand.
            old_position = self.__tiff_file.tell()
            self.__tiff_file.seek(int.from_bytes(offset, self.__byte_order))
            value_as_bytes = self.__tiff_file.read(tag.size)
            self.__tiff_file.seek(old_position)

        tag.value_from_bytes(value_as_bytes, self.__byte_order)

        return tag.name, tag.value

    def __write_ifd(self, img_number):
        """Write an image to file.

        Parameters
        ----------
        img_number : int
            Index in self.images of the image to write.

        Returns
        -------
        None
        """
        # Prepare the tags for this image, including
        # extra information that may be present
        image = self.images[img_number]
        try:
            extras = self.info(img_number)
        except RuntimeError:
            # No extra info for this image
            extras = {}

        (tags,
         strip_offsets,
         total_n_bytes) = self.__make_tags_for_image(image, extras)

        # Here begins the 'header' for this IFD
        # (H.1) number of tags
        n_tags = len(tags)
        self.__tiff_file.write(n_tags.to_bytes(2, self.__byte_order))

        # Prepare to later write the data of each tag. We start after
        # the end of this IFD 'header', i.e., current position + 12
        # bytes per tag, + 4 bytes for next IFD offset
        curr_pos = self.__tiff_file.tell()
        tag_data_offset = curr_pos + 12*n_tags + 4

        # tag_data_offset needs to be even:
        tag_data_offset = ((tag_data_offset + 1) // 2) * 2

        # Set up the data offsets for each tag
        for tag in tags:
            if tag.has_value_outside:
                tag.offset = tag_data_offset
                tag_data_offset += tag.size
                tag_data_offset = ((tag_data_offset + 1) // 2) * 2

        # Now update StripOffsets. This will directly point to the
        # first image data byte, as we use only one strip. It will
        # go right after the tag data.
        strip_offsets.value = tag_data_offset

        # (H.2) Now write all tags...
        for tag in tags:
            self.__tiff_file.write(tag.as_bytes(self.__byte_order))

        # ...store them as image info...
        self.__image_info[img_number] = {tag.name: tag.value for tag in tags}

        # (H.3)...and write the next IFD offset.
        # This concludes the IFD 'header'.
        if img_number == len(self.images) - 1:
            # Last image. Write four zeros to signal the end.
            next_ifd_offset = b'\x00\x00\x00\x00'
        else:
            next_ifd_offset = tag_data_offset + total_n_bytes
            next_ifd_offset = next_ifd_offset.to_bytes(4, self.__byte_order)
        self.__tiff_file.write(next_ifd_offset)

        # Finally, tag data...
        for tag in tags:
            if tag.has_value_outside:
                self.__tiff_file.seek(tag.offset)
                self.__tiff_file.write(tag.value_as_bytes(self.__byte_order))

        # ...and the actual pixels, after
        # converting to the right byte order
        new_dtype = image.dtype.newbyteorder(self.__byte_order)
        if new_dtype != image.dtype:
            image = image.astype(new_dtype)
        self.__tiff_file.seek(strip_offsets.value[0])
        self.__tiff_file.write(image.tobytes())

    def __make_tags_for_image(self, image, extras):
        """Return a list of tags for an image.

        Extra informaiton are also included.

        Parameters
        ----------
        image : numpy.ndarray
            The image for which tags should be returned
        extras : dict
            Extra information for this image in the form
            {<tag_name>: <tag_value>}.

        Returns
        -------
        tags : list
            Sorted list of tags that describe this image.
        strip_offsets : TiffTag
            The "StripOffsets" tag that requires updating while
            other tags are being written to file.
        total_n_bytes : int
            Total number of bytes that this image takes.

        Raises
        -------
        ValueError
            If image has an invalid data type. Acceptable data
            types are signed and unsigned integers, bytes, and
            float.
        """
        # Read some info from the image:
        dtype = image.dtype
        if not dtype.kind in 'ifub':
            raise ValueError(f"{self.__class__.__name__}: Cannot save image "
                             f" with data type {dtype.kind!r}. Only "
                             " integer and float data types are allowed.")
        height, width = image.shape
        n_bytes_px = dtype.itemsize
        total_n_bytes = width*height*n_bytes_px
        if dtype.kind == 'i':
            pixel_type = 2
        elif dtype.kind == 'f':
            pixel_type = 3
        else:  # unsigned integer
            pixel_type = 1

        # TODO: it appears to be better for readers to store
        # images in multiple strips, each roughly 8Kbytes long.
        # For now, store the whole image in a single strip, as
        # it is much simpler. Multiple strips require setting
        # a smaller RowsPerStrip, and storing appropriate
        # StripOffsets and StripByteCounts.

        # Now prepare tags:
        strip_offsets = TiffTag(name='StripOffsets',
                                n_entries=1)  # Value is set later
        tags = [
            TiffTag(name='ImageWidth', value=width),
            TiffTag(name='ImageLength', value=height),
            TiffTag(name='BitsPerSample', value=8*n_bytes_px),
            TiffTag(name='PhotometricInterpretation', value=1),
            strip_offsets,
            TiffTag(name='RowsPerStrip', value=height),
            TiffTag(name='StripByteCounts', value=total_n_bytes),
            TiffTag(name='SampleFormat', value=pixel_type),
            TiffTag(name='Software', value=f"ViPErLEED v{__version__}"),
            ]

        for name, value in extras.items():
            tags.append(TiffTag(name=name, value=value))

        # Tags should be written in increasing ordinal number
        tags.sort(key=lambda tag: tag.index)

        return tags, strip_offsets, total_n_bytes
