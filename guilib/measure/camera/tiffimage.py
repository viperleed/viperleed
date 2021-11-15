"""Module tiffimage of viperleed.guilib.measure.camera.

========================================
   ViPErLEED Graphical User Interface
========================================

Created: 2021-11-15
Author: Michele Riva

This module defines the TiffImage class that is used to read and write
from/to disk uncompressed, gray-scale Tagged Images into/from numpy arrays.
"""
from dataclasses import dataclass
from pathlib import Path
import struct
import typing

import numpy as np


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

    __tag_ids = {
        254: 'NewSubfileType',             # NO
        255: 'SubfileType',                # NO
        256: 'ImageWidth',                 # b'\x01\x00' + b'\x00\x04' (may also be b'\x00\x03') + b'\x00\x00\x00\x01' + bytes(n_cols)
        257: 'ImageLength',  # = height    # b'\x01\x01' + b'\x00\x04' (may also be b'\x00\x03') + b'\x00\x00\x00\x01' + bytes(n_rows)
        258: 'BitsPerSample',              # b'\x01\x02\x00\x03\x00\x00\x00\x01' + bytes(n_bits)
        259: 'Compression',                # NO
        262: 'PhotometricInterpretation',  # b'\x01\x06\x00\x03\x00\x00\x00\x01\x00\x01\x00\x00' (black is zero)
        263: 'Thresholding',               # NO
        269: 'DocumentName',               # NO
        270: 'ImageDescription',           # NO, perhaps later
        271: 'Make',                       # NO, perhaps later
        272: 'Model',                      # NO, perhaps later
        273: 'StripOffsets',               # YES: where do the data strips start (wrt beginning of file)
                                           #      b'\x01\x11' + (word [b'\x00\x03'] or dword [b'\x00\x04']) + (depends on n_strips, should be b'\x00\x00\x00\x01' if we don't split images) + bytes(offset address)
                                           #      offset address needs to be contiguous when no compression, at the end of the tags
        274: 'Orientation',                # NO (defaults to 1st row = top, 1st col = left)
        277: 'SamplePerPixel',             # NO
        278: 'RowsPerStrip',               # in principle not necessary, although it should be < 8K if one wants to speed up access to the file. e.g., n_rows can be used:
                                           # b'\x01\x16' + b'\x00\x03' (may also be b'\x00\x04') + b'\x00\x00\x00\x01' + bytes(n_rows)
        279: 'StripByteCounts',            # b'\x01\x17' + b'\x00\x04' (may also be b'\x00\x03') + b'\x00\x00\x00\x01' + bytes(n_cols*rows_per_strip(=n_rows)*n_bits/8)
        282: 'XResolution',                # NO
        283: 'YResolution',                # NO
        284: 'PlanarConfiguration',        # NO
        286: 'XPosition',                  # NO
        287: 'YPosition',                  # NO
        290: 'GrayResponseUnit',           # NO
        291: 'GrayResponseCurve',          # NO
        292: 'Group3Options',              # NO
        293: 'Group4Options',              # NO
        296: 'ResolutionUnit',             # NO
        297: 'PageNumber',                 # NO, perhaps later for energy?
        301: 'ColorResponseCurves',        # NO
        305: 'Software',                   # NO, perhaps later
        306: 'DateTime',                   # NO, perhaps later
        315: 'Artist',                     # NO
        316: 'HostComputer',               # NO, perhaps later
        317: 'Predictor',                  # NO
        318: 'WhitePoint/ColorImageType',  # NO - depends on whether it is WORD(3)/RATIONAL(5)
        319: 'PrimaryChromaticities/ColorList', # NO - depends if it is RATIONAL(5)/BYTE(1)orWORD(3)
        320: 'ColorMap',                   # NO
        }

    def __init__(self):
        """Initialize object."""
        self.index = -1
        self.__field_type = None
        self.n_entries = 0
        self.__value = None

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
        except IndexError:
            raise ValueError(f"{self.__class__.__name__}: Unknown TiffTag "
                             f"field type ordinal {new_field_type}")

    @property
    def name(self):
        """Return the TiffTag name for self as a string."""
        return self.__tag_ids.get(self.index, str(self.index))

    @property
    def item_size(self):
        """Return the size in bytes of each entry in self."""
        return self.__bytes_per_entry[self.field_type]

    @property
    def size(self):
        """Return the total size of self in bytes."""
        return self.n_entries * self.item_size

    @property
    def value(self):
        """Return the value of this TiffTag."""
        if self.__value is None:
            raise RuntimeError(f"{self.__class__.__name__}: value was never "
                               "read. Call .value_from_bytes(bytes, order).")
        return self.__value

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
            value = struct.unpack(f"{byte_fmt}{self.n_entries}f",
                                  value_as_bytes)
        elif self.field_type == 'double':
            value = struct.unpack(f"{byte_fmt}{self.n_entries}d",
                                  value_as_bytes)
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
                value = tuple(int_from_bytes(v[:4], byte_order)
                              / int_from_bytes(v[4:], byte_order)
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


class TiffImage:
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
        self.__tiff_file = None

    @classmethod
    def from_array(cls, array):
        """Return a TiffImage object from an array.

        Parameters
        ----------
        array : numpy.ndarray
            Image data. Only unsigned 8-, 16- or 32-bit
            arrays can be used.

        Returns
        -------
        image : TiffImage
            The TiffImage object with pixel data from array

        Raises
        ------
        TypeError
            If array data type is invalid
        ValueError
            If array does not have exactly two dimensions
        """
        array = np.asarray(array)
        if not np.issubdtype(array.dtype, np.unsignedinteger):
            raise TypeError(f"{cls.__name__}: invalid data type {array.dtype} "
                            "for TIFF data. Expected unsigned integer.")
        if len(array.shape) < 2 or len(array.shape) > 3:
            raise ValueError(f"{cls.__name__}: array with {len(array.shape)} "
                             "dimensions is not a valid gray-scale image.")
        self = cls()
        if len(array.shape) == 3:
            self.__images = [img for img in array]
        else:
            self.__images = [array]
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
            tiff_id = int.from_bytes(tiff_file.read(2),
                                            self.__byte_order)
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

        if not self.__images:
            raise ValueError(f"{cls.__name__}: {filename.resolve()} "
                             "does not contain any images.")
        return self

    @property
    def images(self):
        """Return the current images as a list of numpy.ndarrays."""
        return self.__images

    def __read_ifd(self):
        """Read one image-file directory.

        Images found are stored in self.__images.

        Returns
        -------
        next_ifd_start : int
            Offset in bytes from the beginning of the file
            at which the next image-file directory starts.
            Will be zero if this was the last IFD.
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

        # Is there a better way to do this with numpy
        # from the bytes directly? Does it make sense
        # to have always dtype='>u2' further below??
        pixel_intensity = []
        bytes_per_px = int(tags['BitsPerSample']//8)
        for b_count, offset in zip(tags['StripByteCounts'],
                                   tags['StripOffsets']):
            self.__tiff_file.seek(offset)
            strip = self.__tiff_file.read(b_count)
            pixel_intensity.extend(
                int.from_bytes(strip[i:i+bytes_per_px], self.__byte_order)
                for i in range(0, b_count, bytes_per_px)
                )
        dtype = f'>u{bytes_per_px}'
        image = np.asarray(pixel_intensity, dtype=dtype).reshape(
            tags['ImageLength'],
            tags['ImageWidth']
            )

        if tags['PhotometricInterpretation'] == 0:
            image = 2**tags['BitsPerSample'] - 1 - image

        self.__images.append(image)

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











