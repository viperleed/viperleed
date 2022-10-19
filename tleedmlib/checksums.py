# -*- coding: utf-8 -*-
"""Module checksums of viperleed.tleedmlib.

@author: Alexander M. Imre

Contains basic information about TensErLEED files and code to check
their checksums before run time compilation. This is supposed to help
avoid security vulnerabilities.

This module can be executed with argument -p "path/to/tensorleed/folder"
to generate the file _checksums.dat which contains an encoded version
of checksums for all TensErLEED Fortran source code files.

We use the common SHA-256 hashing algorithm as available from Pythons
hashlib. The check is toggled by parameter TL_IGNORE_CHECKSUM.
"""

import hashlib
import base64
import ast
from pathlib import Path
import argparse

# TensErLEED Versions
KNOWN_TL_VERSIONS = (
    "1.6",
    "1.61",
    "1.71",
    "1.72",
    "1.73",  # TODO: use Version when available
    "1.8",
)


# sections of TensErLEED - currently unused
KNOWN_TL_SECTIONS = ("ref-calc", "r-factor", "deltas", "search",
                     "superpos", "errors")


class UnknownTensErLEEDVersionError(Exception):
    """Exception for invalid TensErLEED version."""


class InvalidChecksumError(Exception):
    """Exception for invalid checksums."""


class InvalidChecksumFileError(Exception):
    """Exception for invalid checksum file."""


class TLSourceFile:
    """Class that holds information of TensErLEED source files."""

    def __init__(self, name, version, checksums):
        """Initialize TensErLEED source file instance.

        Parameters
        ----------
        name : str, or path-like
            Path to the source file
        version : str
            String with TensErLEED version. Must be part of KNOWN_TL_VERSIONS.
        checksums : Sequence of str
            Valid checksums for this TensErLEED version and file.
            Multiple checksums may be permissible per version to allow
            for minor patches.

        Returns
        -------
        None.
        """
        self.path = Path(name).resolve()
        # Get tensorleed folder:
        #            xxx.f  /  src  / TensErLEED- / tensorleed
        base_path = self.path.parent.parent.parent
        self._name = self.path.relative_to(base_path)

        if version not in KNOWN_TL_VERSIONS:
            raise UnknownTensErLEEDVersionError(f"Unknown TensErLEED version {version}")
        self._version = version

        self._valid_checksums = tuple(checksums)

    def __hash__(self):
        """Return a hash for this instance."""
        return hash((self._name, self.tl_version))

    def __repr__(self):
        """Return a string representation of self."""
        txt = f"TLSourceFile({self.name}, {self.tl_version}, "
        return txt + f"{self.valid_checksums})"

    @property
    def valid_checksums(self):
        """Return valid checksums as a tuple of str.

        Returns
        -------
        checksums : tuple of str
            The known checksums for this source file. The values are
            those given at instantiation, and are usually taken from
            _checksums.dat
        """
        return self._valid_checksums

    @property
    def tl_version(self):
        """Return the version of this source as a string."""
        return self._version

    @property
    def file_subdir(self):
        """Return the subdirectory of this file as a string.

        Returns
        -------
        subdir : str
            Currently this is 'src' or 'lib', depending on
            whether the file is a "PROGRAM" ('src') or if
            it contains a library of subroutines.
        """
        return self.path.parent.name

    @property
    def name(self):
        """Return the qualified name of this source file as a Path."""
        return self._name


def get_tl_version_files(version):
    """Return a tuple of TLSourceFile instances for a given version."""
    version_files = (f for f in TL_INPUT_FILES if f.tl_version == version)
    return tuple(version_files)


def _get_checksums(tl_version, filename):
    """Return a tuple of valid checksums for given version and filename."""
    if tl_version not in KNOWN_TL_VERSIONS:
        raise UnknownTensErLEEDVersionError(
            f"Unrecognized TensErLEED Version: {tl_version}"
        )

    tl_version_files = get_tl_version_files(tl_version)

    applicable_files = tuple(f for f in tl_version_files if f.name == filename)
    if not applicable_files:
        raise ValueError(
            f"Unrecognized filename '{filename}' "
            f"for TensErLEED version {tl_version}."
        )

    nested_valid_checksums = (f.valid_checksums for f in applicable_files)
    valid_checksums = (cs for nest in nested_valid_checksums for cs in nest)
    return tuple(valid_checksums)


def get_file_checksum(file_path):
    """Return the SHA256 hash of the file at file_path.

    Parameters
    ----------
    file_path : str, or path-like
        Path to the file whose checksum should be calculated.

    Returns
    -------
    checksum : str
        Checksum of file.

    Raises
    ------
    FileNotFoundError
        If file_path does not exist or if it is not a file.
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(
            f"Could not calculate checksum of file {file_path}. File not found."
        )
    if not file_path.is_file():
        raise FileNotFoundError(
            f"Could not calculate checksum of file {file_path}. Not a file."
        )
    with file_path.open(mode="rb") as open_file:
        content = open_file.read()
    return hashlib.sha256(content).hexdigest()


def validate_checksum(tl_version, filename):
    """Compare checksum for filename with known checksums.

    The known checksums are stored (encoded) in _checksums.dat.

    Parameters
    ----------
    tl_version : str, or float
        TensErLEED version
    filename : str, or path-like
        Path to the file to be checked

    Raises
    ------
    TypeError
        If tl_version or filename have an invalid type.
    UnknownTensErLEEDVersionError
        If tl_version is not one of the known versions.
    InvalidChecksumFileError
        If no known checksum exists for filename. Consider
        running 'python checksums -p "path/to/tensorleed"' to
        generate a new checksum file.
    InvalidChecksumError
        If checksum does not match any of the known checksums
        for the same file.
    """
    # ensure TL version is valid
    if not isinstance(tl_version, (str, float)):
        raise TypeError("Invalid type for tl_version")
    clean_tl_version = str(tl_version)

    if clean_tl_version not in KNOWN_TL_VERSIONS:
        raise UnknownTensErLEEDVersionError(
            f"Unrecognized TensErLEED version: {clean_tl_version}"
        )

    # ensure filename is valid and cleaned up
    if not isinstance(filename, (str, Path)):
        raise TypeError(
            "Invalid type of filename: "
            f"{type(filename).__name__}. "
            "Allowed are str and Path."
        )

    file_path = Path(filename).resolve()
    base_path = file_path.parent.parent.parent
    filename_clean = file_path.relative_to(base_path)

    # get checksum
    file_checksum = get_file_checksum(file_path)

    # get tuple of reference checksums
    try:
        reference_checksums = _get_checksums(tl_version, filename_clean)
    except ValueError as err:
        raise InvalidChecksumFileError(
            f"Could not find checksum for file {filename_clean}."
        ) from err

    if file_checksum not in reference_checksums:
        raise InvalidChecksumError(
            f"SHA-256 checksum comparison failed for file {filename}."
        )


def validate_multiple_files(files_to_check, logger, calc_part_name, version):
    """Validate multiple files by calling validate_checksum on each.

    Parameters
    ----------
    files_to_check : iterable of str, pathlike, None
        Files to validate. Notice that the iterable will be consumed.
        If element is None, it will be skipped. This way we can deal 
        with optional files.
    logger : logging.Logger
        Logger from logging module to be used.
    calc_part_name : str
        String to be written into log referring to
        the part of the calculation (e.g. "refcalc").
    version : str
        TensErLEED version used. To be taken from rp.TL_VERSION_STR.

    Raises
    ------
    InvalidChecksumError
        If any of the files fails the checksum check.
    """
    problematic = []
    for file_path in files_to_check:
        if file_path is None: # may be passed a None for e.g. muftin -> skip
            continue
        try:
            validate_checksum(version, file_path)
        except (InvalidChecksumError, InvalidChecksumFileError) as err:
            logger.error(
                "Error in checksum comparison of TensErLEED files for "
                f"{calc_part_name}. Could not verify file {file_path}."
                f" Info: {err}"
            )
            problematic.append(str(file_path))

    if problematic:
        txt = ', '.join(problematic)
        raise InvalidChecksumError(
            f"SHA-256 checksum comparison failed for files {txt}."
        )

    # if you arrive here, checksums were successful
    logger.debug(
        f"Checksums of TensErLEED source files for {calc_part_name} validated."
    )


def encode_checksums(source_file_versions):
    """Return a base-64 encoded version of source_file_versions.

    Parameters
    ----------
    source_file_versions : dict
        Keys are TensErLEED versions as strings, values are dicts
        {path:(checksum,)}.

    Returns
    -------
    encoded : bytes
        base-64 encoded version of source_file_versions.
    """
    bytes_source = repr(source_file_versions).encode("utf-8")
    return base64.b64encode(bytes_source)


def decode_checksums(encoded_checksums):
    """Return a dict of source_file_versions from its encoded form.

    This is the inverse of encode_checksums.

    Parameters
    ----------
    encoded_checksums : bytes
        Encoded version of source_file_versions.

    Returns
    -------
    source_file_versions : dict
        Keys are TensErLEED versions as strings, values are dicts
        {path:(checksum,)}.
    """
    encoded_checksums += b"==="  # Prevents padding errors
    bytes_source = base64.b64decode(encoded_checksums)
    try:
        return ast.literal_eval(bytes_source.decode("utf-8"))
    except (TypeError, ValueError,
            MemoryError, RecursionError,
            SyntaxError, UnicodeDecodeError) as err:
        raise InvalidChecksumFileError(
            "Looks like the TensErLEED source checksum file "
            "has been tampered with!"
        ) from err


def read_encoded_checksums(encoded_file_path=None):
    """Return checksums read from encoded_file_path.

    Parameters
    ----------
    encoded_file_path : str, or pathlike, optional
        Optional location of encoded checksum file. If None the default
        location (tleedmlib/_checksums.dat) is assumed. Default is
        None.

    Returns
    -------
    source_file_versions : dict
        Keys are TensErLEED versions as strings, values are dicts
        {path:(checksum,)}.
    """
    if encoded_file_path is None:
        # file should be in tleedmlib/
        encoded_file_path = Path(__file__).resolve().parent
        encoded_file_path /= "_checksums.dat"

    with open(encoded_file_path, "rb") as file:
        return decode_checksums(file.read())


def _write_encoded_checksums(source_file_versions, encoded_file_path=None):
    """Write source_file_versions to encoded_file_path.

    Parameters
    ----------
    source_file_versions : dict
        Keys are TensErLEED versions as strings, values are dicts
        {path:(checksum,)}.
    encoded_file_path : str, or pathlike, optional
        Optional location of encoded checksum file. If None the default
        location (tleedmlib/_checksums.dat) is assumed. Default is
        None.

    Returns
    -------
    None
    """
    if encoded_file_path is None:
        # file should be in tleedmlib/
        encoded_file_path = Path(__file__).resolve().parent
        encoded_file_path /= "_checksums.dat"

    with open(encoded_file_path, "wb") as file:
        file.write(encode_checksums(source_file_versions))


def _generate_checksums_for_dir(path, patterns=("*/GLOBAL", "*/*.f*")):
    """Generate copy-paste-able string with all checksums for a directory.

    This function is intended for tleedm developers. The checksums are
    printed to standard out.

    Parameters
    ----------
    path : str, or path-like
        Directory contaning the files. Typically one of the
        "TensErLEED-vXXX" directories in tensorleed.
    patterns : tuple of str, optional
        File patterns used to select files (e.g., extension).
        Syntax should be the one expected by glob. Default is
        ("*/GLOBAL", "*/*.f*").

    Returns
    -------
    checksum_dict : dict
        Keys are filenames, values are a tuple with the checksum for the
        file contents. The tuple will only contain one item.
    """
    checksum_dict_ = {}
    path = Path(path)
    base_path = path.parent  # /path/to/tensorleed
    for pattern in patterns:
        for file in path.glob(pattern):
            checksum = get_file_checksum(file)
            checksum_dict_[str(file.relative_to(base_path))] = (checksum,)
    return checksum_dict_


if __name__ != "__main__":
    # permissible checksums for various TensErLEED version
    # implemented as dicts for each version - key is filename
    # allowed checksums for each version are given in a tuple
    SOURCE_FILE_VERSIONS = read_encoded_checksums()

    # check version codes are valid
    try:
        assert all(v in KNOWN_TL_VERSIONS for v in SOURCE_FILE_VERSIONS)
    except AssertionError as err:
        raise UnknownTensErLEEDVersionError(
            "Unknown TensErLEED version detected in checksum file."
        ) from err

    # generate set of all files
    TL_INPUT_FILES = set()
    for f_version, input_files in SOURCE_FILE_VERSIONS.items():
        for file_, f_checksums in input_files.items():
            TL_INPUT_FILES.add(TLSourceFile(file_, f_version, f_checksums))

else:  # Write new checksum file when executed as a module
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p", "--tlpath",
        help=("specify TensErLEED source directory"), type=str
    )
    args = parser.parse_known_args()

    if args[0].tlpath:
        tl_base_path = Path(args[0].tlpath)
    else:
        raise RuntimeError("No path specified. Use --tlpath or -p")

    if not tl_base_path.exists():
        raise FileNotFoundError(f"Could not find {tl_base_path}")

    tl_folders = tuple(tl_base_path.glob("TensErLEED*"))
    if not tl_folders:
        raise FileNotFoundError(
            f"No TensErLEED folders found in {tl_base_path}"
        )

    checksum_dict = {}
    for folder in tl_folders:
        version_ = folder.name.split("v")[1]
        if version_ not in KNOWN_TL_VERSIONS:
            raise UnknownTensErLEEDVersionError(
                    f"Unknown TensErLEED version {version_} in {tl_base_path}"
            )
        checksums_ = _generate_checksums_for_dir(path=folder)
        checksum_dict[version_] = checksums_

    # write to file
    _write_encoded_checksums(checksum_dict)

    # try to read to make sure it's ok
    read_checksums = read_encoded_checksums()

    assert read_checksums == checksum_dict
    print("Wrote _checksums.dat successfully!")
