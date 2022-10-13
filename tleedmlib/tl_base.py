# -*- coding: utf-8 -*-
"""Module tl_base of viperleed.tleedmlib.

@author: Alexander M. Imre

Contains basic information about TensErLEED files and code to check
their checksums before run time compilation. This is supposed to help
avoid security vulnerabilities.

We use the common SHA-256 hashing algorithm as available from Pythons
hashlib. The check is toggled by parameter TL_IGNORE_CHECKSUM.
"""

import hashlib
from pathlib import Path

from viperleed.tleedmlib.files.checksums import SOURCE_FILE_VERSIONS


# TensErLEED Versions
KNOWN_TL_VERSIONS = (
    "1.6",
    "1.61",
    "1.71",
    "1.72",
    "1.73",  # TODO: use Version when available
)


# sections of TensErLEED - currently unused
KNOWN_TL_SECTIONS = ("ref-calc", "r-factor", "deltas", "search",
                     "superpos", "errors")


class InvalidChecksumError(Exception):
    """Exception for invalid checksums."""


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
        self._name = self.path.name

        assert version in KNOWN_TL_VERSIONS
        self._version = version

        self._valid_checksums = tuple(checksums)

    def __hash__(self):
        """Return a hash for this instance."""
        return hash((self._name, self.tl_version))

    @property
    def valid_checksums(self):
        """Return valid checksums as a tuple of str.

        Returns
        -------
        checksums : tuple of str
            The known checksums for this source file. The values are
            those given at instantiation, and are usually taken from
            tleedmlib.files.checksums.
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
        """Return the name of this source file as a string."""
        return self._name

    @property
    def file_name_w_subdir(self):
        """Return '<subdir>/<file_name>' as a string."""
        file_name = self.name
        subdir = self.file_subdir
        return str(Path(subdir) / file_name)


# check version codes are valid
assert all(v in KNOWN_TL_VERSIONS for v in SOURCE_FILE_VERSIONS)


# generate set of all files
TL_INPUT_FILES = set()
for f_version, input_files in SOURCE_FILE_VERSIONS.items():
    for file_, f_checksums in input_files.items():
        TL_INPUT_FILES.add(TLSourceFile(file_, f_version, f_checksums))


def get_tl_version_files(version):
    """Return a tuple of TLSourceFile instances for a given version."""
    version_files = (f for f in TL_INPUT_FILES if f.tl_version == version)
    return tuple(version_files)


def _get_checksums(tl_version, filename):
    """Return a tuple of valid checksums for given version and filename."""
    if tl_version not in KNOWN_TL_VERSIONS:
        raise ValueError(f"Unrecognized TensErLEED Version: {tl_version}")

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
            f"Could not calculate checksum of file {file_path}. "
            "File not found."
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

    The known checksums are stored in tleedmlib.files.checksums.

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
    ValueError
        If tl_version is not one of the known versions.
    InvalidChecksumError
        If checksum does not match any of the known checksums
        for the same file.
    """
    # ensure TL version is valid
    if not isinstance(tl_version, (str, float)):
        raise TypeError("Invalid type for tl_version")
    clean_tl_version = str(tl_version)

    if not tl_version in KNOWN_TL_VERSIONS:
        raise ValueError(f"Unrecognized TensErLEED version: {clean_tl_version}")

    # ensure filename is valid and cleaned up
    if not isinstance(filename, (str, Path)):
        raise TypeError(
            "Invalid type of filename: "
            f"{type(filename).__name__}. "
            "Allowed are str and Path."
        )

    file_path = Path(filename).resolve()
    filename_clean = file_path.name  # filename without path

    # get checksum
    file_checksum = get_file_checksum(file_path)

    # get tuple of reference checksums
    reference_checksums = _get_checksums(tl_version, filename_clean)

    if file_checksum not in reference_checksums:
        raise InvalidChecksumError(
            f"SHA-256 checksum comparison failed for file {filename}."
        )


def validate_multiple_files(files_to_check, logger, calc_part_name, version):
    """Validate multiple files by calling validate_checksum on each.

    Parameters
    ----------
    files_to_check : iterable of str of Path
        Files to validate. Notice that the iterable will be consumed.
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
        try:
            validate_checksum(version, file_path)
        except InvalidChecksumError:
            logger.error(
                "Error in checksum comparison of TensErLEED files for "
                f"{calc_part_name}. Could not verify file {file_path}"
                )
            problematic.append(file_path)

    if problematic:
        raise InvalidChecksumError(
            "SHA-256 checksum comparison failed for files "
            f"{', '.join(problematic)}."
        )

    # if you arrive here, checksums were successful
    logger.debug(
        f"Checksums of TensErLEED source files for {calc_part_name} validated."
    )


def _generate_checksums_for_dir(path, patterns=("*/GLOBAL", "*/*.f*")):
    """Generate copy-paste-able string with all checksums for a directory.

    This function is intended for tleedm developers. The checksums are
    printed to standard out.

    Parameters
    ----------
    path : str, or path-like
        Directory contaning the files.
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
    checksum_dict = {}
    for pattern in patterns:
        for file in Path(path).glob(pattern):
            checksum = get_file_checksum(file)
            print(f"    '{str(file)[45:]}':\n        ('{checksum}', ),")
            checksum_dict[file] = (checksum,)
    return checksum_dict
