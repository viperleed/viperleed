"""Checksums for TensErLEED Fortran source code files.

Contains basic information about TensErLEED files and code to check
their checksums before run time compilation. This is supposed to help
avoid security vulnerabilities.

This module can be executed with argument -p "path/to/tensorleed/folder"
to generate the file _checksums.dat which contains an encoded version
of checksums for all TensErLEED Fortran source code files.

We use the common SHA-256 hashing algorithm as available from Pythons
hashlib. The check is toggled by parameter TL_IGNORE_CHECKSUM.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2022-10-12'
__license__ = 'GPLv3+'

import argparse
import ast
import base64
import hashlib
from pathlib import Path
from warnings import warn

# Where encoded checksums are stored
CHECKSUMS_FILE_NAME = '_checksums.dat'


# TensErLEED Versions
KNOWN_TL_VERSIONS = (
    '1.6',
    '1.61',
    '1.71',
    '1.72',
    '1.73',  # TODO: use Version when available
    '1.74',
    '1.75',
    '1.76',
    '1.8',
    )


# sections of TensErLEED - currently unused
KNOWN_TL_SECTIONS = ('ref-calc', 'r-factor', 'deltas', 'search',
                     'superpos', 'errors')


class UnknownTensErLEEDVersionError(ValueError):
    """Exception for invalid TensErLEED version."""

    def __init__(self, version=None, message=''):
        """Initialize exception."""
        self.version = version
        full_message = ('' if version is None
                        else f'Unrecognized TensErLEED version {version}.')
        if message and full_message:
            full_message += ' '
        if message:
            full_message += message
        super().__init__(full_message)


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
        name : str or pathlike
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
        # Get "tensorleed" parent folder:
        #     [2]       [1]         [0] path
        # tensorleed/TensErLEED-vXX/src/xxx.f
        base_path = self.path.parents[2]
        self._name = self.path.relative_to(base_path)

        if version not in KNOWN_TL_VERSIONS:
            raise UnknownTensErLEEDVersionError(version)
        self._version = version

        self._valid_checksums = set(checksums)

    def __hash__(self):
        """Return a hash for this instance."""
        return hash((self._name, self.tl_version))

    def __repr__(self):
        """Return a string representation of self."""
        txt = f'TLSourceFile({self.name}, {self.tl_version}, '
        return txt + f'{self.valid_checksums})'

    @property
    def valid_checksums(self):
        """Return valid checksums as a tuple of str.

        Returns
        -------
        checksums : set of str
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
        raise UnknownTensErLEEDVersionError(tl_version)

    tl_version_files = get_tl_version_files(tl_version)
    applicable_files = tuple(f for f in tl_version_files if f.name == filename)
    if not applicable_files:
        raise ValueError(f'Unrecognized filename {filename!r} '
                         f'for TensErLEED version {tl_version}.')

    nested_valid_checksums = (f.valid_checksums for f in applicable_files)
    valid_checksums = (cs for nest in nested_valid_checksums for cs in nest)
    return tuple(valid_checksums)


def get_file_checksum(file_path):
    """Return the SHA256 hash of the file at file_path.

    Parameters
    ----------
    file_path : str or pathlike
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
    file_path = Path(file_path).resolve()
    if not file_path.exists():
        raise FileNotFoundError('Could not calculate checksum of '
                                f'file {file_path}. File not found.')
    if not file_path.is_file():
        raise FileNotFoundError('Could not calculate checksum of '
                                f'file {file_path}. Not a file.')
    with file_path.open(mode='rb') as open_file:
        content = open_file.read()
    # Make sure we always have '\n' for line endings, and not '\r\n'.
    # The latter seems to appear in files synced with git on Windows:
    # line endings are automatically changed with the default
    # git/GitHub Desktop configuration.
    content = content.replace(b'\r\n', b'\n')
    return hashlib.sha256(content).hexdigest()


def validate_checksum(tl_version, filename):
    """Compare checksum for filename with known checksums.

    The known checksums are stored (encoded) in _checksums.dat.

    Parameters
    ----------
    tl_version : str or float
        TensErLEED version
    filename : str or pathlike
        Path to the file to be checked

    Raises
    ------
    TypeError
        If tl_version or filename have an invalid type.
    UnknownTensErLEEDVersionError
        If tl_version is not one of the known versions.
    InvalidChecksumFileError
        If no known checksum exists for filename. Consider
        running 'python -m checksums -p "path/to/tensorleed"'
        to generate a new checksum file.
    InvalidChecksumError
        If checksum does not match any of the known checksums
        for the same file.
    """
    # Ensure TL version is valid
    if not isinstance(tl_version, (str, float)):
        raise TypeError('Invalid type for tl_version')
    clean_tl_version = str(tl_version)

    if clean_tl_version not in KNOWN_TL_VERSIONS:
        raise UnknownTensErLEEDVersionError(clean_tl_version)

    # Ensure filename is valid and cleaned up
    if not isinstance(filename, (str, Path)):
        raise TypeError('Invalid type of filename: '
                        f'{type(filename).__name__}. '
                        'Allowed are str and Path.')

    file_path = Path(filename).resolve()
    base_path = file_path.parents[2]  # 3 folders up
    filename_clean = file_path.relative_to(base_path)

    # Get checksum
    file_checksum = get_file_checksum(file_path)

    # Get known checksums
    try:
        reference_checksums = _get_checksums(clean_tl_version, filename_clean)
    except ValueError as exc:
        raise InvalidChecksumFileError('Could not find checksum '
                                       f'for file {filename_clean}.') from exc

    if file_checksum not in reference_checksums:
        raise InvalidChecksumError('SHA-256 checksum comparison '
                                   f'failed for file {filename}.')


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
        if file_path is None: # May be None, e.g. for muftin.f -> skip
            continue
        try:
            validate_checksum(version, file_path)
        except (InvalidChecksumError, InvalidChecksumFileError) as exc:
            logger.error(
                'Error in checksum comparison of TensErLEED files for '
                f'{calc_part_name}. Could not verify file {file_path}.'
                f' Info: {exc}'
                )
            problematic.append(str(file_path))
    if problematic:
        txt = ', '.join(problematic)
        raise InvalidChecksumError(
            f'SHA-256 checksum comparison failed for files {txt}.'
            )
    # If you arrive here, checksums were successful
    logger.debug('Checksums of TensErLEED source '
                 f'files for {calc_part_name} validated.')


def encode_checksums(source_file_checksums):
    """Return a base-64 encoded version of source_file_checksums.

    Parameters
    ----------
    source_file_checksums : dict
        Keys are paths to files as strings, values are sets of known
        checksums. Paths are relative to the tensorleed folder.

    Returns
    -------
    encoded : bytes
        base-64 encoded version of source_file_checksums.
    """
    bytes_source = repr(source_file_checksums).encode('utf-8')
    return base64.b64encode(bytes_source)


def decode_checksums(encoded_checksums):
    """Return a dict of source_file_checksums from its encoded form.

    This is the inverse of encode_checksums.

    Parameters
    ----------
    encoded_checksums : bytes
        Encoded version of source_file_checksums.

    Returns
    -------
    source_file_checksums : dict
        Keys are paths to files as strings, values are sets of known
        checksums. Paths are relative to the tensorleed folder.

    Raises
    ------
    InvalidChecksumFileError
        When decoding of encoded_checksums fails.
    """
    encoded_checksums += b'==='  # Prevents padding errors
    bytes_source = base64.b64decode(encoded_checksums)
    try:
        return ast.literal_eval(bytes_source.decode('utf-8'))
    except (TypeError, ValueError, MemoryError,
            RecursionError, SyntaxError) as exc:
        raise InvalidChecksumFileError(
            'Looks like the TensErLEED source '
            'checksum file has been tampered with!'
            ) from exc


def read_encoded_checksums(encoded_file_path=None):
    """Return checksums read from encoded_file_path.

    Parameters
    ----------
    encoded_file_path : str or pathlike, optional
        Optional location of encoded checksum file. If None, the
        default location (viperleed/calc/_checksums.dat) is assumed.
        Default is None.

    Returns
    -------
    source_file_checksums : dict
        Keys are paths to files as strings, values are sets of known
        checksums. Paths are relative to the tensorleed folder.
    """
    if encoded_file_path is None:
        # file should be in calc/
        encoded_file_path = Path(__file__).resolve().parent
        encoded_file_path /= CHECKSUMS_FILE_NAME

    with open(encoded_file_path, 'rb') as file:
        return decode_checksums(file.read())


def _write_encoded_checksums(source_file_checksums, encoded_file_path=None):
    """Write source_file_checksums to encoded_file_path.

    Parameters
    ----------
    source_file_checksums : dict
        Keys are paths to files as strings, values are sets of known
        checksums. Paths are relative to the tensorleed folder.
    encoded_file_path : str, or pathlike, optional
        Optional location of encoded checksum file. If None, the
        default location (viperleed/calc/_checksums.dat) is assumed.
        Default is None.

    Returns
    -------
    None.
    """
    if encoded_file_path is None:
        # file should be in calc/
        encoded_file_path = Path(__file__).resolve().parent
        encoded_file_path /= CHECKSUMS_FILE_NAME

    with open(encoded_file_path, 'wb') as file:
        file.write(encode_checksums(source_file_checksums))


def _add_checksums_for_dir(path,
                           checksum_dict_,
                           patterns=('*/GLOBAL', '*/*.f*')):
    """Add checksums for files in path into checksum_dict_.

    This function is intended for viperleed developers.

    Parameters
    ----------
    path : str, or path-like
        Directory containing the files. Typically one of the
        "TensErLEED-vXXX" directories in tensorleed.
    checksum_dict_ : dict
        Checksums dict to start with. Should be {path: {checksums}}.
        New checksums will be added to existing values.
        Modified in place.
    patterns : tuple of str, optional
        File patterns used to select files (e.g., extension).
        Syntax should be the one expected by glob. Default is
        ("*/GLOBAL", "*/*.f*").

    Returns
    -------
    None.
    """
    path = Path(path).resolve()
    base_path = path.parent  # /path/to/tensorleed
    for pattern in patterns:
        for file in path.glob(pattern):
            checksum = get_file_checksum(file)
            key = str(file.relative_to(base_path).as_posix())
            if key not in checksum_dict_:
                checksum_dict_[key] = set()
            checksum_dict_[key].add(checksum)


def _parse_args(args):
    """Parse command line arguments.

    Parameters
    ----------
    args : argparse.Namespace
        First element returned by ArgumentParser.parse_known_args()

    Returns
    -------
    tl_base_path : Path
        Base path of folder containing TensErLEED source files.
    tl_folders : tuple
        Contains Paths to all subfolders contained in tl_base_path
    checksum_dict : dict
        Dictionary containing {filename:set(known_checksums)}.

    Raises
    ------
    RuntimeError
        If no valid tensor-LEED path was specified as a command-line
        argument.
    FileNotFoundError
        If the path given in args does not exist or does not contain
        source files.
    """
    tl_base_path = _resolve_tensorleed_path_argument(args)
    tl_folders = tuple(tl_base_path.glob('TensErLEED*'))                        # TODO: would need to be changed to include beamgen, etc.
    if not tl_folders:
        raise FileNotFoundError(
            f'No TensErLEED folders found in {tl_base_path}'
            )

    # Check for --no-append flag
    checksum_dict = {}
    if not args.no_append:
        try:
            checksum_dict = read_encoded_checksums()
        except (FileNotFoundError, InvalidChecksumFileError) as exc:
            warn(f'Could not read {CHECKSUMS_FILE_NAME} file. '
                 f'Creating a new one. Info: {exc}')

    return tl_base_path, tl_folders, checksum_dict


def _resolve_tensorleed_path_argument(args):
    """Return a path to the tensor-LEED folder from CLI args."""
    if args.tlpath:
        tl_base_path = Path(args.tlpath).resolve()
    else:
        raise RuntimeError("No path specified. Use --tlpath or -p")

    if not tl_base_path.exists():
        raise FileNotFoundError(f"Could not find {tl_base_path}")
    return tl_base_path


def _add_parser_args(parser):
    """Add CLI arguments to parser."""
    parser.add_argument('-p', '--tlpath',
                        help='Specify TensErLEED source directory',
                        type=str)
    parser.add_argument(
        '-n', '--no-append',
        help=(f'Do not read in existing {CHECKSUMS_FILE_NAME} file and '
              'create a new one instead containing ONLY the current files. '
              'If not specified, the default is to append new checksums '
              'to existing ones'),
        action='store_true'
        )


if __name__ != "__main__":
    # Permissible checksums for various source files
    # in the form {file_path: set(known_checksums)}
    VALID_CHECKSUMS = read_encoded_checksums()

    # Generate set of all files
    TL_INPUT_FILES = set()
    for file_, f_checksums in VALID_CHECKSUMS.items():
        folder = tuple(Path(file_).parents)[-2]
        if not folder.name.startswith("TensErLEED"):
            raise NotImplementedError(
                "Checksums are only implemented "
                "for TensErLEED files at the moment."
                )
        f_version = folder.name.split("v")[1]
        TL_INPUT_FILES.add(TLSourceFile(file_, f_version, f_checksums))

else:  # Write new checksum file when executed as a module
    parser = argparse.ArgumentParser()
    _add_parser_args(parser)
    _args = parser.parse_known_args()

    _tl_base_path, _tl_folders, _checksum_dict = _parse_args(_args[0])

    for folder in _tl_folders:
        version_ = folder.name.split("v")[1]
        if version_ not in KNOWN_TL_VERSIONS:
            raise UnknownTensErLEEDVersionError(
                    f"Unknown TensErLEED version {version_} in {_tl_base_path}"
                    )
        _add_checksums_for_dir(folder, _checksum_dict)

    # write to file
    _write_encoded_checksums(_checksum_dict)

    # try to read to make sure it's ok
    read_checksums = read_encoded_checksums()

    assert read_checksums == _checksum_dict
    print("Wrote _checksums.dat successfully!")
