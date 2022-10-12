# -*- coding: utf-8 -*-
"""
@author: Alexander M. Imre

Contains basic information about TensErLEED files and code to check their 
checksums before run time compilation. This is supposed to help avoid 
security vulnerabilities.

We use the common SHA-256 hashing algorithm as available from Pythons hashlib.
The check is toggled by parameter TL_IGNORE_CHECKSUM.
"""

import hashlib
from pathlib import Path
from viperleed.tleedmlib.files.checksums import SOURCE_FILE_VERSIONS

# TensErLEED Versions
KNOWN_TL_VERSIONS = (
    '1.6', '1.61', '1.71', '1.72', '1.73',
)

# sections of TensErLEED - currently unused
KNOWN_TL_SECTIONS = (
    'ref-calc', 'r-factor', 'deltas', 'search', 'superpos', 'errors'
)

class TLSourceFile:
    """Class that holds information of TensErLEED source files.
    """
    def __init__(self, name, versions, checksums):
        self.path = Path(name)
        self.name = self.path.name
        
        assert version in KNOWN_TL_VERSIONS
        self._version = versions
        
        self._valid_checksums = checksums
    
    @property
    def valid_checksums(self):
        return self._valid_checksums
    
    @property
    def TL_version(self):
        return self._version
    
    @property
    def file_subdir(self):
        # subdir lib or src
        return self.path.parent.name
    
    @property
    def file_name(self):
        return self.path.name
    
    @property
    def file_name_w_subdir(self):
        file_name = self.file_name()
        subdir = self.file_subdir()
        return str(Path(subdir) / Path(file_name))

# generate set of all files
TL_INPUT_FILES = set()
for version, input_files in SOURCE_FILE_VERSIONS.items():
    for file, checksums in input_files.items():
        version_files = (TLSourceFile(file, version, checksums))
        TL_INPUT_FILES.add(version_files)
    
# check version codes are valid
assert all(v in KNOWN_TL_VERSIONS for v in SOURCE_FILE_VERSIONS.keys())

def get_TL_version_files(version):
    version_files = tuple(file for file in TL_INPUT_FILES if file.TL_version == version)
    return version_files

def _get_checksums(tl_version, filename):
    """Returns valid checksums for given TL version and filename as a tuple."""
    
    if tl_version not in KNOWN_TL_VERSIONS:
        raise ValueError(f"Unrecognized TensErLEED Version: {tl_version}")
    
    tl_version_files = get_TL_version_files(tl_version)
    
    applicable_files = tuple(file for file in tl_version_files if file.name == filename)
    if not applicable_files:
        raise ValueError(f"Unrecognized filename '{filename}' "
                         f"for TensErLEED version {tl_version}.")
        
    nested_valid_checksums = (file.valid_checksums for file in applicable_files)
    valid_checksums = tuple(cs for nest in nested_valid_checksums for cs in nest)
    
    return valid_checksums

def get_path_checksum(file_path):
    """Calculates and returns the SHA256 hash of the file at file_path.

    Args:
        file_path (_type_): _description_

    Returns:
        str: Checksum of file.
    """
    
    try:
        with file_path.open(mode='rb') as open_file:
            content = open_file.read()
    except FileNotFoundError:
        raise FileNotFoundError(f"Could not calculate checksum of file {file_path}. "
                                "File not found.")

    # get checksum
    file_checksum = hashlib.sha256(content).hexdigest()
    return file_checksum

def validate_checksum(tl_version, filename):
    """Compares checksum for filename with checksums stored 
    in files/checksums.py.

    Args:
        tl_version (str, float): TensErLEED version
        filename (str, path): filename of file to be checked
        
    Raises:
        ValueError: If tl_version or filename have an invalid type.
        RuntimeError: If checksum does not match any stored checksums for 
        that file.
    """
    # ensure TL version is valid
    if not (isinstance(tl_version, str) or isinstance(tl_version, float)):
        raise ValueError("Invalid type for tl_version")
    clean_tl_version = str(tl_version)
    
    if not tl_version in KNOWN_TL_VERSIONS:
        raise ValueError(f"Unrecognized TensErLEED version: {clean_tl_version}")
    
    
    # ensure filename is valid and cleaned up
    if not (isinstance(filename, str) or isinstance(filename, Path)):
        raise ValueError(f"Invalid type of filetype: {type(filename)}."
                         "Allowed are str and Path.")
    
    file_path = Path(filename)
    filename_clean = file_path.name # filename without path
    
    # get checksum
    file_checksum = get_path_checksum(file_path)
    
    # get tuple of reference checksums
    reference_checksums = _get_checksums(tl_version, filename_clean)
    
    if not file_checksum in reference_checksums:
        raise RuntimeError("SHA-256 checksum comparison failed for "
                           f"file {filename}.")
    return

def validate_multiple_files(files_to_check, logger, calc_part_name, version):
    """Validates multiple files by calling validate_checksum on each.

    Args:
        files_to_check (iterable of str of Path): Files to be validated.
        logger (logger): Logger from logging module to be used.
        calc_part_name (str): String to be written into log referring to 
            the part of the calculation (e.g. "refcalc").
        version (str): TensErLEED version used. To be taken from rp.TL_VERSION_STR.
    """
    for file_path in files_to_check:
        try:
            validate_checksum(version, file_path)
        except RuntimeError:
            logger.error("Error in checksum comparison of TensErLEED files for "
                        f"{calc_part_name}. Could not verify file {file_path}")
            raise
    # if you arrive here, checksums were successful
    logger.debug(f"Checksums of TensErLEED source files for {calc_part_name} validated.")
    return

def _generate_checksums_for_dir(path, patterns = ("*/GLOBAL", "*/*.f*")):
    """Function for tleedm developers. Generates copy-paste-able
    string with all checksums for a directory.

    Args:
        path (pathlib path): Directory contaning the files.
        pattern (str, optional): File pattern, e.g. extension. Defaults to "*.f*".
    """
    for pattern in patterns:
        for file in Path(path).glob(pattern):
            checksum = get_path_checksum(file)
            print(f"    '{str(file)[45:]}':\n        ('{checksum}', ),")