# -*- coding: utf-8 -*-
"""
@author: Alexander M. Imre

Contains basic information about TensErLEED files and code to check their 
checksums before run time compilation. This is supposed to help avoid 
security vulnerabilities.

We use the common SHA-256 hashing algorithm as available from 
Pythons hashlib.
The check is toggled by parameter TL_IGNORE_CHECKSUM.
"""

import hashlib
from pathlib import Path
from viperleed.tleedmlib.files.checksums import INPUT_FILE_VERSIONS

# TensErLEED Versions
KNOWN_TL_VERSIONS = (
    '1.6', '1.61', '1.71', '1.72', '1.73',
)

KNOWN_TL_SECTIONS = (
    'refcalc', 'deltas', 'search', 'superpos'
)

# generate set of all files
TL_INPUT_FILES = set()
for version, input_files in INPUT_FILE_VERSIONS.items():
    (TLInputFile(file.key(), version, file.value()) for file in input_files)
    TL_INPUT_FILES.update()
    
# check version codes are valid
assert all(v in KNOWN_TL_VERSIONS for v in INPUT_FILE_VERSIONS.keys())

class TLInputFile:
    def __init__(self, name, versions, checksums):
        self.path = Path(name)
        self.file_name= self.path.name
        
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

def get_TL_version_files(version):
    version_files = (file for file in TL_INPUT_FILES if file.version == version)
    return version_files

def _generate_checksums_for_dir(path, pattern = "*.f*"):
    """Function for tleedm developers. Generates copy-paste-able
    string with all checksums for a directory.

    Args:
        path (pathlib path): Directory contaning the files.
        pattern (str, optional): File pattern, e.g. extension. Defaults to "*.f*".
    """
    
    for file in Path(path).glob(pattern):
        checksum = get_path_checksum(file)
        print(f"    '{str(file)[45:]}':\n        ('{checksum}', ),")