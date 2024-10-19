"""Module meta of viperleed.calc.bookkeeper.history.

Defines classes to handle the metadata files stored to history together
with files from the root folder. Metadata files contain information on
what was archived to history at which point in time, especially which
history folders were created simultaneously.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-10-13'
__license__ = 'GPLv3+'

import hashlib
from configparser import ConfigParser
from pathlib import Path

from viperleed.calc.lib.checksums import with_posix_line_endings
from viperleed.calc.lib.string_utils import strip_comments


_CHUNK_SIZE = 4096  # bytes to read from file
_HEADER = '''
# This file was automatically generated by the bookkeeper utility of
# viperleed.calc. It contains information that should NOT BE EDITED.
# Editing this file may impair the ability of the bookkeeper of
# interpreting or reconstructing the history of this calculation.
'''
_METADATA_NAME = 'bookkeeper_meta'
_EMPTY = '_not_defined_'
_SECTIONS = {
    'archived' : {
        # Contains synthetic information about what was archived in a
        # history directory and which other directories were archived
        # with it.
        'hash': _EMPTY,
        'with': _EMPTY,  # hash of 'main' directory
        },
    }


# Hashing inspired by https://stackoverflow.com/questions/24937495, but
# uses posix-style paths relative to the root one instead of bare file
# or directory names to hash.


class BookkeeperMetaFile:
    """A handler for bookkeeper_meta files."""

    def __init__(self, path):
        """Initialize a metadata file for `path`."""
        self._hash = None
        self._path = Path(path)
        self._parser = ConfigParser()
        self._parser.read_dict(_SECTIONS)

    @property
    def file(self):
        """Return the path to the metadata file handled."""
        return self.path / _METADATA_NAME

    @property
    def hash_(self):
        """Return the hash of the history folder."""
        if self._hash is None:
            raise AttributeError('No hash_ available. '
                                 'Call compute_hash first.')
        try:
            return self._hash.hexdigest()
        except AttributeError:
            assert isinstance(self._hash, str)
            return self._hash

    @property
    def parent(self):
        """Return the hash of the folder created together with this one."""
        parent = self._parser['archived']['with']
        return None if parent == _EMPTY else parent

    @property
    def path(self):
        """Return the path where this metadata file is stored."""
        return self._path

    def collect_from(self, other):
        """Store information read from another metadata file in here.

        Parameters
        ----------
        other : BookkeeperMetaFile
            The metadata file of another history subfolder from which
            the following information is collected:
            - hash:
                This BookkeeperMetaFile will be considered as created
                together with `other`.
        """
        self._parser['archived']['with'] = other.hash_

    def compute_hash(self):
        """Compute and store a hash from all the files in path."""
        if not any(self._hashable_contents(self.path)):
            return
        self._hash = hashlib.md5(self.path.name.encode())
        self._update_hash_from_folder(self.path)
        self._parser['archived']['hash'] = self.hash_

    def read(self):
        """Read the metadata file in the root folder."""
        if not self.file.is_file():
            raise FileNotFoundError(f'No {self.file.name} file at {self.path}')
        with self.file.open('r', encoding='utf-8') as meta_file:
            # Skip header lines
            lines = (line for line in meta_file if strip_comments(line))
            self._parser.read_file(lines)
        self._hash = self._parser['archived']['hash']

    def write(self):
        """Store the information to file."""
        if self._hash is None:
            raise TypeError('Cannot write to file without a hash. '
                            'Call compute_hash first.')
        with self.file.open('w', encoding='utf-8') as meta_file:
            meta_file.write(_HEADER)
            self._parser.write(meta_file)

    def _hashable_contents(self, folder_path):
        """Collect all contents of `folder_path` used for hashing."""
        assert folder_path.is_dir()
        contents = (p for p in folder_path.iterdir()
                    # Keep only files and directories
                    if (p.is_file() or p.is_dir())
                    # Skip the metadata file itself
                    and p.name != _METADATA_NAME)
        return contents

    def _update_hash_from_file(self, file_path):
        """Update self._hash with the contents of a file."""
        assert file_path.is_file()
        with file_path.open('rb') as file:
            for chunk in _read_chunked(file):
                self._hash.update(with_posix_line_endings(chunk))

    def _update_hash_from_folder(self, folder_path):
        """Recursively update self._hash with the contents of a folder."""
        sorted_contents = sorted(self._hashable_contents(folder_path),
                                 key=lambda p: p.name.lower())
        for path in sorted_contents:
            relative_path = path.relative_to(self.path).as_posix()
            self._hash.update(str(relative_path).encode())
            if path.is_file():
                self._update_hash_from_file(path)
                continue
            # A folder
            self._update_hash_from_folder(path)



def _read_chunked(file, chunk_size=_CHUNK_SIZE):
    """Yield bytes read from `file` in `chunk_size`-long bits."""
    while True:
        chunk = file.read(chunk_size)
        if not chunk:
            return
        yield chunk
