"""Module manifest of viperleed.calc.files.

Defines the ManifestFile class for handling files and directories
created by viperleed.calc during execution and that should be
copied back to the folder in which viperleed.calc was originally
executed.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-02-24'
__license__ = 'GPLv3+'

from pathlib import Path
import re

_HEADER_LINE_RE = re.compile(r'\[(?P<label>.*)at (?P<path>.*)\]')
_MANIFEST_NAME = 'manifest'


class ManifestFileError(Exception):
    """Base exception for manifest-related errors."""


class InconsistentPathsError(ManifestFileError):
    """The path of a content does not match the one of its header."""


class ManifestFile:
    """A helper to handle calc-generated files/folders.

    A manifest file has the following structure:
        root_file_or_folder_1
        ...

        [folder_1_label at folder_1_path]
        folder_1_path/file_or_folder_1
        ...

        [folder_2_label at folder_2_path]
        ...

    If folder_X is a subfolder of the main file's path, folder_X_path
    is a relative path. Notice that header lines for (sub)folders are
    always present, irrespective of whether they have any contents.
    """

    def __init__(self, *contents, root_path=None):
        """Initialize an instance with some contents."""
        self._root = Path(root_path or '').resolve()
        self._sections = {self._root: set()}  # {path: {str_contents}}
        self._labels = {self._root: 'root'}
        for content in contents:
            if isinstance(content, ManifestFile):
                self += content
            else:
                self.add(content)

    @property
    def has_absolute_paths(self):
        """Return whether any of the sub-manifests is not a subfolder."""
        for path in self.sections[1:]:
            try:
                path.relative_to(self.root)
            except ValueError:
                return True
        return False

    @property
    def name(self):
        """Return the default name of the manifest file."""
        return _MANIFEST_NAME

    @property
    def path(self):
        """Return the path to this manifest file."""
        return self.root / self.name

    file = path

    @property
    def root(self):
        """Return the path to the folder containing this ManifestFile."""
        return self._root

    @property
    def sections(self):
        """Return a tuple of all the paths in this ManifestFile."""
        return tuple(self._sections)

    def __contains__(self, item):
        """Return whether `item` is contained in this ManifestFile."""
        return item in self._sections[self.root]

    def __iadd__(self, other):
        """Add the contents of another manifest file to this one."""
        if not isinstance(other, ManifestFile):
            return NotImplemented
        for path, contents in other.iter_sections():
            if path in self._sections:
                self._sections[path].update(contents)
            else:
                self._sections[path] = contents  # Same object!
        return self

    def add(self, item):
        """Add an `item` to the top-level path."""
        self._sections[self.root].add(item)

    def add_manifest(self, other, label='folder'):
        """Add another manifest file with an optional label."""
        self += other
        if other.root != self.root:  # Don't override root
            self._labels[other.root] = label

    def iter_sections(self, relative=False):
        """Yield sections and their contents.

        Parameters
        ----------
        relative : bool, optional
            Whether paths to the sections should be relative to
            self.root.

        Yields
        ------
        path : Path
            Path to each (sub)section of this ManifestFile.
        contents : set of str
            Names of files/folders at path.
        """
        contents = self._sections.items()
        if relative and self.has_absolute_paths:
            raise ManifestFileError('relative paths are supported only if all '
                                    'paths are subfolders of the root one.')
        if relative:
            contents = ((p.relative_to(self.root), c) for p, c in contents)
        yield from contents

    def read(self):
        """Read the contents of a manifest file in the current directory."""
        if not self.file.is_file():
            return
        section = self.root
        rel_path = None
        with self.file.open('r', encoding='utf-8') as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue
                try:
                    section, rel_path = self._read_header_line(line)
                except ValueError:  # It's some contents
                    self._read_contents_line(line, section, rel_path)

    def write(self):
        """Write the contents of this ManifestFile to file."""
        with self.file.open('w', encoding='utf-8') as file:
            root_contents = self._sections[self.root]
            if root_contents:
                self._write_contents('', root_contents, file)
            for path, contents in self.iter_sections():
                if path is self.root:
                    continue
                try:
                    relative_path = path.relative_to(self.root)
                except ValueError:  # Not a subfolder, write it in full
                    relative_path = path
                label = self._labels.get(path, '')
                if label:
                    label += ' '
                file.write(f'[{label}at {relative_path.as_posix()}]\n')
                self._write_contents(relative_path, contents, file)

    def _read_contents_line(self, line, section, rel_path):
        """Store the contents of `line` into a `section`."""
        if not rel_path:  # Skip check for root files
            content = line
        elif not line.startswith(rel_path):
            raise InconsistentPathsError(f'{line} should begin with '
                                         f'{rel_path}.') from None
        else:
            *_, content = line.split(rel_path + '/', maxsplit=1)
        self._sections[section].add(str(content))

    def _read_header_line(self, line):
        """Return the path of a subsection and its relative version as str."""
        line = line.strip()
        match_ = _HEADER_LINE_RE.fullmatch(line)
        if not match_:
            raise ValueError(f'Not a header line: {line!r}.')
        path_str = match_['path'].strip()
        path = Path(path_str).resolve()
        self._labels[path] = match_['label'].strip()
        self._sections[path] = set()
        return path, str(Path(path_str).as_posix())

    @staticmethod
    def _write_contents(path, contents, file):
        """Write the `contents` of a section to `file`."""
        if path:
            contents = ((path/item).as_posix() for item in contents)
        lines = (f'{item}\n' for item in contents)
        file.writelines(lines)
        file.write('\n')
