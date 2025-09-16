"""Module add_dunders for viperleed files.

Adds the __copyright__ and __license__ module dunders for given files.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-03-22'
__license__ = 'GPLv3+'

import ast
from collections import Counter
from pathlib import Path
import re


DUNDER_RE = re.compile(r'^__[a-z]+__$')
DUNDER_EQUALS_RE = re.compile(r'^(?P<dunder>__[a-z]+__) = ')
ACCEPTABLE_DUNDER_TYPES = (
    ast.Tuple,
    ast.Constant,
    )
IMPORT_RE = re.compile(r'^\s*(from\s+\w+\s+)?import\s+\w+')


class DunderError(Exception):
    """Base exception for dunder-related errors."""


class MissingDunderError(DunderError):
    """Did not find a module dunder."""

    def __init__(self, dunder=None):
        """Initialize exception."""
        self.dunder = dunder
        message = 'Could not find '
        if dunder:
            message += dunder
        else:
            message += 'any module dunder'
        super().__init__(message)


class MisplacedDunderError(DunderError):
    """A dunder is found after imports."""

    def __init__(self, *dunders):
        """Initialize exception."""
        super().__init__('Found dunders after imports: '
                         + ', '.join(dunders))


class DuplicateDundersError(DunderError):
    """Dunders were found multiple times."""

    def __init__(self, dunders_and_counts):
        """Initialize exception."""
        msg = ('Found duplicate dunders: '
               ', '.join(f'{d} ({c})' for d, c in dunders_and_counts))
        super().__init__(msg)


class DunderSortingError(DunderError):
    """Dunders are not alphabetically sorted."""

    def __init__(self):
        """Initialize exception."""
        super().__init__('Found dunders sorted not alphabetically')


class DunderFormatError(DunderError):
    """Something is wrong with the format of a dunder."""


def ast_parse_dunder_lines(contents):
    """Return a dict of {dunder_name: ast.Node} from contents."""
    dunder_start_lines = find_dunders(contents)
    first_dunder = next(iter(dunder_start_lines.values()))
    first_import = find_first_import_line(contents)
    dunder_parsed = ast.parse(''.join(contents[first_dunder:first_import]))
    dunder_nodes = [dunder_node for dunder_node in dunder_parsed.body
                    if isinstance(dunder_node, ast.Assign)]
    if any(len(dunder_node.targets) > 1 for dunder_node in dunder_nodes):
        raise DunderFormatError('More than one dunder in a line found')
    return {dunder_node.targets[0].id: dunder_node.value
            for dunder_node in dunder_nodes
            if DUNDER_RE.match(dunder_node.targets[0].id)}


def find_dunders(contents):
    """Return a dict of {dunder: line} for all the module dunders found."""
    return {line: i for i, line in enumerate(contents)
            if DUNDER_EQUALS_RE.match(line)}


def find_dunder_line_ranges(contents):
    """Return {dunder: (line_start, line_end)} for dunders in contents."""
    dunders_start_lines = find_dunders(contents)
    dunder_items = tuple(dunders_start_lines.items())
    dunder_line_numbers = {}
    for i, (dunder_line, start_line) in enumerate(dunder_items):
        dunder = DUNDER_EQUALS_RE.match(dunder_line)['dunder']
        try:
            _, next_line = dunder_items[i+1]
        except IndexError:  # Last dunder. Use first empty line
            try:
                next_line = start_line + next(
                    i for i, line in enumerate(contents[start_line:])
                    if not line.strip()
                    )
            except StopIteration:
                next_line = len(contents)
        dunder_line_numbers[dunder] = (start_line, next_line - 1)
    return dunder_line_numbers


def find_imports(contents):
    """Return line indices that contain imports (except __future__)."""
    return [i for i, line in enumerate(contents)
            if IMPORT_RE.match(line) and 'from __future__' not in line]


def find_first_import_line(contents):
    """Return the first import line (except __future__)."""
    try:
        return min(find_imports(contents))
    except ValueError:  # No imports
        pass
    return len(contents)


def check_dunder_postion(contents):
    """Raise if dunders are after imports."""
    dunders = find_dunders(contents)
    if not dunders:
        raise MissingDunderError()

    first_import = find_first_import_line(contents)
    faulty_dunders = [dunder for dunder, line in dunders.items()
                      if line > first_import]
    if faulty_dunders:
        raise MisplacedDunderError(*faulty_dunders)


def check_duplicate_dunders(contents):
    """Raise if a dunder is found more then once."""
    dunder_lines = find_dunders(contents)
    dunders = {lineno: DUNDER_EQUALS_RE.match(line)['dunder']
               for line, lineno in dunder_lines.items()}
    dunder_counts = Counter(dunders.values())
    duplicates = {d: c for d, c in dunder_counts.items() if c>1}
    if duplicates:
        raise DuplicateDundersError(duplicates)


def check_dunders_sorted(contents):
    """Raise if dunders are not sorted alphabetically."""
    dunders = [DUNDER_EQUALS_RE.match(line)['dunder']
               for line in find_dunders(contents)]
    if dunders != sorted(dunders):
        raise DunderSortingError


def check_dunders_format(contents):
    """Ensure all dunders have the expected format."""
    ast_dunders = ast_parse_dunder_lines(contents)
    check_dunder_types(ast_dunders)
    dunder_line_ranges = find_dunder_line_ranges(contents)                      # TODO: check single-quote strings; tuples have one string element per line, ends with comma
    check_no_empty_line_in_dunders(dunder_line_ranges, contents)
    check_dunders_contiguous(dunder_line_ranges)


def check_dunder_types(ast_dunders):
    """Check that all ast_dunders have ACCEPTABLE_DUNDER_TYPES."""
    invalid = [(d, node) for d, node in ast_dunders.items()
               if not isinstance(node, ACCEPTABLE_DUNDER_TYPES)]
    if invalid:
        err_ = (
            'Found non-tuple/string dunders: '
            + ', '.join(f'{d}: {type(node).__name__}' for d, node in invalid)
            )
        raise DunderFormatError(err_)


def check_no_empty_line_in_dunders(dunder_line_ranges, contents):
    """Raise if a dunder has an empty line in its definition."""
    invalid = []
    for dunder, (start, end) in dunder_line_ranges.items():
        if end == start:
            continue
        dunder_def = contents[start:end+1]
        if any(not line.strip() for line in dunder_def):
            invalid.append(dunder)
    if invalid:
        raise DunderFormatError('Empty line(s) found in dunder definition: '
                                + ', '.join(invalid))


def check_dunders_contiguous(dunder_line_ranges):
    """Raise if there is any empty line between two dunders."""
    invalid = []
    dunder_ranges = tuple(dunder_line_ranges.values())
    dunder_range_pairs = zip(dunder_ranges[:-1], dunder_ranges[1:])
    for (_, end_first), (begin_second, _) in dunder_range_pairs:
        if begin_second != end_first + 1:
            invalid.append(f'L{end_first}-{begin_second}')
    if invalid:
        raise DunderFormatError('Non contiguous dunders: '+ ', '.join(invalid))


def check_one_file(filename):
    """Check that one file is OK with its dunders."""
    with filename.open('r', encoding='utf-8') as file:
        contents = file.readlines()
    if not contents:  # Don't bother with empty files
        return
    check_dunder_postion(contents)
    check_duplicate_dunders(contents)
    check_dunders_sorted(contents)
    check_dunders_format(contents)


def fix_one_file(filename):
    """Fix dynamic dunders in filename."""
    added_copyright = add_copyright_dunder(filename)
    added_license = add_license_dunder(filename)
    return added_copyright or added_license


def add_one_dunder(filename, new_dunder, value):
    """Add a new_dunder with value if not present in filename."""
    with filename.open('r', encoding='utf-8') as file:
        contents = file.readlines()
    if not contents:
        return
    current_dunders = find_dunder_line_ranges(contents)
    try:
        lineno, _ = current_dunders[new_dunder]
    except KeyError:  # Not there
        pass
    else:  # Found                                                              # TODO: edit if necessary
        return False
    lineno = find_new_dunder_position(current_dunders, new_dunder)
    contents.insert(lineno, f'{new_dunder} = {value!r}\n')
    with filename.open('w', encoding='utf-8') as file:
        file.writelines(contents)
    return True


def find_new_dunder_position(current_dunders, new_dunder):
    """Find the right line for new_dunder."""
    new_dunders = sorted((*current_dunders, new_dunder))
    new_dunder_idx = new_dunders.index(new_dunder)
    try:  # Should it go before an old one?
        old_dunder_after = new_dunders[new_dunder_idx+1]
    except IndexError:  # Should go after
        pass
    else:
        old_dunder_start, _ = current_dunders[old_dunder_after]
        return old_dunder_start
    # Must go after an old one
    old_dunder_before = new_dunders[new_dunder_idx-1]
    _, old_dunder_end = current_dunders[old_dunder_before]
    return old_dunder_end + 1


def add_license_dunder(filename):
    """Add a __license__ = 'GPLv3+' dunder."""
    return add_one_dunder(filename, '__license__', 'GPLv3+')


def add_copyright_dunder(filename):
    """Add a __copyright__ dunder."""
    _copyright = 'Copyright (c) 2019-2024 ViPErLEED developers'
    return add_one_dunder(filename, '__copyright__', _copyright)


def check_and_fix_viperleed_repo(package):
    """Check viperleed repository, and fix files that should be fixed."""
    base = Path.cwd().parents[1]
    repo = base / package
    print(f'\nChecking {repo}')
    for file in repo.glob('**/*.py'):
        try:
            check_one_file(file)
        except DunderError as exc:
            print(f'Problem in {file.relative_to(base)}: {exc}')
            continue
        fixed = fix_one_file(file)
        if fixed:
            print(f'Fixed {file.relative_to(base)}')


if __name__ == '__main__':
    check_and_fix_viperleed_repo('src/viperleed')
    check_and_fix_viperleed_repo('tests')
    check_and_fix_viperleed_repo('build/scripts')
