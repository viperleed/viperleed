"""Script check_import_order.

Provides functionality for checking the sort order of imports in all
viperleed files.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-03-24'
__license__ = 'GPLv3+'

import ast
from collections import Counter, defaultdict
from dataclasses import dataclass
from enum import IntEnum
from operator import itemgetter
from pathlib import Path

from isort import place_module
from isort.sections import STDLIB


# TODO: complain if there are too many imports in a single statement
#       suggesting to make multiple statements instead
# TODO: check exactly one empty line between groups
# TODO: check no __dunder__ module import

_ALL_IMPORTS = 'all'  # Just a dict key for all the imports
HELP_ME_I_CANT_FIGURE_OUT_THE_RIGHT_ORDER = False


class _BaseImportError(Exception):
    """Base exception for all import-related errors."""

    _default_message =''

    def __init__(self, message=''):
        """Initialize exception."""
        super().__init__(message or self._default_message)


class DuplicateImportError(_BaseImportError):
    """A name was imported more than once."""

    def __init__(self, *duplicates):
        """Initialize instance."""
        super().__init__('Found duplicate imports: ' + ', '.join(duplicates))


class MisplacedImports(_BaseImportError):
    """Found some imports in the wrong group."""

    _default_message = ('Imports are not sorted as stdlib, third-party, '
                        'viperleed, relative')


class SortOrderError(_BaseImportError):
    """Found some sorting fuck-up."""

    def __init__(self, group, message=''):
        """Initialize exception."""
        message = message or (f'{group}: {self._default_message}')
        super().__init__(message)


class ImportSortOrderError(SortOrderError):
    """Stuff imported from a module is not sorted."""

    _default_message = 'Imported names are not alphabetically sorted'


class ModulesSortOrderError(SortOrderError):
    """Modules are in the wrong sort order."""

    _default_message = 'Imported modules are not alphabetically sorted'


class RelativeImportsLevelSortError(_BaseImportError):
    """Relative imports are not sorted from highest to lower in the tree."""

    _default_message = ('Relative imports not sorted from higher '
                        'to lower level (e.g., "..." before "..")')


class RelativeImportsNameAfterModuleError(_BaseImportError):
    """Found 'from .module import name' before 'from . import module'."""

    _default_message = ('"from .module import name" must '
                        'be after "from . import module"')


class StarImportError(_BaseImportError):
    """Found 'import *'."""

    _default_message = 'Found "import *"'


class TooManyImportsInOneLine(_BaseImportError):
    """Found a line as 'import x, y'."""

    _default_message = ('Should not have multiple comma-'
                        'separated imports in one line')


class ImportType(IntEnum):  # Sorted as expected
    """A classification for import lines."""
    STDLIB = 1
    THIRD_PARTY = 2
    VIPERLEED = 3
    RELATIVE = 4


@dataclass
class ImportStatement:
    """Collection of information concerning an import line."""
    module : str             # Where stuff is imported from
    imported : tuple = ()    # What is imported (only for 'from').
    level: int = 0           # How many 'dots'


def check_duplicate_import(module_name_pairs):
    """Check if there is any duplication among the imported stuff."""
    counts = Counter(module_name_pairs)
    # print(counts)
    duplicates = [f'{module}.{name}' if module != name else name
                  for (module, name), repeats in counts.items()
                  if repeats > 1]
    if duplicates:
        raise DuplicateImportError(*duplicates)


def check_one_file(filename):
    """Check sort order of imports in filename."""
    with filename.open('r', encoding='utf-8') as file:
        contents = file.read()
    if not contents:  # Don't bother with empty files
        return
    imports = classify_imports(contents)
    all_names = get_all_imported_names(imports)
    check_star_import(all_names)
    check_duplicate_import(flat_imports(imports))
    check_sort_order_among_groups(imports)
    for import_type, imports_by_type in imports.items():
        if import_type == _ALL_IMPORTS:
            continue
        if import_type is ImportType.RELATIVE:
            check_relative_imports_sort_order(imports_by_type)
            continue
        check_simple_group_sort_order(import_type.name, imports_by_type)


def check_relative_imports_sort_order(relative_imports):
    """Raise if relative imports are not correctly sorted."""
    # First descending level
    levels = [i.level for i in relative_imports]
    if not is_sorted(levels, reverse=True):
        raise RelativeImportsLevelSortError

    # Then sorting between same-level imports
    by_level = defaultdict(list)  # {level: [imports]}
    for import_line in relative_imports:
        by_level[import_line.level].append(import_line)

    for level, imports in by_level.items():
        # Split into "from . import module" and "from .module import".
        # The first ones should come earlier than the latter ones
        modules = [stmt for stmt in imports if not stmt.module]
        names = [i for i in imports if i.module]
        check_relative_module_imports_before_name_imports(imports,
                                                          modules,
                                                          names)
        # Then sorted within each group
        check_simple_group_sort_order('"from ' + '.'*level + '" imports',
                                      modules)
        check_simple_group_sort_order('"from ' + '.'*level + 'module" imports',
                                      names)


def check_relative_module_imports_before_name_imports(imports, modules, names):
    """Raise if any of the imported names comes before the modules."""
    if not imports or not modules or not names:
        return
    if imports.index(modules[-1]) > imports.index(names[0]):
        raise RelativeImportsNameAfterModuleError


def check_simple_group_sort_order(group, imports):                              # TODO: treat single-leading-underscore modules as if they had no underscore
    """Check sort order of a group of imports."""
    modules = [statement.module for statement in imports]
    if not is_sorted(modules):
        raise ModulesSortOrderError(group)
    # Now collate all the imports from one module,
    # and make sure things are alphabetically ordered
    imported = defaultdict(list)
    for statement in imports:
        imported[statement.module].extend(statement.imported)
    invalid_modules = []
    for module, all_imports in imported.items():
        if not is_sorted(all_imports):
            invalid_modules.append(module)
            if HELP_ME_I_CANT_FIGURE_OUT_THE_RIGHT_ORDER:
                print(f'   {module}: {sorted(all_imports)}')
    if invalid_modules:
        raise ImportSortOrderError(', '.join(invalid_modules) or group)


def check_sort_order_among_groups(imports):
    """Raise if imports are not sorted in groups as expected."""
    # Sort by the first item, i.e., the name type
    if not is_sorted(imports[_ALL_IMPORTS], key=itemgetter(0)):
        raise MisplacedImports


def check_star_import(imported_names):
    """Raise if there is 'import *'."""
    if '*' in imported_names:
        raise StarImportError


def classify_imports(contents):
    """Return imports in contents, split by type."""
    all_imports = {t: [] for t in ImportType}  # Correct group order
    all_imports[_ALL_IMPORTS] = []        # For sorting among groups
    ast_tree = ast.parse(contents)
    for node in ast_tree.body:
        if isinstance(node, ast.ImportFrom):
            import_line = _handle_import_from_node(node)
        elif isinstance(node, ast.Import):
            import_line = _handle_simple_import_node(node)
        else:  # Not an import
            continue
        if import_line.level:
            key = ImportType.RELATIVE
        elif is_stdlib(import_line.module):
            key = ImportType.STDLIB
        elif 'viperleed' in import_line.module:
            key = ImportType.VIPERLEED
        else:
            key = ImportType.THIRD_PARTY
        all_imports[key].append(import_line)
        all_imports[_ALL_IMPORTS].append((key, import_line))
    return all_imports


def flat_imports(classified_imports):
    """Yield pairs of (module, name) for all the single imports."""
    for _, import_line in classified_imports[_ALL_IMPORTS]:
        imported = import_line.imported or (import_line.module,)
        for name in imported:
            yield import_line.module, name


def get_all_imported_names(classified_imports):
    """Return a list of all the imported in classified_imports."""
    all_names = []
    for _, import_line in classified_imports[_ALL_IMPORTS]:
        all_names.extend(import_line.imported)
    return all_names


def _handle_simple_import_node(node):
    """Return an ImportStatement from a node of type 'import xxx'."""
    alias, *other_names = node.names
    if other_names:
        raise TooManyImportsInOneLine
    return ImportStatement(alias.name)  # Don't care about 'as'


def _handle_import_from_node(node):
    """Return an ImportStatement from a node of type 'from xxx import yyy'."""
    names = tuple(alias.name for alias in node.names)
    return ImportStatement(node.module or '', names, node.level)


def is_sorted(input_list, **kwargs):
    """Return whether input_list is sorted."""
    return input_list == sorted(input_list, **kwargs)


def is_stdlib(module):
    """Return whether module is a stdlib one."""
    return place_module(module) == STDLIB


def check_sort_order_in_repo(package):
    """Check sort order in all files in repo."""
    base = Path.cwd().parents[1]
    repo = base / package
    print(f'\nChecking {repo}')
    for file in repo.glob('**/*.py'):
        try:
            check_one_file(file)
        except _BaseImportError as exc:
            print(f'Problem in {file.relative_to(base)}: {exc}')
            if HELP_ME_I_CANT_FIGURE_OUT_THE_RIGHT_ORDER:
                print()


if __name__ == '__main__':
    check_sort_order_in_repo('viperleed')
    check_sort_order_in_repo('tests')
    check_sort_order_in_repo('build/scripts')