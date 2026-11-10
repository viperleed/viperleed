"""Module qt_codemod for viperleed files.

Rewrites Python code written for PyQt5 so that it runs with both PyQt5
and PyQt6 via viperleed.gui.qt_compat. In particular, it:
- imports Qt modules from viperleed.gui.qt_compat rather than from
  PyQt5 or PyQt6;
- turns short enumeration members into fully scoped ones, e.g.,
  qtc.Qt.AlignLeft -> qtc.Qt.AlignmentFlag.AlignLeft;
- takes the classes that PyQt6 moved from QtWidgets to QtGui (e.g.,
  QAction) from QtGui, adding an import of QtGui if needed;
- renames exec_ to exec, both where it is used and where it is
  defined (e.g., in subclasses that override it), and replaces
  QtWidgets.qApp with QtWidgets.QApplication.instance().

Only attribute chains that start at a name imported from a Qt module
(e.g., qtc, qtw, Qt, QMessageBox) can be handled, as well as short
enumeration members accessed via self in methods of classes that
derive from Qt classes (e.g., self.Ok in a QMessageBox subclass).
Everything else that needs attention (e.g., removed classes, ambiguous
enumeration members, lines that became too long) is reported, but left
untouched. The rest of each file is kept byte by byte, including its
line endings. Running the script again on its output changes nothing.
With --strings, strings mentioning PyQt5, PyQt6, or Qt5 (e.g., in
docstrings or in calls to importlib.import_module) are also reported.

Requires Python 3.8 or later. The fully scoped names are taken from
PyQt6, which must be installed. Use --verify in an environment with
PyQt5 (or PyQt6) to check that all the Qt attribute chains in the
files resolve with that binding.

Usage:
    python qt_codemod.py [--check | --verify] [--strings] PATH [...]

PATH may be a file or a folder (searched recursively for *.py files).
Files of the viperleed.gui.qt_compat package are skipped. The exit
code is 1 if --check finds files to change, or if anything is
reported (including chains that do not resolve with --verify). It is
0 otherwise.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-09'
__license__ = 'GPLv3+'

import argparse
import ast
from collections import namedtuple
import enum
import importlib
import os
from pathlib import Path
import sys

COMPAT = 'viperleed.gui.qt_compat'
BINDINGS = ('PyQt5', 'PyQt6')
QTCORE = 'QtCore'
QTGUI = 'QtGui'
QTWIDGETS = 'QtWidgets'
QT_MODULES = (QTCORE, QTGUI, 'QtNetwork', 'QtSerialPort', QTWIDGETS)
QTGUI_ALIAS = 'qtg'
MOVED_TO_QTGUI = (
    'QAction',
    'QActionGroup',
    'QFileSystemModel',
    'QShortcut',
    'QUndoCommand',
    'QUndoGroup',
    'QUndoStack',
    )
QAPP = 'qApp'
RENAMED_METHODS = {'exec_': 'exec'}
QT_STRING_MARKERS = ('PyQt5', 'PyQt6', 'qt5agg', 'Qt5Agg')
MAX_LINE_LENGTH = 79
CRLF = '\r\n'

# A replacement of the text between two (line, column) positions
Edit = namedtuple('Edit', ('start', 'end', 'text'))

# An attribute chain rooted at a Name. nodes has the Name node followed
# by the Attribute nodes, from the innermost to the outermost. names
# are the identifiers along the chain.
Chain = namedtuple('Chain', ('nodes', 'names'))

# A Qt attribute chain that does not resolve. path are the attribute
# names starting from the Qt module, index is the position in path of
# the first missing attribute, owner the object that misses it.
Missing = namedtuple('Missing', ('chain', 'qt_module', 'path', 'index',
                                 'owner', 'n_prefix'))


class QtFile:
    """A Python file that may use Qt."""

    def __init__(self, path, qt_modules):
        """Initialize instance.

        Parameters
        ----------
        path : Path
            The file to be processed.
        qt_modules : dict
            Keys are names of Qt modules (e.g., 'QtCore'), values the
            corresponding modules that are used to resolve names.
        """
        self.path = path
        self.qt_modules = qt_modules
        self.source = path.read_bytes().decode('utf-8')
        self.tree = ast.parse(self.source, filename=str(path))
        self.edits = []
        self.reports = []
        # {name in file: (qt module, (attributes, ...))}
        self.roots = {}

    @property
    def line_ending(self):
        """Return the line ending used in this file."""
        return CRLF if CRLF in self.source else '\n'

    @property
    def lines(self):
        """Return the lines of the source, with line endings."""
        return self.source.splitlines(keepends=True)

    def report(self, node, msg):
        """Store a message about `node`."""
        self.reports.append(f'{self.path}:{node.lineno}: {msg}')

    # ---------------------------- imports ----------------------------
    def find_imports(self, rewrite):
        """Collect names imported from Qt modules, rewriting imports."""
        for node in ast.walk(self.tree):
            if isinstance(node, ast.ImportFrom):
                self._handle_import_from(node, rewrite)
            elif isinstance(node, ast.Import):
                self._handle_import(node, rewrite)

    def _handle_import_from(self, node, rewrite):
        """Process a 'from X import Y' statement."""
        module = node.module or ''
        if module in (*BINDINGS, COMPAT):        # from PyQt5 import QtCore
            for alias in node.names:
                if alias.name in self.qt_modules:
                    name = alias.asname or alias.name
                    self.roots[name] = (alias.name, ())
            if rewrite and module != COMPAT:
                self._replace_module_name(node, module)
            return
        package, _, qt_module = module.rpartition('.')
        if package not in (*BINDINGS, COMPAT):
            return
        if qt_module not in self.qt_modules:
            if package != COMPAT:     # Ignore, e.g., qt_compat.binding
                self.report(node, f'Unknown Qt module {module}')
            return
        for alias in node.names:                 # from PyQt5.QtCore import Qt
            self.roots[alias.asname or alias.name] = (qt_module,
                                                      (alias.name,))
        if rewrite and package != COMPAT:
            self._replace_module_name(node, package)

    def _handle_import(self, node, rewrite):
        """Process an 'import X' statement."""
        for alias in node.names:
            package, _, qt_module = alias.name.rpartition('.')
            if alias.name in BINDINGS or (package in BINDINGS
                                          and not alias.asname):
                self.report(node, f'Import of {alias.name} without an '
                                  'alias: please fix manually')
                continue
            if package not in BINDINGS:
                continue
            self.roots[alias.asname] = (qt_module, ())
            if not rewrite:
                continue
            if len(node.names) > 1:
                self.report(node, 'Multiple imports in one statement: '
                                  'please fix manually')
                continue
            # Rewrite, e.g., the PyQt5 QtCore import as an import of
            # QtCore from qt_compat, keeping the same alias
            self._replace_node(
                node,
                f'from {COMPAT} import {qt_module} as {alias.asname}'
                )

    def _replace_module_name(self, node, old):
        """Replace `old` in the first line of `node` with COMPAT."""
        line = self.lines[node.lineno - 1]
        col = line.find(old)
        if col < 0:
            self.report(node, f'Could not find {old!r} on the first line '
                              'of the import statement')
            return
        self.edits.append(Edit((node.lineno, col),
                               (node.lineno, col + len(old)),
                               COMPAT))

    def _replace_node(self, node, text):
        """Replace the whole source of `node` with `text`."""
        self.edits.append(Edit(self._start(node), self._end(node), text))

    # ------------------------- node positions ------------------------
    def _to_char_col(self, line_no, byte_col):
        """Return the character column of a UTF-8 byte offset."""
        line = self.lines[line_no - 1]
        return len(line.encode('utf-8')[:byte_col].decode('utf-8'))

    def _start(self, node):
        """Return the (line, column) where `node` starts."""
        return node.lineno, self._to_char_col(node.lineno, node.col_offset)

    def _end(self, node):
        """Return the (line, column) where `node` ends."""
        return (node.end_lineno,
                self._to_char_col(node.end_lineno, node.end_col_offset))

    # ------------------------ attribute chains ------------------------
    def _iter_chains(self, tree=None):
        """Yield a Chain for each maximal chain rooted at a Name."""
        tree = tree or self.tree
        inner_attributes = {
            id(node.value) for node in ast.walk(tree)
            if isinstance(node, ast.Attribute)
            }
        for node in ast.walk(tree):
            if (not isinstance(node, ast.Attribute)
                    or id(node) in inner_attributes):
                continue
            nodes = [node]
            while isinstance(nodes[-1], ast.Attribute):
                nodes.append(nodes[-1].value)
            if not isinstance(nodes[-1], ast.Name):
                continue
            nodes.reverse()
            names = [nodes[0].id, *(n.attr for n in nodes[1:])]
            yield Chain(nodes, names)

    def rename_methods(self):
        """Rename methods such as exec_, regardless of the object.

        Both uses (e.g., obj.exec_()) and definitions (e.g., in
        subclasses overriding exec_) are renamed.
        """
        for node in ast.walk(self.tree):
            if (isinstance(node, ast.Attribute)
                    and node.attr in RENAMED_METHODS):
                end = self._end(node)
                start = (end[0], end[1] - len(node.attr))
                self.edits.append(Edit(start, end,
                                       RENAMED_METHODS[node.attr]))
            elif (isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                    and node.name in RENAMED_METHODS):
                self._rename_definition(node)

    def _rename_definition(self, node):
        """Rename the function defined by `node`."""
        line = self.lines[node.lineno - 1]
        col = line.find(f'def {node.name}(')
        if col < 0:
            self.report(node, f'Could not rename definition of '
                              f'{node.name}: please fix manually')
            return
        col += len('def ')
        self.edits.append(Edit((node.lineno, col),
                               (node.lineno, col + len(node.name)),
                               RENAMED_METHODS[node.name]))

    def fix_chains(self, rewrite):
        """Make Qt attribute chains valid, or report them.

        Parameters
        ----------
        rewrite : bool
            If True, chains are fixed with PyQt6 names. If False, the
            chains are only checked against the Qt modules.
        """
        for chain in self._iter_chains():
            missing = self._find_missing(chain, rewrite)
            if not missing:
                continue
            if rewrite:
                self._fix_missing(missing)
                continue
            dotted = '.'.join(missing.path[:missing.index + 1])
            self.report(chain.nodes[0],
                        f'{missing.qt_module}.{dotted} does not exist')

    def fix_self_enums(self, rewrite):
        """Make short enum members accessed via self fully scoped.

        Only methods of classes whose bases are Qt classes (directly,
        or via other classes defined in the same file) are considered.
        Attributes that this file assigns to self are left alone.
        """
        assigned = {
            node.attr for node in ast.walk(self.tree)
            if isinstance(node, ast.Attribute)
            and isinstance(node.ctx, ast.Store)
            }
        for class_node, qt_base in self._iter_qt_classes():
            for chain in self._iter_chains(class_node):
                if chain.names[0] not in ('self', 'cls'):
                    continue
                attr = chain.names[1]
                if attr in assigned or hasattr(qt_base, attr):
                    continue
                enum_name = self._find_enum(qt_base, attr, chain.nodes[1])
                if not enum_name:
                    continue
                if not rewrite:
                    self.report(chain.nodes[1],
                                f'{chain.names[0]}.{attr} does not exist')
                    continue
                self_end = self._end(chain.nodes[0])
                self.edits.append(Edit(self_end, self_end,
                                       f'.{enum_name}'))

    def _iter_qt_classes(self):
        """Yield (ClassDef, Qt base) for each subclass of a Qt class."""
        qt_bases = {}   # {class name: Qt base class}
        for node in ast.walk(self.tree):
            if not isinstance(node, ast.ClassDef):
                continue
            for base in node.bases:
                qt_base = self._resolve_base(base, qt_bases)
                if qt_base is not None:
                    qt_bases[node.name] = qt_base
                    yield node, qt_base
                    break

    def _resolve_base(self, base, qt_bases):
        """Return the Qt class that `base` refers to, or None."""
        if isinstance(base, ast.Name) and base.id in qt_bases:
            return qt_bases[base.id]
        names = []
        while isinstance(base, ast.Attribute):
            names.append(base.attr)
            base = base.value
        if not isinstance(base, ast.Name) or base.id not in self.roots:
            return None
        qt_module, prefix = self.roots[base.id]
        obj = self.qt_modules[qt_module]
        for attr in (*prefix, *reversed(names)):
            obj = getattr(obj, attr, None)
        return obj if isinstance(obj, type) else None

    def _find_missing(self, chain, rewrite):
        """Return a Missing for `chain`, or None if it resolves."""
        root = chain.names[0]
        if root not in self.roots:
            return None
        qt_module, prefix = self.roots[root]
        path = [*prefix, *chain.names[1:]]
        owner = self.qt_modules[qt_module]
        for index, attr in enumerate(path):
            if attr in RENAMED_METHODS and not rewrite:
                return None  # Will be renamed anyway
            try:
                owner = getattr(owner, attr)
            except AttributeError:
                return Missing(chain, qt_module, path, index, owner,
                               len(prefix))
        return None

    def _fix_missing(self, missing):
        """Fix a chain with a missing attribute, or report it."""
        nodes = missing.chain.nodes
        attr = missing.path[missing.index]
        # nodes[0] is the root Name. The node whose attribute is
        # missing comes after the len(prefix) attributes in the root.
        node_index = missing.index - missing.n_prefix + 1
        node = nodes[max(node_index, 0)]
        dotted = '.'.join(missing.path[:missing.index + 1])
        dotted = f'{missing.qt_module}.{dotted}'
        is_module_attr = not missing.index and node_index == 1
        if missing.qt_module == QTWIDGETS and attr in MOVED_TO_QTGUI:
            if not is_module_attr:
                self.report(node, f'{dotted} moved to QtGui: please fix '
                                  'manually')
                return
            self._replace_node(nodes[0], self._ensure_qtgui_import())
            return
        if missing.qt_module == QTWIDGETS and attr == QAPP:
            if not is_module_attr:
                self.report(node, f'{dotted}: please fix manually')
                return
            self._replace_node(nodes[1],
                               f'{nodes[0].id}.QApplication.instance()')
            return
        enum_name = self._find_enum(missing.owner, attr, node)
        if not enum_name or node_index < 1:
            self.report(node, f'{dotted} does not exist in PyQt6')
            return
        # Insert '.EnumName' right after the owner's node
        owner_end = self._end(nodes[node_index - 1])
        self.edits.append(Edit(owner_end, owner_end, f'.{enum_name}'))

    def _find_enum(self, owner, attr, node):
        """Return the name of the enumeration of `owner` with `attr`."""
        candidates = []
        for name in dir(owner):
            if name.startswith('_'):
                continue
            value = getattr(owner, name, None)
            if (isinstance(value, type) and issubclass(value, enum.Enum)
                    and attr in value.__members__):
                candidates.append(name)
        if len(candidates) > 1:
            owner_name = getattr(owner, '__name__', owner)
            enums = ', '.join(candidates)
            self.report(node, f'{owner_name}.{attr} is ambiguous among '
                              f'{enums}: please fix manually')
            return None
        return candidates[0] if candidates else None

    def _ensure_qtgui_import(self):
        """Return the name used for QtGui, importing QtGui if needed."""
        for name, (qt_module, prefix) in self.roots.items():
            if qt_module == QTGUI and not prefix:
                return name
        if QTGUI_ALIAS in self.roots:
            self.reports.append(f'{self.path}: name {QTGUI_ALIAS!r} is '
                                'used for something else. Please import '
                                'QtGui manually')
            return QTGUI_ALIAS
        # Insert right before the import of QtWidgets
        for node in self.tree.body:
            if (isinstance(node, (ast.Import, ast.ImportFrom))
                    and any(self.roots.get(a.asname or a.name)
                            == (QTWIDGETS, ()) for a in node.names)):
                start = (node.lineno, 0)
                import_qtgui = f'from {COMPAT} import QtGui as {QTGUI_ALIAS}'
                self.edits.append(Edit(start, start,
                                       import_qtgui + self.line_ending))
                self.roots[QTGUI_ALIAS] = (QTGUI, ())
                return QTGUI_ALIAS
        self.reports.append(f'{self.path}: could not add an import of '
                            'QtGui. Please do it manually')
        return QTGUI_ALIAS

    # ----------------------------- output -----------------------------
    def find_qt_strings(self):
        """Report strings that mention PyQt5/PyQt6 or Qt5."""
        for node in ast.walk(self.tree):
            if (not isinstance(node, ast.Constant)
                    or not isinstance(node.value, str)):
                continue
            for i, line in enumerate(node.value.splitlines()):
                marker = next((m for m in QT_STRING_MARKERS if m in line),
                              None)
                if marker:
                    self.reports.append(
                        f'{self.path}:{node.lineno + i}: string mentions '
                        f'{marker}, check manually: {line.strip()!r}'
                        )

    def new_source(self):
        """Return the source with all edits applied."""
        offsets = [0]
        for line in self.lines:
            offsets.append(offsets[-1] + len(line))

        def _absolute(position):
            line_no, col = position
            return offsets[line_no - 1] + col

        source = self.source
        for edit in sorted(set(self.edits),
                           key=lambda e: (_absolute(e.start),
                                          _absolute(e.end)),
                           reverse=True):
            start, end = _absolute(edit.start), _absolute(edit.end)
            source = source[:start] + edit.text + source[end:]
        return source

    def report_long_lines(self, new_source):
        """Report lines of `new_source` that became too long."""
        old_lines = set(self.source.splitlines())
        for line_no, line in enumerate(new_source.splitlines(), start=1):
            if len(line) > MAX_LINE_LENGTH and line not in old_lines:
                self.reports.append(f'{self.path}:{line_no}: line longer '
                                    f'than {MAX_LINE_LENGTH} characters')


def _iter_files(paths):
    """Yield all Python files in `paths`, except those of qt_compat."""
    compat_folder = COMPAT.replace('.', '/')
    for path_str in paths:
        path = Path(path_str)
        files = sorted(path.rglob('*.py')) if path.is_dir() else (path,)
        yield from (f for f in files
                    if compat_folder not in f.resolve().as_posix())


def _import_qt_modules(package):
    """Return a dict of the Qt modules in `package`."""
    modules = {}
    for name in QT_MODULES:
        try:
            modules[name] = importlib.import_module(f'{package}.{name}')
        except ImportError as exc:
            sys.exit(f'Cannot import {package}.{name}: {exc}')
    return modules


def _parse_args(argv):
    """Return the command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Make PyQt5 code work with PyQt5 and PyQt6.'
        )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--check', action='store_true',
                      help='Do not write files. Report what would change.')
    mode.add_argument('--verify', action='store_true',
                      help=('Do not write files. Check that Qt attribute '
                            'chains resolve with the installed binding.'))
    parser.add_argument('--strings', action='store_true',
                        help=('Also report strings (e.g., docstrings) '
                              'that mention PyQt5, PyQt6, or Qt5.'))
    parser.add_argument('paths', nargs='+', metavar='PATH')
    return parser.parse_args(argv)


def _process(path, args, qt_modules):
    """Process one file. Return whether it was (or would be) changed."""
    rewrite = not args.verify
    qt_file = QtFile(path, qt_modules)
    qt_file.find_imports(rewrite)
    if args.strings:
        qt_file.find_qt_strings()
    if rewrite:
        qt_file.rename_methods()
    if qt_file.roots:
        qt_file.fix_chains(rewrite)
        qt_file.fix_self_enums(rewrite)
    new_source = qt_file.new_source() if rewrite else qt_file.source
    qt_file.report_long_lines(new_source)
    for line in qt_file.reports:
        print(line)
    changed = new_source != qt_file.source
    if changed and args.check:
        print(f'{path}: would be changed')
    elif changed:
        path.write_bytes(new_source.encode('utf-8'))
    return changed, bool(qt_file.reports)


def main(argv=None):
    """Run the script. Return the exit code."""
    if sys.version_info < (3, 8):
        sys.exit('qt_codemod requires Python 3.8 or later')
    args = _parse_args(argv)
    if args.verify:
        qt_modules = _import_qt_modules(COMPAT)
    else:
        os.environ.setdefault('QT_API', 'pyqt6')
        qt_modules = _import_qt_modules('PyQt6')
    results = [_process(path, args, qt_modules)
               for path in _iter_files(args.paths)]
    n_changed = sum(changed for changed, _ in results)
    if not args.verify:
        action = 'would change' if args.check else 'changed'
        print(f'{action} {n_changed} file(s)', file=sys.stderr)
    any_reports = any(reported for _, reported in results)
    return int(any_reports or (args.check and n_changed > 0))


if __name__ == '__main__':
    sys.exit(main())
