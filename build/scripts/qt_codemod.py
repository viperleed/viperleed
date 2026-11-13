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
Short enumeration members accessed via other objects (e.g., palette.Text
or bar.tabBar().SelectLeftTab) are reported, as the type of the object
is unknown. With --fix-instance-enums, they are also fixed, guessing the
enumeration from the name of the object if the member is ambiguous
(e.g., palette.Text -> palette.ColorRole.Text, not OpenModeFlag). These
changes should be reviewed. A few members removed in Qt 6 (e.g.,
QPalette.Background) are replaced as well. Everything else that needs
attention (e.g., removed classes, ambiguous enumeration members, lines
that became too long) is reported, but left untouched. The rest of
each file is kept byte by byte, including its line endings. Running
the script again on its output changes nothing.
With --strings, strings mentioning PyQt5, PyQt6, or Qt5 (e.g., in
docstrings or in calls to importlib.import_module) are also reported.

Requires Python 3.8 or later. The fully scoped names are taken from
PyQt6, which must be installed. Use --verify in an environment with
PyQt5 (or PyQt6) to check that all the Qt attribute chains in the
files resolve with that binding. With PyQt5, --verify also reports
short members of PyQt5 enumerations accessed via objects, including
those that PyQt6 removed.

Use --compare REV in an environment with PyQt5 to check that the
rewritten code means the same as it did at git revision REV (e.g.,
the commit before running the script): all Qt attribute chains that
changed must resolve to the same objects (e.g., Qt.AlignLeft and
Qt.AlignmentFlag.AlignLeft). Short enum members accessed via objects
(see --fix-instance-enums) are not checked.

Usage:
    python qt_codemod.py [--check | --verify | --compare REV]
                         [--strings] [--fix-instance-enums] PATH [...]

PATH may be a file or a folder (searched recursively for *.py files).
Files of the viperleed.gui.qt_compat package are skipped. The exit
code is 1 if --check finds files to change, or if anything is
reported (including chains that do not resolve with --verify, and
chains that changed meaning with --compare). It is 0 otherwise.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-09'
__license__ = 'GPLv3+'

import argparse
import ast
from collections import Counter, namedtuple
import difflib
import enum
import importlib
import os
from pathlib import Path
import subprocess
import sys

BINDINGS = ('PyQt5', 'PyQt6')
COMPAT = 'viperleed.gui.qt_compat'
CRLF = '\r\n'
MAX_LINE_LENGTH = 79
# Shortest name of an object used to guess its Qt class (e.g., 'pal'
# would match QPalette, while 'p' would match too many classes)
MIN_GUESS_LENGTH = 3
MISSING = object()  # Marks Qt attributes that do not exist
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
QTCORE = 'QtCore'
QTGUI = 'QtGui'
QTGUI_ALIAS = 'qtg'
QTWIDGETS = 'QtWidgets'
QT_MODULES = (QTCORE, QTGUI, 'QtNetwork', 'QtSerialPort', QTWIDGETS)
QT_STRING_MARKERS = ('PyQt5', 'PyQt6', 'qt5agg', 'Qt5Agg')
# Enumeration members removed in Qt 6, with their replacement
Removed = namedtuple('Removed', ('qt_class', 'enum', 'replacement'))
REMOVED_MEMBERS = {
    'Background': Removed('QPalette', 'ColorRole', 'Window'),
    'Foreground': Removed('QPalette', 'ColorRole', 'WindowText'),
    'MidButton': Removed('Qt', 'MouseButton', 'MiddleButton'),
    }
RENAMED_METHODS = {'exec_': 'exec'}
SIP_ENUM_TYPE = 'enumtype'  # The metaclass of PyQt5 enumerations


def main(argv=None):
    """Run the script. Return the exit code."""
    if sys.version_info < (3, 8):
        sys.exit('qt_codemod requires Python 3.8 or later')
    args = _parse_args(argv)
    if args.compare:
        qt_modules = _import_qt_modules(COMPAT)
        messages = [message
                    for path in _iter_files(args.paths)
                    for message in _compare(path, args.compare, qt_modules)]
        for message in messages:
            print(message)
        return int(bool(messages))
    if args.verify:
        qt_modules = _import_qt_modules(COMPAT)
    else:
        os.environ.setdefault('QT_API', 'pyqt6')
        qt_modules = _import_qt_modules('PyQt6')
    enum_index = _build_enum_index(qt_modules)
    results = [_process(path, args, qt_modules, enum_index)
               for path in _iter_files(args.paths)]
    n_changed = sum(changed for changed, _ in results)
    if not args.verify:
        action = 'would change' if args.check else 'changed'
        print(f'{action} {n_changed} file(s)', file=sys.stderr)
    any_reports = any(reported for _, reported in results)
    return int(any_reports or (args.check and n_changed > 0))


# An attribute chain rooted at a Name. nodes has the Name node followed
# by the Attribute nodes, from the innermost to the outermost. names
# are the identifiers along the chain.
Chain = namedtuple('Chain', ('nodes', 'names'))

# A replacement of the text between two (line, column) positions
Edit = namedtuple('Edit', ('start', 'end', 'text'))


# A Qt attribute chain that does not resolve. path are the attribute
# names starting from the Qt module, index is the position in path of
# the first missing attribute, owner the object that misses it.
Missing = namedtuple('Missing', ('chain', 'qt_module', 'path', 'index',
                                 'owner', 'n_prefix'))


class QtFile:
    """A Python file that may use Qt."""

    def __init__(self, path, qt_modules, source=None):
        """Initialize instance.

        Parameters
        ----------
        path : Path
            The file to be processed.
        qt_modules : dict
            Keys are names of Qt modules (e.g., 'QtCore'), values the
            corresponding modules that are used to resolve names.
        source : str, optional
            The source code to process. If not given, it is read
            from `path`.
        """
        self.path = path
        self.qt_modules = qt_modules
        if source is None:
            source = path.read_bytes().decode('utf-8')
        self.source = source
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

    # ---------------------------- imports ----------------------------
    def find_imports(self, rewrite):
        """Collect names imported from Qt modules, rewriting imports."""
        for node in ast.walk(self.tree):
            if isinstance(node, ast.ImportFrom):
                self._handle_import_from(node, rewrite)
            elif isinstance(node, ast.Import):
                self._handle_import(node, rewrite)

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

    def fix_instance_enums(self, enum_index, rewrite):
        """Report (or fix) short enum members accessed via instances.

        These are, e.g., palette.Text or widget.tabBar().SelectLeftTab:
        the type of the object is unknown, so any attribute named like
        a short member of a Qt enumeration is suspicious. Attributes of
        Qt module names (handled by fix_chains), of self and cls (see
        fix_self_enums), of names imported or of classes defined in
        this file, as well as attributes this file assigns, are skipped.

        Parameters
        ----------
        enum_index : dict
            Keys are the names of members of Qt enumerations, values
            sets of (Qt class name, enumeration name) for each Qt
            enumeration that has such a member.
        rewrite : bool
            If True, '.EnumName' is inserted before the member if all
            candidates have the same enumeration name. Otherwise, the
            member is only reported.
        """
        skipped_roots = {*self.roots, 'self', 'cls',
                         *self._imported_names(), *self._class_names()}
        assigned = {
            node.attr for node in ast.walk(self.tree)
            if isinstance(node, ast.Attribute)
            and isinstance(node.ctx, ast.Store)
            }
        for node in ast.walk(self.tree):
            if (not isinstance(node, ast.Attribute)
                    or node.attr not in enum_index
                    or node.attr in assigned
                    or isinstance(node.ctx, ast.Store)):
                continue
            candidates = enum_index[node.attr]
            enum_names = sorted({enum for _, enum in candidates})
            _, removed_enum, replacement = REMOVED_MEMBERS.get(
                node.attr, Removed(None, None, None)
                )
            scoped = (isinstance(node.value, ast.Attribute)
                      and node.value.attr in enum_names)
            if ((scoped and not replacement)
                    or self._root_name(node.value) in skipped_roots):
                continue
            if replacement:
                enum_names = [removed_enum]
            elif len(enum_names) > 1:
                enum_names = self._guess_enum(node.value, candidates)
            if rewrite and len(enum_names) == 1:
                self._scope_member(node, enum_names[0], scoped, replacement)
                continue
            obj = ast.get_source_segment(self.source, node.value)
            owners = ', '.join(sorted(f'{cls}.{enum}.{node.attr}'
                                      for cls, enum in candidates))
            if replacement:
                hint = f'Removed in Qt 6: use {removed_enum}.{replacement}'
            elif len(enum_names) == 1:
                hint = f'Likely {enum_names[0]}'
            else:
                hint = 'Ambiguous: please fix manually'
            self.report(node, f'{obj}.{node.attr} may be a short enum '
                              f'member. Candidates: {owners}. {hint}.')

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

    def qt_chains(self):
        """Return the chains rooted at Qt names, in source order."""
        chains = [c for c in self._iter_chains() if c.names[0] in self.roots]
        return sorted(chains, key=lambda c: (c.nodes[0].lineno,
                                             c.nodes[0].col_offset))

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

    def report(self, node, msg):
        """Store a message about `node`."""
        self.reports.append(f'{self.path}:{node.lineno}: {msg}')

    def report_long_lines(self, new_source):
        """Report lines of `new_source` that became too long."""
        old_lines = set(self.source.splitlines())
        for line_no, line in enumerate(new_source.splitlines(), start=1):
            if len(line) > MAX_LINE_LENGTH and line not in old_lines:
                self.reports.append(f'{self.path}:{line_no}: line longer '
                                    f'than {MAX_LINE_LENGTH} characters')

    def resolve(self, chain):
        """Return the object that a Qt `chain` refers to, or MISSING."""
        qt_module, prefix = self.roots[chain.names[0]]
        obj = self.qt_modules[qt_module]
        for attr in (*prefix, *chain.names[1:]):
            obj = getattr(obj, attr, MISSING)
            if obj is MISSING:
                break
        return obj

    def _class_names(self):
        """Yield the names of the classes defined in this file."""
        for node in ast.walk(self.tree):
            if isinstance(node, ast.ClassDef):
                yield node.name

    def _end(self, node):
        """Return the (line, column) where `node` ends."""
        return (node.end_lineno,
                self._to_char_col(node.end_lineno, node.end_col_offset))

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

    @staticmethod
    def _guess_enum(node, candidates):
        """Return the enumeration names whose class matches `node`.

        The match is by name: e.g., palette or self.palette() match
        QPalette, tab_bar or tabBar() match QTabBar.

        Parameters
        ----------
        node : ast.AST
            The expression whose attribute is a short enum member.
        candidates : set
            (Qt class name, enumeration name) for all the enumerations
            with that member.

        Returns
        -------
        enum_names : list
            The sorted names of the enumerations whose class name is
            similar to the name of `node`. All candidate enumerations
            if none matches.
        """
        if isinstance(node, ast.Call):
            node = node.func
        if isinstance(node, ast.Attribute):
            name = node.attr
        elif isinstance(node, ast.Name):
            name = node.id
        else:
            name = ''
        name = name.replace('_', '').lower()

        def _matches(cls_name):
            cls_name = cls_name[1:].lower()  # Drop the leading 'Q'
            return len(name) >= MIN_GUESS_LENGTH and (
                cls_name.startswith(name) or name.startswith(cls_name)
                )
        matching = {enum for cls, enum in candidates if _matches(cls)}
        return sorted(matching or {enum for _, enum in candidates})

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

    def _imported_names(self):
        """Yield the names bound by import statements in this file."""
        for node in ast.walk(self.tree):
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                for alias in node.names:
                    yield (alias.asname
                           or alias.name.partition('.')[0])

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

    @staticmethod
    def _root_name(node):
        """Return the identifier at the root of an attribute chain."""
        while isinstance(node, ast.Attribute):
            node = node.value
        return node.id if isinstance(node, ast.Name) else None

    def _scope_member(self, node, enum_name, scoped, replacement):
        """Make the member that `node` accesses fully scoped.

        Parameters
        ----------
        node : ast.Attribute
            The node accessing a short enumeration member.
        enum_name : str
            The name of the enumeration to insert before the member.
        scoped : bool
            Whether `node` already accesses the member via `enum_name`.
        replacement : str or None
            The new name of the member, if it was removed in Qt 6.
        """
        if not replacement:
            value_end = self._end(node.value)
            self.edits.append(Edit(value_end, value_end, f'.{enum_name}'))
            return
        end = self._end(node)
        start = (end[0], end[1] - len(node.attr))
        new_text = replacement if scoped else f'{enum_name}.{replacement}'
        self.edits.append(Edit(start, end, new_text))

    def _start(self, node):
        """Return the (line, column) where `node` starts."""
        return node.lineno, self._to_char_col(node.lineno, node.col_offset)

    # ------------------------- node positions ------------------------
    def _to_char_col(self, line_no, byte_col):
        """Return the character column of a UTF-8 byte offset."""
        line = self.lines[line_no - 1]
        return len(line.encode('utf-8')[:byte_col].decode('utf-8'))


def _build_enum_index(qt_modules):
    """Return the members of the enumerations of Qt classes.

    Parameters
    ----------
    qt_modules : dict
        Keys are names of Qt modules, values the modules.

    Returns
    -------
    enum_index : dict
        Keys are names of members of Qt enumerations (only those that
        start with an uppercase letter), values sets of (Qt class name,
        enumeration name) for each Qt enumeration that has the member.
    """
    enum_index = {
        member: {(cls_name, enum_name)}
        for member, (cls_name, enum_name, _) in REMOVED_MEMBERS.items()
        }
    # Names such as QIcon are also members (of QMetaType.Type) but
    # are far more likely to refer to the classes. Skip them.
    module_names = {name for module in qt_modules.values()
                    for name in dir(module)}
    for module in qt_modules.values():
        for cls_name in dir(module):
            cls = getattr(module, cls_name)
            if not isinstance(cls, type) or issubclass(cls, enum.Enum):
                continue
            for name in dir(cls):
                if name.startswith('_'):
                    continue
                for member, enum_type in _iter_enum_members(cls, name):
                    if member[0].isupper() and member not in module_names:
                        # Use the class that defines the enumeration
                        # rather than those inheriting it (e.g., all
                        # QIODevice subclasses have QIODeviceBase's
                        # OpenModeFlag).
                        qualname = enum_type.__qualname__
                        owner, _, enum_name = qualname.rpartition('.')
                        enum_index.setdefault(member, set()).add(
                            (owner or cls_name, enum_name)
                            )
    return enum_index


def _compare(path, rev, qt_modules):
    """Return messages for Qt chains of `path` that changed since `rev`.

    The Qt attribute chains of the file at `rev` and of the current
    file are resolved with the installed binding, and the sequences
    of objects compared. Objects that only appear in one version are
    reported, unless they just moved within the file (e.g., into a
    local variable). Chains that changed but resolve to the same
    object (e.g., Qt.AlignLeft -> Qt.AlignmentFlag.AlignLeft) are not.
    """
    old_source = _read_at_revision(path, rev)
    if old_source is None:
        print(f'{path}: not in {rev}, skipped', file=sys.stderr)
        return []
    removed, added, examples = _find_unmatched(
        _get_qt_chains(path, qt_modules, old_source),
        _get_qt_chains(path, qt_modules),
        )
    messages = []
    for key in removed - added:
        chain = examples[key]
        messages.append(f'{path}: {_dotted(chain)} (line '
                        f'{chain.nodes[0].lineno} in {rev}) has no '
                        'equivalent any longer')
    for key in added - removed:
        chain = examples[key]
        messages.append(f'{path}:{chain.nodes[0].lineno}: '
                        f'{_dotted(chain)} is new since {rev}')
    return messages


def _dotted(chain):
    """Return the dotted name of `chain`."""
    return '.'.join(chain.names)


def _find_unmatched(old, new):
    """Return the keys of chains in only one of `old` and `new`.

    Parameters
    ----------
    old, new : tuple
        The chains and their keys, as returned by _get_qt_chains.

    Returns
    -------
    removed, added : Counter
        The keys of chains only in `old` and only in `new`.
    examples : dict
        One chain for each key, from `old` for removed keys.
    """
    matcher = difflib.SequenceMatcher(a=old[1], b=new[1], autojunk=False)
    blocks = matcher.get_matching_blocks()
    matched = (
        {b.a + i for b in blocks for i in range(b.size)},
        {b.b + i for b in blocks for i in range(b.size)},
        )
    counters = Counter(), Counter()
    examples = {}
    for (chains, keys), counter, matched_ids in zip((old, new), counters,
                                                    matched):
        for i, (chain, key) in enumerate(zip(chains, keys)):
            if i not in matched_ids:
                counter[key] += 1
                examples.setdefault(key, chain)
    return (*counters, examples)


def _get_qt_chains(path, qt_modules, source=None):
    """Return the Qt chains in `path` and keys of their objects."""
    qt_file = QtFile(path, qt_modules, source)
    qt_file.find_imports(rewrite=False)
    chains = qt_file.qt_chains()
    return chains, [_object_key(qt_file, chain) for chain in chains]


def _import_qt_modules(package):
    """Return a dict of the Qt modules in `package`."""
    modules = {}
    for name in QT_MODULES:
        try:
            modules[name] = importlib.import_module(f'{package}.{name}')
        except ImportError as exc:
            sys.exit(f'Cannot import {package}.{name}: {exc}')
    return modules


def _iter_enum_members(cls, name):
    """Yield (member name, enumeration) for attribute `name` of `cls`.

    In PyQt6, enumerations are Enum subclasses: if attribute `name` is
    one of them, yield all its members. In PyQt5, enumerations are sip
    enumtype instances and members are attributes of `cls`: if `name`
    is one of these members, yield only that.
    """
    value = getattr(cls, name, None)
    if isinstance(value, type) and issubclass(value, enum.Enum):
        yield from ((member, value) for member in value.__members__)
    elif type(type(value)).__name__ == SIP_ENUM_TYPE:
        yield name, type(value)


def _iter_files(paths):
    """Yield all Python files in `paths`, except those of qt_compat."""
    compat_folder = COMPAT.replace('.', '/')
    for path_str in paths:
        path = Path(path_str)
        files = sorted(path.rglob('*.py')) if path.is_dir() else (path,)
        yield from (f for f in files
                    if compat_folder not in f.resolve().as_posix())


def _object_key(qt_file, chain):
    """Return a hashable key for the object that `chain` refers to.

    Two chains have the same key if they refer to the same object
    with the installed binding (e.g., Qt.AlignLeft and
    Qt.AlignmentFlag.AlignLeft). Enumeration members are compared
    by type and value, as PyQt5 may give equal but distinct objects.
    qApp and QApplication.instance (that replaces it) cannot be
    resolved to the same object without a running QApplication, and
    have the same special key.
    """
    if QAPP in chain.names or chain.names[-2:] == ['QApplication',
                                                    'instance']:
        return (QAPP,)
    obj = qt_file.resolve(chain)
    if obj is MISSING:
        return (MISSING, _dotted(chain))
    if isinstance(obj, int):
        return (int, type(obj).__module__, type(obj).__qualname__,
                int(obj))
    try:
        hash(obj)
    except TypeError:
        return (id, id(obj))
    return (object, obj)


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
    mode.add_argument('--compare', metavar='REV',
                      help=('Do not write files. Check that the Qt '
                            'attribute chains that changed since git '
                            'revision REV resolve to the same objects '
                            'with the installed binding.'))
    parser.add_argument('--strings', action='store_true',
                        help=('Also report strings (e.g., docstrings) '
                              'that mention PyQt5, PyQt6, or Qt5.'))
    parser.add_argument('--fix-instance-enums', action='store_true',
                        help=('Also fix (rather than only report) short '
                              'enum members accessed via instances, e.g., '
                              'palette.Text. Review the changes!'))
    parser.add_argument('paths', nargs='+', metavar='PATH')
    args = parser.parse_args(argv)
    if (args.verify or args.compare) and args.fix_instance_enums:
        parser.error('--fix-instance-enums cannot be used with --verify '
                     'or --compare')
    return args


def _process(path, args, qt_modules, enum_index):
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
    qt_file.fix_instance_enums(enum_index, args.fix_instance_enums)
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


def _read_at_revision(path, rev):
    """Return the source of `path` at git revision `rev`, or None."""
    result = subprocess.run(  # noqa: S603
        ('git', 'show', f'{rev}:./{path.name}'),  # noqa: S607
        cwd=path.parent,
        capture_output=True,
        check=False,
        )
    return result.stdout.decode('utf-8') if not result.returncode else None


if __name__ == '__main__':
    sys.exit(main())
