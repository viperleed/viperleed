"""Module cli_base of viperleed.

Defines classes and functions that can (and should) be used to add new
command-line utilities to viperleed.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-03-26'
__license__ = 'GPLv3+'

import argparse
from collections import defaultdict
from functools import partial
import importlib
import inspect
import logging
from pathlib import Path
import pkgutil
import sys

from viperleed import GLOBALS
from viperleed.calc.lib.string_utils import parent_name


# Default name of the submodules of package utilities
CLI_MODULE_NAME = 'cli'


class NotAViPErLEEDCLIError(Exception):
    """Object is not a ViPErLEEDCLI, or module does not define one."""


def _to_float(value_str):
    """Raise ArgumentTypeError if `value_str` is not a float."""
    try:
        return float(value_str)
    except ValueError:
        raise argparse.ArgumentTypeError('Must be a floating-point '
                                         'value') from None


def float_in_zero_one(value_str):
    """Raise ArgumentTypeError if `value_str` is not in range (0, 1).

    This can be used as a type checker for parser.add_argument.

    Parameters
    ----------
    value_str : str
        The raw value supplied via the command-line.

    Raises
    ------
    ArgumentTypeError
        If `value_str` cannot be converted to float or if it is
        not in range (0, 1).
    """
    value = _to_float(value_str)
    if not 0 < value < 1:
        raise argparse.ArgumentTypeError('Must be between zero '
                                         'and one (excluded)')
    return value


def length_choices(*choices):
    """Return an action that forces `choices` as the only acceptable nargs."""

    if len(choices) > 2:  # pylint: disable=magic-value-comparison
        fmt_choices = ', '.join(str(c) for c in choices[:-1])
        fmt_choices += f', or {choices[-1]}'
    else:
        fmt_choices = ' or '.join(str(c) for c in choices)

    # pylint: disable-next=too-few-public-methods
    class _LengthChoices(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            n_values = len(values)
            if n_values not in choices:
                err_ = f'argument {self.dest!r} requires {fmt_choices} values'
                parser.error(err_)
            setattr(args, self.dest, values)

    return _LengthChoices


def positive_float(value_str):
    """Raise ArgumentTypeError if `value_str` is not strictly positive.

    This can be used as a type checker for parser.add_argument.

    Parameters
    ----------
    value_str : str
        The raw value supplied via the command-line.

    Raises
    ------
    ArgumentTypeError
        If `value_str` cannot be converted to float or if it is <= 0.
    """
    value = _to_float(value_str)
    if value <= 0:
        raise argparse.ArgumentTypeError('Must be strictly positive')
    return value


def required_length(n_min=None, n_max=None):
    """Return an action ensuring there are between n_min and n_max items."""
    # Code adapted from https://stackoverflow.com/questions/4194948

    # pylint: disable-next=too-few-public-methods
    class _RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            n_values = len(values)
            err_ = f'argument {self.dest!r} requires '
            if n_min is not None and n_values < n_min:
                err_ += f'at least {n_min} values'
                parser.error(err_)
            if n_max is not None and n_values > n_max:
                err_ += f'at most {n_max} values'
                parser.error(err_)
            setattr(args, self.dest, values)

    return _RequiredLength


minimum_length = partial(required_length, n_max=None)
maximum_length = partial(required_length, n_min=None)


def strip_cli_module(module_name):
    """Return `module_name` without a terminating '.cli' (if it has one)."""
    return parent_name(module_name, remove=CLI_MODULE_NAME)


class ViPErLEEDCLI:
    """Base class for any viperleed command-line interface.

    How to create a new viperleed CLI/utility
    -----------------------------------------
    Subclass `ViPErLEEDCLI`, optionally passing a name.
        The subclass should be defined in a script-like module (small
        single-module utilities) or in a cli.py module of a utility
        package. Subclasses whose name starts with an underscore are
        considered to be base implementations to be further extended
        and will not be collected as valid children. Normally, each
        utility module should define only a single ViPErLEEDCLI
        subclass. For info, see
        >>> help(ViPErLEEDCLI.__init_subclass__)
    If needed, override the __call__ method.
        This should be the main functionality of the CLI tool. For
        detailed info, see
        >>> help(ViPErLEEDCLI.__call__)
    Extend the `add_parser_arguments` method if the utility has special
        arguments. For detailed info, see
        >>> help(ViPErLEEDCLI.add_parser_arguments)
    Checking of command-line arguments should be done in
        `parse_cli_args`, and, only in exceptional cases
        in `__call__`. For detailed info, see
        >>> help(ViPErLEEDCLI.parse_cli_args)
    A utility can have children sub-utilities. If this is the case,
        call `register_child`, typically in `__init__`. Notice that
        the order in which children are registered is the one with
        which they will later appear. To automatically filter out
        invalid modules/packages (i.e., those that do not define a
        subclass of `ViPErLEEDCLI`) use `register_valid_children.`
        See also the `ViPErLEEDCLIWithAutoChildren` subclass of this
        module that automatically finds all valid children.
    Children can be given aliases too (with `add_child_aliases`).
        They will appear as alternative command names.
    """

    @classmethod
    def __init_subclass__(cls, *,
                          cli_name=None,
                          description=None,
                          help_=None,
                          **kwargs):
        """Create a ViPErLEED command-line interface subclass.

        Parameters
        ----------
        cli_name : str, optional
            Typically the __name__ of the module that implements the
            children utility. If not given, cls.__module__ (stripped
            of the last .cli entry) is used instead. Default is None.
        description : str, optional
            Description text that will be added to the CLI.
            Default is None.
        help_ : str, optional
            Help text to be displayed when this utility is used as a
            child of another one. If not given, a default help text
            is generated from `cli_name`. Default is None.
        **kwargs : dict
            Other optional arguments, passed on to type.

        Returns
        -------
        None.
        """
        super().__init_subclass__(**kwargs)
        # For each class, store a {module_name: cls} register of all
        # the ViPErLEEDCLI classes that would like to be a child
        # of this utility. module_name is the cli_name of the module
        # in which the CLI utility is defined. The cli_name excludes,
        # however a terminating '.cli', as sub-utilities should
        # normally define their command-line interface either in
        # the main module itself (for a single-module utility) or
        # in a cli.py module.
        # Notice that each new subclass MUST have its own dedicated
        # _children dictionary to avoid infinite recursion when
        # sub-classing twice the same subclass of ViPErLEEDCLI.
        cls._children = {}
        cls._cli_name = cli_name or strip_cli_module(cls.__module__)
        cls._description = description
        cls._help = help_
        cls._aliases = defaultdict(set)  # {child: {aliases}}

    def __call__(self, args=None):
        """Call the main functionality of this ViPErLEEDCLI.

        This method must be overridden in concrete subclasses that do
        not have any children, i.e., those that perform any action.
        The overridden method must execute the primary functionality
        of the CLI utility. On the other hand, CLIs with children do
        not need (and typically must not) modify this method to ensure
        proper resolution of the children's functions.

        Parameters
        ----------
        args : Sequence or argparse.Namespace or None, optional
            Command-line arguments for this ViPErLEEDCLI. When
            overriding __call__, use self.parse_cli_args to properly
            resolve both global and specific arguments. If not given
            or None, sys.argv is used instead, as implemented in
            argparse.ArgumentParser. Default is None.

        Returns
        -------
        exit_code : int
            Commonly the following meaning for exit codes is assumed:
            0   -> success
            1   -> clean exit (e.g., via KeyboardInterrupt)
            > 1 -> error
        """
        args = self.parse_cli_args(args)
        try:
            command = args.func
        except AttributeError:  # Called without arguments
            command = None
        # About disable: we really want to compare command
        # with the callable object, NOT with the result.
        # pylint: disable-next=comparison-with-callable
        if command == self.__call__:  # Avoid infinite recursion
            command = None
        if not command:
            self.parser.parse_args(['--help'])
            return 0  # pragma: no cover  # Help does sys.exit
        return command(args)

    @property
    def children(self):
        """Return all the children of this utility."""
        return self._children

    @property
    def cli_name(self):
        """Return the name of this CLI utility."""
        return self._cli_name

    @property
    def help(self):
        """Return the short help text for this CLI. Used by parents."""
        return self._help

    @property
    def parser(self):
        """Return an ArgumentParser for this CLI with all arguments added."""
        parser = argparse.ArgumentParser(prog=self.cli_name,
                                         description=self._description)
        self.add_parser_arguments(parser)
        return parser

    @classmethod
    def find_sub_utilities(cls):
        """Look up the file system for children of this CLI/utility.

        Returns
        -------
        sub_utilities : tuple
            Fully qualified names of the packages or modules that are
            (potentially) sub-utilities of this one. Packages appear
            before modules. Only those packages that have a cli.py
            module are returned. No check is performed on whether the
            modules (or the cli.py modules of packages) do subclass
            `ViPErLEEDCLI`.
        """
        module_name = cls.__module__
        if not module_name.endswith(f'.{CLI_MODULE_NAME}'):
            # No children for a module-only utility
            return tuple()
        base_path = Path(importlib.util.find_spec(module_name).origin).parent
        all_modules = [module
                       # str(Path) due to cpython/issues/88227
                       for module in pkgutil.iter_modules((str(base_path),))
                       if module.name not in {'__main__', CLI_MODULE_NAME}]
        # All the modules are potentially OK
        modules = (module.name
                   for module in all_modules if not module.ispkg)
        # The packages, however, must have a cli.py module
        packages = []
        for package in all_modules:
            if not package.ispkg:
                continue
            package_path = Path(package.module_finder.path) / package.name
            sub_modules = pkgutil.iter_modules((str(package_path),))
            try:
                cli_module = next(f'{package.name}.{m.name}'
                                  for m in sub_modules
                                  if m.name == CLI_MODULE_NAME)
            except StopIteration:
                continue  # No cli module
            packages.append(cli_module)
        base_name = strip_cli_module(module_name)
        return tuple(f'{base_name}.{module_name}'
                     for module_name in (*packages, *modules))

    @classmethod
    def get_logger(cls):
        """Return a logging.Logger for the correct module/package.

        This method should be used especially by base sub-classes
        of ViPErLEEDCLI that would like to access the same logger
        as the one returned by logging.getLogger(__name__), where
        __name__ is the module that implements the concrete subclass
        of the base class.

        Returns
        -------
        logger: logging.Logger
            A logger for the module in which `cls` is defined
            (except a trailing .cli for package utilities).
        """
        return logging.getLogger(strip_cli_module(cls.__module__))

    @classmethod
    def register_child(cls, child):
        """Register a sub-CLI of this one.

        Parameters
        ----------
        child : str or ViPErLEEDCLI or type
            If a string, it is interpreted as the name of a module that
            defines a ViPErLEEDCLI subclass. Otherwise an instance of a
            ViPErLEEDCLI subclass or a ViPErLEEDCLI subclass itself are
            also acceptable.

        Raises
        ------
        ValueError
            If `child` is a string but errors occur while trying to
            import it as a module.
        NotAViPErLEEDCLIError
            If `child` is a module name but the module does not define
            a concrete `ViPErLEEDCLI` subclass, i.e., one whose name
            does not start with an underscore.
        NotAViPErLEEDCLIError
            If `child` is neither a `ViPErLEEDCLI` subclass nor an
            instance of a `ViPErLEEDCLI` subclass.
        """
        if isinstance(child, str):  # Assume it's a module's name
            child = cls._child_class_from_module_name(child)
        if inspect.isclass(child) and issubclass(child, ViPErLEEDCLI):
            child_cls = child
        elif isinstance(child, ViPErLEEDCLI):
            child_cls = type(child)
        else:
            raise NotAViPErLEEDCLIError('child must be a ViPErLEEDCLI')
        child_name = strip_cli_module(child_cls.__module__)
        cls._children[child_name] = child_cls

    @staticmethod
    def _child_class_from_module_name(module_name):
        """Return a concrete ViPErLEEDCLI subclass from `module_name`."""
        try:
            module = importlib.import_module(module_name)
        except ImportError as exc:
            raise ValueError(f'Cannot register {module_name}.') from exc
        for name, member in inspect.getmembers(module):
            if not inspect.isclass(member):
                continue
            if member in (ViPErLEEDCLI, ViPErLEEDCLIWithAutoChildren):
                continue
            if name.startswith('_'):
                continue
            if issubclass(member, ViPErLEEDCLI):
                return member
        raise NotAViPErLEEDCLIError(f'No ViPErLEEDCLI in module {module_name}')

    @classmethod
    def register_valid_children(cls, children):
        """Register the `children` that subclass ViPErLEEDCLI."""
        for child in children:
            try:
                cls.register_child(child)
            except NotAViPErLEEDCLIError:
                pass

    @classmethod
    def run_as_script(cls):
        """Run this CLI as stand-alone script, reading argv.

        This class method is intended to be used in __main__.py
        files or in "if __name__ == '__main__':"-fenced blocks
        in scripts. It is a safe 'main' function that gracefully
        handles KeyboardInterrupt.

        Raises
        ------
        SystemExit
            With a suitable exit code. If unhandled, this will
            TERMINATE THE INTERPRETER.
        """
        cli = cls()
        try:
            exit_code = cli()
        except KeyboardInterrupt:
            print('Terminated by keyboard interrupt')
            exit_code = 1
        sys.exit(exit_code)

    def add_child_aliases(self, child_cli_name, *aliases):
        """Add `aliases` for `child_cli_name`."""
        self._aliases[child_cli_name].update(aliases)

    def add_global_arguments(self, parser):
        """Add arguments to `parser`. This method is called on all children."""
        parser.add_argument('--version',
                            help='print version number',
                            action='version',
                            version=GLOBALS['version_message'])

    def add_parser_arguments(self, parser):
        """Add CLI arguments for this utility to `parser`.

        This method must be **extended** in subclasses that want to
        introduce new arguments: call super().add_parser_arguments
        in the extended method.

        The base implementation does `add_global_arguments` to parser,
        and adds one sub-parser for each of the children utilities of
        this one. All the added parsers have global arguments as well
        as their own arguments (via their own `add_parser_arguments`).
        This ensures that children utilities are added recursively.

        Arguments
        ---------
        parser : ArgumentParser
            The parser to which new arguments should be added.

        Returns
        -------
        None.
        """
        self.add_global_arguments(parser)
        if not self.children:
            return
        subparsers = parser.add_subparsers(dest='command')
        for child_cls in self.children.values():
            child = child_cls()
            help_txt = child.help or f'call {child.cli_name}'
            aliases = sorted(self._aliases[child.cli_name])
            child_parser = subparsers.add_parser(child.cli_name,
                                                 help=help_txt,
                                                 aliases=aliases)
            child_parser.set_defaults(func=child.__call__)
            # Explicitly add all the child's arguments to the parser
            # created by subparser.add_parser. We would normally do
            # this via child.parser, but we have to do it explicitly
            # here in order to get all the correct sub-children already
            # at this stage. This is important to obtain a proper help
            # message when calling the child without arguments or with
            # the usual -h/--help
            child.add_parser_arguments(child_parser)

    @staticmethod
    def add_verbose_option(parser):
        """Add --verbose flag to `parser`."""
        parser.add_argument('-v', '--verbose',
                            help='increase output verbosity',
                            action='store_true')

    def parse_cli_args(self, args):
        """Return a version of `args` parsed with this CLI's parser.

        Subclasses that need to perform some specific checking of
        arguments can proceed in two ways:
        1. Arguments that do not depend on other arguments should
            be given a dedicated checker in add_parser_arguments
            using the functionality provided by argparse. Examples
            include:
            - parser.add_mutually_exclusive_group;
            - the `type` kwarg of parser.add_argument using a dedicated
                callable, like those implemented in this module (e.g.,
                float_in_zero_one, positive_float)
            - the `action` kwarg of parser.add_argument also using a
                dedicated action, like those implemented in this
                module (e.g., length_choices, required_length,
                minimum_length, maximum_length). Notice that actions
                should use parser.error(msg) to properly report error
                conditions rather than raising ArgumentTypeError.
        2. Arguments that require consistency checks with other parsed
            arguments should be checked in an **extended** version of
            parse_cli_args, or, in rare cases, in a reimplementation
            of __call__ (if some non-CLI-related processing is needed
            for some of the arguments before). These modified methods
            must use self.parser.error(msg) to properly report argument
            errors in the format of argparse rather than raising
            exceptions.

        Parameters
        ----------
        args : Sequence or argparse.Namespace or None
            The command-line arguments to be parsed using self.parser.
            If None, sys.argv is used, as per argparse.ArgumentParser
            implementation.

        Returns
        -------
        parsed_args : argparse.Namespace
        """
        if isinstance(args, argparse.Namespace):
            return args  # Parsed already by someone else
        return self.parser.parse_args(args)


class ViPErLEEDCLIWithAutoChildren(ViPErLEEDCLI):
    """A command-line interface that automatically adds all children."""

    def __init_subclass__(cls, *_, **kwargs):
        """Register all valid sub-utilities found."""
        super().__init_subclass__(**kwargs)
        cls.register_valid_children(cls.find_sub_utilities())
