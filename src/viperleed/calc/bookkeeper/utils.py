"""Module utils of viperleed.calc.bookkeeper.

Collects functions used in multiple spots in the bookkeeper package.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-10-14'
__license__ = 'GPLv3+'

import filecmp
from functools import wraps
from operator import attrgetter
import shutil

from .log import LOGGER


def discard_files(*file_paths):
    """Delete files at `file_paths`. Log if they can't be deleted."""
    for file in file_paths:
        if not file.exists():
            continue
        if file.is_file():
            try:
                file.unlink()
            except OSError:
                LOGGER.error(f'Failed to discard file {file.name}.')
            continue  # Marked as uncovered by pytest, but it is
        assert file.is_dir()  # Should be a directory
        try:
            shutil.rmtree(file)
        except OSError:
            LOGGER.error(f'Failed to discard directory {file.name}.')


def file_contents_identical(file_one, file_two):
    """Return whether two files have the same contents."""
    try:
        return filecmp.cmp(file_one, file_two, shallow=False)
    except FileNotFoundError:
        return False


def make_property(attr, needs_update=False, updater='update_from_cwd'):
    """Return a property getter that looks up `attr`."""
    getter, property_name = _get_attr_or_dict_item(attr)
    if needs_update:
        decorator = needs_update_for_attr(attr,
                                          attr_name=property_name,
                                          updater=updater)
        getter = decorator(getter)
    return property(getter)


def needs_update_for_attr(attr, attr_name=None, updater='update_from_cwd'):
    """Return a decorator that complains if `attr` is None.

    Parameters
    ----------
    attr : str
        The attribute of self to be looked up. If the value is None,
        the decorated function raises AttributeError suggesting to
        call update_from_cwd beforehand. `attr` may also have the form
        'dict_name[dict_item]'. In this case, self.dict_name[dict_item]
        is looked up instead.
    attr_name : str, optional
        The name of the decorated function. Automatically fetched from
        the decorated function if not given. Default is None.
    updater : str, optional
        Which method should appear as the "updater" method o be called
        to make `attr_name` available. Default is 'update_from_cwd'.

    Returns
    -------
    _decorator : callable
        A decorator to apply to a function that raises AttributeError
        before calling the function if `attr` is None.
    """
    _getattr, *_ = _get_attr_or_dict_item(attr)
    def _decorator(func):
        func_name = attr_name or func.__name__
        @wraps(func)
        def _wrapper(self, *args, **kwargs):
            try:
                value = _getattr(self)
            except KeyError:
                value = None
            if value is None:
                raise AttributeError(
                    f'{type(self).__name__} has no {func_name} yet. '
                    f'Call {updater}() beforehand.'
                    )
            return func(self, *args, **kwargs)
        return _wrapper
    return _decorator


def _get_attr_or_dict_item(attr):
    """Return an attribute getter for 'attr_name' or 'dict_name[item_name]'.

    Parameters
    ----------
    attr : str
        The attribute of self to be looked up. If the value is None,
        the decorated function raises AttributeError suggesting to
        call update_from_cwd beforehand. `attr` may also have the form
        'dict_name[dict_item]'. In this case, self.dict_name[dict_item]
        is looked up instead.

    Returns
    -------
    getter : callable
        A function that returns `attr` from its first argument.
    cleaned_attr : str
        'attr_name' or 'item_name', depending on the form of `attr`.
    """
    if '[' in attr:  # pylint: disable=magic-value-comparison
        dict_, attr = attr.split('[')
        attr = attr.replace(']', '')
        def _getattr(self):
            container = attrgetter(dict_)(self)
            return container[attr]
    else:
        _getattr = attrgetter(attr)
    return _getattr, attr
