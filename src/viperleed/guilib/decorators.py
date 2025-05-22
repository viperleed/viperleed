"""Module decorators of guilib.

======================================
  ViPErLEED Graphical User Interface
======================================

Author: Michele Riva
Created: 2020-01-29

This module provides Qt-independent decorators based on the
wrapt Python module.
"""
import inspect
import cProfile
import pstats
from time import perf_counter as timer

import wrapt
try:
    from line_profiler import LineProfiler
except ImportError:
    pass


def ensure_decorates_class(superclass=object):
    """Make sure a @wrapt.decorator decrates a class.

    This decorator is meant to be applied to a @wrapt.decorator.
    Raises exceptions if the @wrapt.decorator on which it acts is
    not applied to a subclass of superclass.

    Parameters
    ----------
    superclass: class
                default: object

    Returns
    -------
    decorator : callable
        The original @wrapt.decorator if the @wrapt.decorator
        is applied to a subclass of superclass.
    Raises
    ------
    RuntimeError
        If the decorator was not applied to a subclass of superclass.
    """

    @wrapt.decorator
    def _wrapper(decorator, _, args, kwargs):
        wrapped = args[0]
        if not inspect.isclass(superclass):
            raise TypeError("@ensure_decorates_class() arg 1 must be a class")

        if not inspect.isclass(wrapped):
            raise RuntimeError(f"Decorator @{decorator.__name__} must be "
                               "applied to a class, not to a "
                               f"<{wrapped.__class__.__name__}>")

        if not issubclass(wrapped, superclass):
            raise RuntimeError(f"Decorator @{decorator.__name__} must be "
                               f"applied to a {superclass.__name__} subclass, "
                               f"not to a {wrapped.__name__} class.")

        return decorator(*args, **kwargs)
    return _wrapper


# The following decorators are useful for development, but should be removed
# otherwise when normally running


def profile_calls(sort_args=('cumulative',), print_args=(10,)):
    """Run cProfile on the decorated function.

    Notice that cProfile exclusively profiles function calls.
    """
    profiler = cProfile.Profile()

    def decorator(func):
        def inner(*args, **kwargs):
            # Disable pylint warning as there is no other sensible
            # way to actually profile the function.
            # pylint: disable=too-many-try-statements
            result = None
            try:
                profiler.enable()
                result = func(*args, **kwargs)
                profiler.disable()
            finally:
                stats = pstats.Stats(profiler)
                print("################",
                      f"Profiling funtion {func.__name__}",
                      "################")
                stats.sort_stats(*sort_args).print_stats(*print_args)
            return result
        return inner
    return decorator


def profile_lines(func):
    """Profile execution time of each line of the decorated function."""
    try:
        LineProfiler
    except NameError:
        raise RuntimeError("No line_profiler module found.") from None

    def _wrapper(*args, **kwargs):
        """Execute and profile function."""
        profiler = LineProfiler()
        profiled_func = profiler(func)
        try:
            result = profiled_func(*args, **kwargs)
        finally:
            profiler.print_stats()
        return result
    _wrapper.__name__ = func.__name__
    return _wrapper


def exec_time(func):
    """Measure execution time of the wrapped function."""
    def _wrapper(*args, **kwargs):
        """Execute function and print the execution time."""
        start_time = timer()
        result = func(*args, **kwargs)
        elapsed = timer() - start_time
        if elapsed > 1:
            mult = 1
            unit = ''
        elif elapsed > 1e-3:
            mult = 1e3
            unit = 'm'
        elif elapsed > 1e-6:
            mult = 1e6
            unit = 'u'
        else:
            mult = 1e9
            unit = 'n'
        fname = str(func).replace('<','').split(' at 0x')[0]
        print(f"Execution time of {fname}: {elapsed*mult:.1f} {unit}s")
        return result
    return _wrapper


def print_call(func):
    """Print the class and the name of the function that is executed."""
    def _wrapper(*args, **kwargs):
        fname = str(func).replace('<','').split(' at 0x')[0]
        print("Called", fname)
        return func(*args, **kwargs)

    return _wrapper
