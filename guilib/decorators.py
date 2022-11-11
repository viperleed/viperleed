"""
======================================
  ViPErLEED Graphical User Interface
======================================
   *** Module guilib.decorators ***

Author: Michele Riva
Created: 2020-01-29

This module provides Qt-independent decorators based on the wrapt Python module
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
    """
    Decorator to be applied to a @wrapt.decorator. Raises exceptions if the
    @wrapt.decorator on which it acts is not applied to a subclass of
    superclass.

    Parameters
    ----------
    superclass: class
                default: object

    Returns
    -------
    The original @wrapt.decorator if the @wrapt.decorator is applied to a
    subclass of superclass. Raises RuntimeError otherwise.
    """

    @wrapt.decorator
    def _wrapper(decorator, always_none, args, kwargs):
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


def profile_calls(sort_args=['cumulative'], print_args=[10]):
    """
    Runs a profiler on the wrapped function, using cProfile. This profiles
    exclusively function calls.
    """
    profiler = cProfile.Profile()

    def decorator(func):
        def inner(*args, **kwargs):
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
        profiler = LineProfiler()
        profiled_func = profiler(func)
        try:
            result = profiled_func(*args, **kwargs)
        finally:
            profiler.print_stats()
        return result
    return _wrapper


def exec_time(func):
    """
    Measures execution time of the wrapped function
    """
    def _wrapper(*args, **kwargs):
        t0 = timer()
        result = func(*args, **kwargs)
        dt = timer()-t0
        if dt > 1:
            mult = 1
            unit = ''
        elif dt > 1e-3:
            mult = 1e3
            unit = 'm'
        elif dt > 1e-6:
            mult = 1e6
            unit = 'u'
        else:
            mult = 1e9
            unit = 'n'
        fname = str(func).replace('<','').split(' at 0x')[0]
        print(f"Execution time of {fname}: {dt*mult:.1f} {unit}s")
        return result
    return _wrapper

