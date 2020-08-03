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

import wrapt


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

