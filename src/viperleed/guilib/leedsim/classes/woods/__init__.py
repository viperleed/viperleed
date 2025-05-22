"""Package woods of guilib.leedsim.classes.

======================================
  ViPErLEED Graphical User Interface
======================================
 *** package guilib.leedsim.classes.woods ***

Author: Michele Riva (@michele-riva)
Created: 2024-03-01

Defines the Woods class, used for conversion of a Wood's notation
string into a matrix and vice versa, as well as for formatting. It
also defines exceptions related to Wood's notation interpretation.

This package was refactored from the former woods.py module, and
exposes the same API.
"""

from ._woods import Woods
from .errors import MatrixIncommensurateError
from .errors import WoodsError
from .errors import WoodsInvalidForBasisError
from .errors import WoodsNotRepresentableError
from .errors import WoodsSyntaxError
