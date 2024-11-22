"""
======================================
  ViPErLEED Graphical User Interface
======================================
 *** module guilib.leedsim.classes ***

Contains Qt-independent classes that are used by widgets to display the real
space lattices and the LEED pattern. Used to be a single module classes.py

Author: Michele Riva
Created: 2021-03-13
"""

from viperleed.guilib.leedsim.classes.woods import Woods
from viperleed.guilib.leedsim.classes.realspace import RealSpace
from viperleed.guilib.leedsim.classes.beams import LEEDEquivalentBeams
from viperleed.guilib.leedsim.classes.symdomains import LEEDSymmetryDomains
from viperleed.guilib.leedsim.classes.leedpattern import (LEEDStructuralDomains,
                                                          # LEEDPattern,
                                                          # LEEDsubpattern
                                                          )
from viperleed.guilib.leedsim.classes.oldleedpatterns import (LEEDPattern,  # These will be replaced once the new LEEDPattern/LEEDsubpattern work as they should
                                                              LEEDsubpattern)


