"""Package classes of viperleed.gui.leedsim.

Contains Qt-independent classes that are used by widgets to display the
real space lattices and the LEED pattern. Used to be a single module
classes.py.

Modules
-------
beams_old
    Older implementation of classes for determining equivalence between
    LEED spots. Only supports a single structure at normal incidence.
equivalent_beams
    Classes for determining equivalence between LEED spots, including
    both multiple structures and non-normal beam incidence.
leedparameters
    Containers of information needed to construct a LEED pattern from
    one or more structural domains.
leedparser
    A ConfigParser for reading/writing the input files for the LEED
    Pattern Simulator.
leedpattern
    A full LEED pattern, with multiple structural and symmetry-induced
    domains. Potentially at non-normal incidence.
leedsubpattern
    A portion of a LEED pattern for one (structural, symmetry) domain.
    Used for plotting.
oldleedpatterns
    Old version of a LEED pattern (and its sub-patterns). Supports only
    one structural domain and normal incidence.
realspace
    Container of information of a real-space structure.
structdomains
    Container of information on multiple coexisting structural domains.
symdomains
    Containers of information on bulk-symmetry-induced domains.
woods_old
    Class for interpreting a Wood's notation.

Packages
--------
woods
    Class for interpreting a Wood's notation.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-03-13'
__license__ = 'GPLv3+'
