"""Package widgets of viperleed.gui.leedsim.

Contains QWidgets used within the pattern simulator
plug-in of the ViPErLEED Graphical User Interface.

Modules
-------
bulkinput
    A widget for users to specify the periodicity and symmetry of
    the (common) bulk of the structure(s).
domainsblock
    A widget to toggle display of symmetry-induced domains and
    displaying the corresponding superlattice matrices.
editablematrix
    A square matrix of widgets that allow, for example, to specify
    a SUPERLATTICE matrix.
energyblock
    A widget for setting the current LEED energy.
hoverannot
    An arrow with a corresponding text box for labeling selected spots.
latticeinput
    A widget for users to specify the periodicity and symmetry of
    a 2D lattice.
leedcanvas
    A matplotlib canvas for displaying a LEED pattern.
matricespopup
    A popup displaying superlattice matrices.
realcanvas
    A matplotlib canvas for displaying a real-space lattice.
rotationblock
    A widget for setting the current view angle.
surfaceinput
    A widget for users to specify the periodicity and symmetry of
    a single structural domain.
togglebutton
    A button with a larger-than-normal size.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-06-01'
__license__ = 'GPLv3+'
