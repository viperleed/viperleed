.. _searchstart:

SEARCH_START
============

SEARCH_START defines the configuration of the search parameters at the beginning of the search.

See Ref. `[1 <SEARCH_CONVERGENCE#ref1>`__] for an explanation of the search algorithm used by TensErLEED.

**Default:** crandom

**Allowed values:** random, centered, control, crandom

**Syntax:**

::

   SEARCH_START = control

-  ``random``: initializes all search parameters to a random index
-  ``centered``: sets all search parameters to the index that is closest to zero displacement, i.e. represents no change from the input structure
-  ``control``: reads starting configuration from control.chem
-  ``crandom``: initializes *one* member of the population centered, all others random

If ``control`` is selected but no control.chem file is found, SEARCH_START will default to ``crandom`` initialization.

M. Kottcke and K. Heinz, *A New Approach to Automated Structure Optimization in LEED Intensity Analysis* `Surf. Sci. 376, 352 (1997) <http://dx.doi.org/10.1016/S0039-6028(96)01307-6>`__.
