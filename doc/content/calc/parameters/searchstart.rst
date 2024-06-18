.. _searchstart:

SEARCH_START
============

SEARCH_START defines the configuration of the search parameters at the
beginning of the search.

See :cite:t:`kottckeNewApproachAutomated1997` for an explanation of the
search algorithm used by TensErLEED.

**Default:** crandom

**Allowed values:** random, centered, control, crandom

**Syntax:**

::

   SEARCH_START = control

-  ``random``: initializes all search parameters to a random index
-  ``centered``: sets all search parameters to the index that is closest to
   zero displacement, i.e. represents no change from the input structure
-  ``control``: reads starting configuration from control.chem
-  ``crandom``: initializes *one* member of the population centered, all
   others random

If ``control`` is selected but no control.chem file is found, SEARCH_START
will default to ``crandom`` initialization.
