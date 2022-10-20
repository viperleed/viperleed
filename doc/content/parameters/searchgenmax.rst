.. _searchgenmax:

SEARCH_MAX_GEN
==============

SEARCH_MAX_GEN defines the maximum number of generations that the search will go through before it stops.

See Ref. `[1 <SEARCH_CONVERGENCE#ref1>`__] for an explanation of the search algorithm used by TensErLEED.

**Default:** 1e5

**Allowed values:** Positive integer

**Syntax:**

::

   SEARCH_MAX_GEN = 5e6

When automatic refinement of `RMUT </protected/surface/LEEDIV/PARAMETERS/RMUT>`__ is active, SEARCH_MAX_GEN will nevertheless count *all* steps and stop after the specified total number, instead of resetting whenever the Gaussian width changes (see :ref:`SEARCH_CONVERGENCE<SEARCH_CONVERGENCE>`). SEARCH_MAX_GEN is therefore the parameter most suited for setting an upper limit to the search time, if convergence is poor.

M. Kottcke and K. Heinz, *A New Approach to Automated Structure Optimization in LEED Intensity Analysis* `Surf. Sci. 376, 352 (1997) <http://dx.doi.org/10.1016/S0039-6028(96)01307-6>`__.
