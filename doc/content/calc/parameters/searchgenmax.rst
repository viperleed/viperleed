.. _searchgenmax:

SEARCH_MAX_GEN
==============

SEARCH_MAX_GEN defines the maximum number of generations that the search will
go through before it stops.

See :cite:t:`kottckeNewApproachAutomated1997` for an explanation of the search
algorithm used by TensErLEED.

**Default:** 1e5

**Allowed values:** Positive integer

**Syntax:**

::

   SEARCH_MAX_GEN = 5e6

When automatic refinement of :ref:`RMUT` is active, SEARCH_MAX_GEN will
nevertheless count *all* steps and stop after the specified total number,
instead of resetting whenever the Gaussian width changes (see
:ref:`SEARCH_CONVERGENCE`). SEARCH_MAX_GEN is therefore the
parameter most suited for setting an upper limit to the search
time, if convergence is poor.
