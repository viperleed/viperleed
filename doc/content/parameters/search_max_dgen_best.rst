.. _search_max_dgen_best:

SEARCH_MAX_DGEN_BEST
====================

SEARCH_MAX_DGEN_BEST defines a convergence criterion for the search.
If SEARCH_MAX_DGEN_BEST generations pass without any changes to the
*best* configuration (full convergence), the search will be stopped.

**Default:** not active, search is limited only by
:ref:`SEARCH_MAX_GEN<SEARCHGENMAX>`

**Allowed values:** Positive integer

**Syntax:**

::

   SEARCH_MAX_DGEN_BEST = 1000
