.. _searchgenmaxdelta:

SEARCH_MAX_DGEN
===============

SEARCH_MAX_DGEN defines a convergence criterion for the search. If
SEARCH_MAX_DGEN generations pass without any changes for any structure,
the search will be stopped.

.. note::
   This is a very naive convergence parameter, and others (e.g.,
   :ref:`SEARCH_MAX_DGEN_BEST<search_max_dgen_best>` ) might be
   more suitable for most applications.

**Default:** not active, search is limited only by
:ref:`SEARCH_MAX_GEN<SEARCHGENMAX>`

**Allowed values:** Positive integer

**Syntax:**

::

   SEARCH_MAX_DGEN = 2000
