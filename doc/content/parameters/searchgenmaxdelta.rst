.. _searchgenmaxdelta:

SEARCH_MAX_DGEN
===============

SEARCH_MAX_DGEN defines a convergence criterion for the search. If SEARCH_MAX_DGEN generations pass without any changes for any structure, the search will be stopped.

**Default:** not active, search is limited only by :ref:`SEARCH_MAX_GEN<SEARCHGENMAX>` 

**Allowed values:** Positive integer

**Syntax:**

::

   SEARCH_MAX_DGEN = 2000

**ToDo - move comment: This is only a very naive convergence parameter, others might be better.**
