.. include:: /substitutions.rst

.. _backend:

=======
BACKEND
=======

.. warning:: The use of the viperleed-jax backend is experimental.

BACKEND selects the backend used by viperleed.calc for calculations. Currently
only the :ref:sec_search: supports different backends (TensErLEED and
viperleed-jax).

**Syntax**:

::

   BACKEND search = viperleed-jax ! use the viperleed-jax backend for the search


Search Backends
---------------

The search backend can be selected with ``BACKEND search = ...`` to choose which
backend is used for the structure optimization. Note that this will also affect
parsing of the :ref:`displacements` file.
