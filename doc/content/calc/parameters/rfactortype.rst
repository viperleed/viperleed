.. include:: /substitutions.rst

.. _rfactortype:

R_FACTOR_TYPE
=============

R_FACTOR_TYPE determines what definition of the |R factor| is used in
|R-factor| calculations, including during the search. For details, see
|R-factor| :ref:`calculation<r-factor_calculation>`.

**Default:** pendry

**Allowed values:** pendry (1), r2 (2), zj (3), smooth (4)

**Syntax:**

::

   R_FACTOR_TYPE = smooth   ! equivalent: R_FACTOR_TYPE = 4
   

.. warning::
    The "smooth" |R-factor| |RS| is currently only implemented for the
    :ref:`sec_search` using the viperleed-jax :ref:`BACKEND`. All segments
    except for the :ref:`sec_search` will default to |RP| when |RS| is picked,
    printing a warning.

.. todo::
    Update once |RS| is supported in old TensErLEED.

Changelog
---------

.. versionchanged:: 0.15.0
   Added support for |RS|
.. versionchanged:: 0.15.0
   Accepts string definitions in addition to integers
   
