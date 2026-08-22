.. include:: /substitutions.rst

.. _rfactortype:

R_FACTOR_TYPE
=============

R_FACTOR_TYPE determines what definition of the |R factor| is used in
|R-factor| calculations, including during the search. For details, see
|R-factor| :ref:`calculation<r-factor_calculation>`.

Note that var(R) in :ref:`error calculations<error_calculation>` is only supported
for Pendry and Smooth |R-factor|, and Zanazzi-Jona is not supported in the
TensErLEED :ref:`structure search<sec_search>`.

**Default:** pendry

**Allowed values:** Pendry (1), R2 (2), Zanazzi-Jona (zj, 3), Smooth (4)

**Syntax:**

::

   R_FACTOR_TYPE = smooth



.. versionchanged:: 0.15.0
   R_FACTOR_TYPE now accepts definition by keywords. Earlier versions only
   accept integer numbers 1 (Pendry) and 2 (R2). |R-factor| type Smooth is also
   only supported since version 0.15.0.
