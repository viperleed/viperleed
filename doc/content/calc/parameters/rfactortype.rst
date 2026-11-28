.. include:: /substitutions.rst

.. _rfactortype:

R_FACTOR_TYPE
=============

R_FACTOR_TYPE determines what definition of the |R factor| is used in
|R-factor| calculations, including during the search. For details, see
|R-factor| :ref:`calculation<r-factor_calculation>`.

**Default:** pendry

**Allowed values:** One of the four |R-factor| types, given by name, by a
synonym, or by its integer code:

-  Pendry |R factor| |RP|: ``pendry``, ``pendry_r``, ``rp``, ``pend``,
   ``p``, ``r_pendry``, ``pendry-r-factor``, or ``1``
-  :math:`R_2`: ``r2``, ``r-squared``, ``r2-factor``, or ``2``
-  Zanazzi-Jona |R factor|: ``zj``, ``zanazzi-jona``, ``zannazi``,
   ``zannazi-jona``, ``z``, ``zj-factor``, or ``3``
-  Smooth |R factor| |RS|: ``smooth``, ``smoothed``, ``smooth-pendry``,
   ``schmid``, ``s``, ``rs``, ``r_smooth``, or ``4``

Names are not case sensitive, and underscores, hyphens, and ``^`` are
ignored (e.g., ``Zanazzi_Jona``, ``R^2``, and ``R²`` are also accepted).
The value must be a single word: use, e.g., ``r-squared`` rather than
``r squared``.

**Syntax:**

::

   R_FACTOR_TYPE = smooth

Which implementation calculates the |R factor|
----------------------------------------------

Two implementations of the |R factor| exist in ViPErLEED: a python-native
one and the one of TensErLEED. Which one is used depends on where the
|R factor| is calculated:

-  The |R-factor| calculation after a :ref:`reference calculation<ref-calc>`
   or a :ref:`superposition calculation<super_pos>`, for
   :ref:`error curves<error_calculation>`, and in the
   :ref:`full-dynamic optimization<fdoptimization>` uses the
   implementation selected by :ref:`R_FACTOR_LEGACY<rfactorlegacy>`
   (by default, the TensErLEED one; see there also for the restrictions
   of the python-native one). Both implementations support all four
   |R-factor| types. The TensErLEED |R factor| supports |RS| only with
   TensErLEED v2.1.0 or later (see :ref:`TL_VERSION<tl_version>`): with
   older versions, a warning is issued and |RP| is calculated instead.
-  The TensErLEED :ref:`structure search<sec_search>` always uses the
   TensErLEED |R factor|, whatever the value of R_FACTOR_LEGACY. The
   search supports |RP|, :math:`R_2`, and |RS| (the latter only with
   TensErLEED v2.1.0 or later). The search stops with an error for the
   Zanazzi-Jona |R factor|, and for |RS| with older TensErLEED versions.

For how the implementations differ in :math:`R_2`, see
:ref:`R_FACTOR_LEGACY<rfactorlegacy>`.

Note that var(R) in :ref:`error calculations<error_calculation>` is only
supported for the Pendry and smooth |R factors|.

.. versionchanged:: 0.15.0
   R_FACTOR_TYPE accepts names and synonyms of the |R-factor| types (see
   above). Earlier versions only accept the integer numbers 1 (Pendry) and
   2 (:math:`R_2`). The integer codes 3 (Zanazzi-Jona) and 4 (smooth) and
   the smooth and Zanazzi-Jona |R-factor| types are new in version 0.15.0.
