.. include:: /substitutions.rst

.. _rfactorsmooth:

R_FACTOR_SMOOTHING
==================

R_FACTOR_SMOOTHING defines the number of times that experimental |IV| curves
will be smoothed by 3-point smoothing during |R-factor| calculations, both
after the reference calculation and during the search.

.. note::
    It is recommended to smooth the experimental curves beforehand, such
    that smoothing during TensErLEED operation is not necessary. See the
    |R-factor| :ref:`page <r-factor_calculation>` for details.
    
    If you do use smoothing through the R_FACTOR_SMOOTHING parameter, check the
    experimental curves in the :ref:`Rfactor_analysis.pdf<rfactoranalysis>`
    file. This plots the curves after smoothing, as well as the Y-function used
    for calculating the |R factor|.


.. versionchanged:: 0.15.0
   Parameter renamed from ``R_FACTOR_SMOOTH`` to ``R_FACTOR_SMOOTHING``.
   The old name is still accepted as an alias for backwards compatibility.

**Default:** 0

**Allowed values:** 0 or positive integer

**Syntax:**

::

   R_FACTOR_SMOOTHING = 4

.. warning::
   This parameter is only applied in the old TensErLEED-based fortran codes,
   and not in the newer |R-factor| implementations. If you use the newer,
   python-based |R-factor| implementations, for example by setting 
   :ref:`R_FACTOR_LEGACY<rfactorlegacy>` to ``False``, this parameter will have
   no effect.