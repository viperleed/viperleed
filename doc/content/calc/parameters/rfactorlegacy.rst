.. include:: /substitutions.rst

.. _rfactorlegacy:

R_FACTOR_LEGACY
===============

R_FACTOR_LEGACY toggles between using the old (run-time compiled) TensErLEED
|R factor| and a new python-native |R-factor| calculation implemented in
ViPErLEED. The new implementation uses spline interpolation, which yields
more stable (second) derivatives. However, this may result in slightly different
|R factor| values even when using the same :ref:`rfactortype`.

By default, the TensErLEED |R factor| is used. The default will change
to the python-native |R factor| once ViPErLEED supports structure
searches that do not use TensErLEED. R_FACTOR_LEGACY only selects the
|R factor| reported after a reference calculation or a superposition
calculation (and for :ref:`error curves<error_calculation>`). The
TensErLEED structure search always uses the TensErLEED |R factor|,
whatever the value of R_FACTOR_LEGACY.

The new implementation supports only cubic interpolation. Combining
R_FACTOR_LEGACY = ``False`` with an INTPOL_DEG other than 3 is rejected
when the parameters are read. Set R_FACTOR_LEGACY to ``True`` to use other
interpolation degrees.

.. versionchanged:: 0.15.0
   INTPOL_DEG other than 3 is rejected if R_FACTOR_LEGACY is ``False``.

**Default:** True

**Allowed values:** True (old TensErLEED), False (new)

**Syntax:**

::

   R_FACTOR_LEGACY = True
