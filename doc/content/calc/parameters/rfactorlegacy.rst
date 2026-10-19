.. include:: /substitutions.rst

.. _rfactorlegacy:

R_FACTOR_LEGACY
===============

.. versionchanged:: 0.15.0
   Since version 0.15.0, the default is to use the new python-native |R factor|
   calculation.


R_FACTOR_LEGACY toggles between using the old (run-time compiled) TensErLEED
|R factor| and a new python-native |R-factor| calculation implemented in
ViPErLEED. The new implementation uses spline interpolation, which yields
more stable (second) derivatives. However, this may result in slightly different
|R factor| values even when using the same :ref:`rfactortype`.

**Default:** False

**Allowed values:** True (old TensErLEED), False (new)

**Syntax:**

::

   R_FACTOR_LEGACY = True
