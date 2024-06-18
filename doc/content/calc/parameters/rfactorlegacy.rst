.. include:: /substitutions.rst

.. _rfactorlegacy:

R_FACTOR_LEGACY
===============

R_FACTOR_LEGACY toggles between using the old (run-time compiled) TensErLEED
|R factor| and a new experimental ViPErLEED |R-factor| calculation. By default,
the old TensErLEED |R factor| is used while the new version is undergoing
tests. It is recommended to keep the default value until further notice.

**Default:** True

**Allowed values:** True (old TensErLEED), False (new)

**Syntax:**

::

   R_FACTOR_LEGACY = True
