.. include:: /substitutions.rst

.. _ivshiftrange:

IV_SHIFT_RANGE
==============

IV_SHIFT_RANGE defines how much calculated and experimental |IV| curves are
allowed to shift (in electronvolts) during calculation of the |R factor|, and
(optionally) the step width.

**Default**: from -3 to +3 eV, with step from experimental or calculated curves
(whichever is smaller)

**Syntax:**

::

   IV_SHIFT_RANGE = -5 5         ! range -5 to +5 eV, with step determined from data
   IV_SHIFT_RANGE = 0 0          ! suppress shifting
   IV_SHIFT_RANGE = -4 4 0.1     ! range -4 to +4 eV, with 0.1 eV step
   IV_SHIFT_RANGE = _ _ 0.2      ! default range, but use 0.2 eV step

**Acceptable values:** Either 2 or 3 float values. Step has to be positive.
Underscore is interpreted as default value.

**See also**: :ref:`FILAMENT_WF<FILWF>`  parameter (filament work function,
shifts the electron energy with respect to the Fermi level).
