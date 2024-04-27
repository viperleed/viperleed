.. _symmetry_nomove:

SYMMETRIZE_INPUT
================

SYMMETRIZE_INPUT can be used to suppress recalculation of atom coordinates to perfectly match a given symmetry. How much atoms might otherwise move is defined by the tolerance :ref:`SYMMETRY_EPS<sym_eps>`.

**Default**: SYMMETRIZE_INPUT = True (atom positions might be changed during symmetry search)

**Syntax**:

::

   SYMMETRIZE_INPUT = False

**Acceptable values**: True, False, true, false, T, F, t, f (not case sensitive)

SYMMETRIZE_INPUT might be used when a system does not entirely fulfill a given symmetry (and should not be "corrected" to that symmetry), but the user still wants to use symmetry detection with a large :ref:`SYMMETRY_EPS<sym_eps>`  in order to automatically move almost-symmetry-equivalent atoms in concert. Similarly, SYMMETRIZE_INPUT could be used to avoid atoms snapping to a nearby rotation axis or mirror plane. **Use of SYMMETRIZE_INPUT is discouraged, unless you know what you are doing!**

Note that SYMMETRIZE_INPUT does **not** prevent atoms from being symmetry-linked during the search. If that is your goal, the parameter you want is :ref:`SYMMETRY_FIX<ISYM>`. If the symmetry is reduced or deactivated by :ref:`SYMMETRY_FIX<ISYM>`, symmetrization of atom positions will only be performed for the reduced symmetry even if SYMMETRIZE_INPUT is True.
