.. _vibdelta:

================================
Vibrational-amplitude variations
================================

Vibration amplitude variations to be used during the
search must be specified in a block starting with the

::

   = VIB_DELTA

header flag, followed by a list of displacements for each of the atoms
for which the vibrational amplitude should be varied during the search.

Example
=======

..  code-block:: none

   = VIB_DELTA
   O 1 = -0.05 0.05 0.02         ! Vibrational amplitude of oxygen atom 1 (and symmetry-equivalent atoms) will be varied over the range [-0.05, 0.05] with step 0.02
   Ir_top = -0.05 0.05 0.01      ! Vibrational amplitude of Iridium atoms in Ir_top sites will be varied over the range [-0.05, 0.05] with step 0.01

How atoms are addressed on the left is described on the main
:ref:`DISPLACEMENTS<DISPLACEMENTS>`  page. The values on the
right will be interpreted relative to the atom's vibrational
amplitude as defined in the :ref:`VIBROCC file<vibrocc>`.

For some applications, it can be useful to apply a static displacement,
without re-doing the reference calculation. For this purpose, the
VIB_DELTA block also accepts single-value input on the right:

..  code-block:: none

   = VIB_DELTA
   O 1 = 0.02         ! Vibrational amplitude of oxygen atom 1 (and symmetry-equivalent atoms) will be offset from the value in VIBROCC by 0.02

When multiple searches are executed consecutively or in a loop, the
displacement ranges are per default centered around the optimized
vibrational amplitude from previous searches per default. However,
if you give a single-valued (static) displacement for an atom, the
optimized vibrational amplitude from previous searches is *discarded*
instead, and the static displacement is applied to its original position
(to avoid "displacement creep" when the search is repeated). If you want
to center the displacement range around the original vibrational amplitude
of the atom, you can also clear the offset manually:

..  code-block:: none

   O 1 offset = 0         ! center range around original vibrational amplitude, instead of the optimized vibrational amplitude resulting from previous searches
   O 1 offset = clear     ! equivalent
   O 1 offset = original  ! equivalent


.. note::
   -  The range of vibrational amplitudes defined with the command will be
      applied to all symmetry-equivalent atoms, unless you turn off symmetry
      constraints via :ref:`SYMMETRY_FIX<ISYM>`  or :ref:`SYM_DELTA<SYMDELTA>`.
   -  Vibrational amplitudes must be *positive* values, so the program will
      throw an error if the combination of the given range with the center
      point defined in the :ref:`VIBROCC file<vibrocc>`  results in negative
      values.
   -  If your combination of ``start``, ``end`` and ``step`` values does not
      yield an odd number of steps, the interval will be expanded symmetrically
      around the midpoint.
