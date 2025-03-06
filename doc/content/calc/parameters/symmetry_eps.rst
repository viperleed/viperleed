.. _sym_eps:

SYMMETRY_EPS
============

SYMMETRY_EPS defines the minimum distance (in angstrom) that is used to
determine whether two atoms occupy symmetry-equivalent positions. If only
one value is given, this value is used both for in-plane comparison and
for determining which atoms are in the same plane. If two values are given,
the first is used for in-plane comparison, while the second is used for Z
comparison.

**Default**: SYMMETRY_EPS = 0.1 (i.e., 0.1 angstrom)

**Syntax**:

::

   SYMMETRY_EPS = 0.2
   SYMMETRY_EPS = 0.1 0.05

**Acceptable values**: One or two floating point values greater than 0.
A warning will be displayed for values > 1.0

For atoms that are recognized as symmetry-equivalent within SYMMETRY_EPS,
atomic positions will be averaged during initialization to fully reflect
that symmetry. The choice of SYMMETRY_EPS, in combination with :ref:`ISYM`,
will therefore determine how strongly the atom positions in the :ref:`POSCAR`
file are modified during initialization. You can use :ref:`SYMMETRIZE_INPUT`
to suppress symmetrization.
