.. _symmetry_cell_transform:

SYMMETRY_CELL_TRANSFORM
=======================

SYMMETRY_CELL_TRANSFORM defines the base unit cell on which symmetry
operations should be found and used for linking atom movements. All
other atoms are assumed to be translationally symmetric, and will be
linked to the base cell accordingly. The relationship between the
symmetry base cell and the :ref:`POSCAR<POSCAR>`  unit cell is
defined by SYMMETRY_CELL_TRANSFORM in the same way as the relationship
between the bulk unit cell and the :ref:`POSCAR<POSCAR>`  unit cell
is defined by :ref:`SUPERLATTICE<SUPERLATTICE>`.

.. note::
    SYMMETRY_CELL_TRANSFORM may be updated automatically during initialization
    if a reducible supercell is detected.

**Default:** (1x1)

**Syntax examples:** see :ref:`SUPERLATTICE page<SUPERLATTICE>`

SYMMETRY_CELL_TRANSFORM is relevant in the context of
:ref:`domain calculations<domain_calculation>`, where the calculations for the
structures in the different domains must all be run with the same supercell,
even when their own unit cell is smaller. For example, when fitting a mixture
of a (1x1) and a (2x1) structure, the (1x1) structure has to be re-calculated
on a (2x1) unit cell. However, during this calculation, it makes sense to still
base the symmetry search for the (1x1) on the base (1x1) unit cell instead of
the (2x1), which may lower symmetry. In addition, atoms that are completely
equivalent on the (1x1) should be linked by translational symmetry on the (2x1)
supercell. SYMMETRY_CELL_TRANSFORM performs this function, essentially defining
a 'real' surface unit cell which should be used as a basis for atom linking
(which may differ from the bulk *and* from the POSCAR unit cell). Note that
this effectively makes the supercell calculation identical to a calculation
on the smaller cell, i.e. the resulting supercell structure is reducible to
the base cell defined by SYMMETRY_CELL_TRANSFORM without any conflicts or
loss of information. Therefore, outside of the context of domain calculations,
it is usually better to use a smaller POSCAR unit cell rather than the
SYMMETRY_CELL_TRANSFORM parameter, because the results will be the same,
but calculating on a larger cell is much less computationally efficient.

Note that the base unit cell defined by SYMMETRY_CELL_TRANSFORM may still
be automatically reshaped to a higher-symmetry / lower-circumference form,
in the same way as during :ref:`optimization of the POSCAR unit cell<POSCAR>`.

**TODO Florian (@Michele, Alex: Florian said this is on his TODO list)**: add
figures to help
