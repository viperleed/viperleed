.. _symmetrybulk:

SYMMETRY_BULK
=============

SYMMETRY_BULK allows you to manually set the symmetry to be used for the bulk,
ignoring the automatic detection during initialization. This symmetry is *only*
used to determine the appropriate beam averaging, in case of terraces with
different orientations (these are assumed to cover the same area fraction).

**Default**: Unused; bulk symmetry will be detected automatically from the
input structure.

**Syntax:**

::

   SYMMETRY_BULK = p6          ! use plane group p6 for the bulk. Only allowed for hexagonal bulk cells.
   SYMMETRY_BULK = cm[1 1]     ! use plane group cm, with mirrors along the short diagonal. Only allowed for rhombic, square or hexagonal bulk cells.
   SYMMETRY_BULK = cm[1 1] r4  ! same as above, with an additional 4-fold rotation symmetry. This combination is only allowed for square cells.
   SYMMETRY_BULK = p2 m[0 1]   ! use plane group p2, with an additional mirror along b. Only allowed for rectangular or square cells.

**Accepted values**: At least one entry must be a valid plane group (see
:ref:`SYMMETRY_FIX<ISYM>` and :ref:`plane symmetry groups<planegroups>`).
Additional rotational symmetry can be specified by ``ri``, with ``i`` the
order for rotation. Additional mirror planes can be specified using syntax
``m[i1 i2]``, where ``i1`` and ``i2`` specify a direction for the mirror
as ``i1``\ ×\ **a**\ :sub:`bulk` + ``i2``\ ×\ **b**\ :sub:`bulk` (with
the bulk unit cell vectors **a**\ :sub:`bulk`, **b**\ :sub:`bulk`).

The symmetry group must be a valid group for the bulk unit cell type
(e.g., a rectangular cell cannot be p3). Additional mirrors and
rotations must also be valid operations for the bulk unit cell type.

Since the bulk cannot be modified during the search, there is no symmetry
linking, and the bulk symmetry is only ever used for determining which beams
are equivalent. There are two major use cases for SYMMETRY_BULK:

-  To reduce symmetry: When two growth directions of an overlayer are possible
   on the given substrate, but only one is actually present in experiment, bulk
   rotation should be switched off to avoid automatic averaging.
-  To increase symmetry: If a thick film is used for LEED experiments, the
   substrate may not be present in the POSCAR (since it does not contribute
   to scattering). However, a higher-symmetry substrate would still results
   in multiple orientations of the film. These higher symmetries due to a
   substrate not present in the POSCAR can be introduced by SYMMETRY_BULK.
   Note however that this currently does not work for arbitrary substrate/film
   relationships. Ideally, the unit cell type of the substrate should be the
   same as that of the film bulk.
