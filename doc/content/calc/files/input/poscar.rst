.. include:: /substitutions.rst

.. |POSCAR_bulk|     replace:: :file:`POSCAR_bulk`
.. |POSCAR_mincell|  replace:: :file:`POSCAR_mincell`
.. |POSCAR_ori|      replace:: :file:`POSCAR_ori`
.. |POSCAR_oricell|  replace:: :file:`POSCAR_oricell`
.. |POSCAR_user|     replace:: :file:`POSCAR_user`
.. |POSCAR_vacuum|   replace:: :file:`POSCAR_vacuum_corrected`

.. _poscar:

POSCAR
======

A |POSCAR| file describes the structure: The unit cell, the atom types (and the
number of atoms for each type), and their coordinates. Note that the atom names
in a |POSCAR| file ("POSCAR elements") need not be actual chemical elements (as
given in the periodic table).
See the page about :ref:`element name collision<ElementNameCollision>` for
the distinction between POSCAR elements and chemical elements, and the
:ref:`ELEMENT_MIX` and :ref:`ELEMENT_RENAME` parameters for mapping between
POSCAR-element names and chemical elements.

The POSCAR file format is the same as used in :program:`VASP`, and can be
exported directly from :program:`VESTA`
:cite:p:`mommaVESTAThreedimensionalVisualization2011`
(\ :guilabel:`File -> Export Data...`, pick :guilabel:`VASP (POSCAR;*.vasp)`
in the :guilabel:`Save as type` dropdown menu, then select
:guilabel:`Fractional coordinates`).
The atoms in the |POSCAR| file are ordered by element, but otherwise, their
order does not matter. In the data lines containing coordinates, only the
first three columns are read by ViPErLEED, so anything to the right will
be ignored. |POSCAR| files with ``Cartesian`` coordinates will be accepted,
but ViPErLEED uses fractional (``Direct``) coordinates in all output files,
so it is recommended to follow this convention in the input as well.

.. seealso::
    `The VASP wiki <https://www.vasp.at/wiki/index.php/POSCAR>`__ describes
    in detail the POSCAR format.

ViPErLEED has some requirements concerning the orientation of the structure
in the |POSCAR| file. The first two unit-cell vectors |a| and |b| must lie in
the surface plane. This means that their |z| component (perpendicular to the
surface) must be zero. The third unit-cell vector |c| must have a non-zero
component in the |z| direction, but does not necessarily need to be
perpendicular to the surface. Slabs must be asymmetric, with +\ |z| pointing
away from the surface into vacuum, and the lowest-lying layers (i.e., smallest
|z| coordinates) bulk-like.

After the |calc| :ref:`initialization<initialization>` is run for the first
time, the original |POSCAR| file will be copied to |POSCAR_user|.
The following changes are then made to |POSCAR|:

-  The unit cell will be simplified to a higher-symmetry / shortest-perimeter
   form, if possible. If this happens, a warning will appear in the log
   together with the transformation matrix.
-  Rhombic, oblique, and hexagonal cells will be transformed from acute to
   obtuse.
-  The origin will be shifted to a high-symmetry point, that is, the
   highest-order rotation axis, or — if no rotation axis is found — a
   mirror or glide plane. See the "conventional origin" mark in
   :ref:`here<planegroups>` for which position is used for each
   plane group.
-  For atoms that were recognized as symmetry-equivalent within
   :ref:`sym_eps`, the atomic positions will be averaged to fully
   correspond to the system's symmetry (using either an automatically
   determined plane group, or the one defined in :ref:`ISYM`). This
   behavior can be altered with the :ref:`SYMMETRIZE_INPUT` parameter.
-  Atoms that lie within :ref:`sym_eps`  of a rotation
   axis or mirror plane will be moved onto that axis or plane to fully
   correspond to the system's symmetry. This can also be prevented via
   :ref:`SYMMETRIZE_INPUT`.
-  If a symmetry reduction that requires rotation of the unit cell has
   been set in the :ref:`ISYM` parameter, the unit cell will be rotated
   in the |POSCAR| file.
-  Comments will be added in the |POSCAR| file, which predict the behavior
   of the system during the subsequent calculations (see
   :ref:`below <poscar_comments>`).

This modified |POSCAR| file can also be found in the |OUT| folder after
:ref:`initialization`.

At the end of each |calc| execution, the |POSCAR| file given as input for that
run is renamed to |POSCAR_ori| (this is identical to |POSCAR_user| after the
very first invocation), while the edited file is copied to the root directory
(from |OUT|) as a new |POSCAR| file. See also :ref:`poscar_out`. This ensures
that further invocations of |calc| will automatically use the output of
previous executions as an input. You can manually call the |bookkeeper|
utility after a specific |calc| run if this behavior is not desirable.
See the :ref:`bookkeeper` page for more details.

.. note::
    A non-halted execution (i.e., one where :ref:`halting` was set to a value
    larger than the default) that includes a structure optimization will
    overwrite the :ref:`poscar_out` file created during :ref:`initialization`
    with the one found by the (last) optimization step. This is also the
    case for the |POSCAR| file found in the root directory a the end of
    such a run.


.. _poscar_comments:

Comment lines
-------------

The |POSCAR| file contains the following **comment lines** after
initialization:

Plane group
    Planar symmetry group of the slab. If :ref:`ISYM` is used to select a
    certain group, this will be indicated, while the full symmetry of the
    slab is mentioned in brackets. In the |POSCAR_oricell| file, this is
    prepended by stars, since the original cell might not display the
    correct symmetry. See also the :ref:`list of plane groups <planegroups>`.

The atoms are then listed one per line, grouped by element.
For each atom the following information is given:

N
    Consecutive numbering of the atoms. Same as atom number in :program:`VESTA`
    :cite:p:`mommaVESTAThreedimensionalVisualization2011`. Atom numbering is
    conserved from the original |POSCAR|. This numbering convention is applied
    everywhere in |calc|.

SiteLabel
    ``element_sitetype``, as determined from :ref:`SITEDEF`.
Layer
    The layer that the atom is in, as determined from :ref:`LAYER_CUTS`.

Linking
    Progressive label that indicates which atoms are related to one another by
    the **Plane group**. Atoms belonging to the same equivalence group share
    the same number. When one of the atoms from an equivalence group is moved
    via the :ref:`DISPLACEMENTS`, its equivalent ones will be also moved such
    that the symmetry is conserved (see the :ref:`DISPLACEMENTS` file for
    further details).

FreeDir
    Allowed in-plane movement direction for the atom during LEED optimization.
    Will be ``locked`` if the atom is on a rotation axis, and ``[i j]`` if the
    atom is on a mirror plane, where the allowed direction is
    ``i``\ |a| + ``j``\ |b|. This column is not displayed in the
    |POSCAR_oricell| file, since the cell (and therefore the unit
    vectors) might be different. Bulk atoms will be labelled as
    ``bulk`` in this column, since they cannot be moved during
    optimization.


.. _poscar_out:

OUT/POSCAR
----------

After running the :ref:`initialization`, a |POSCAR| file can be found in
the |OUT| folder. This is an edited version of the user-given |POSCAR|,
as described :ref:`above<poscar>`.

After executing a structure optimization (i.e., :ref:`search<search>` or
:ref:`fdoptimization`), the |POSCAR| file in |OUT| corresponds to the one
that realizes the best (i.e., smallest) |R factor|. It has the same format
as the one after initialization.

At the end of each |calc| run, the :file:`OUT/POSCAR` file is copied to
the root folder as a new |POSCAR| file. The |POSCAR| file given as input
for that run is renamed to |POSCAR_ori|. This is such that the next |calc|
run will use as input the output of the previous one. You can call the
|bookkeeper| utility after a specific |calc| run if this behavior is
not desirable. See the :ref:`bookkeeper` page for more details.

.. versionchanged:: 0.13.0
    In earlier versions of |calc|, the automatically edited |POSCAR| file would
    only appear in the root directory after :ref:`initialization`, and only the
    |POSCAR| file resulting from a structural optimization would be stored in
    |OUT|. This file used to be named :file:`POSCAR_OUT`.


.. _poscar_oricell:

POSCAR_oricell
--------------

A separate |POSCAR_oricell| file is created (see |SUPP| folder), which contains
comments and corrections of atomic positions, but with the same orientation and
position of the unit cell as in the original |POSCAR|. This can be used for
direct comparison (e.g., in :program:`VESTA`
:cite:p:`mommaVESTAThreedimensionalVisualization2011`) with the original
file, and can be useful to judge whether the :ref:`sym_eps` value chosen
is appropriate.


.. _poscar_bulk:

POSCAR_bulk
-----------

In addition, a |POSCAR_bulk| file is created (see |SUPP| folder) based
on the :ref:`LAYER_CUTS`, :ref:`N_BULK_LAYERS`, :ref:`BULK_REPEAT`
and :ref:`SUPERLATTICE` parameters. The structure in this file represents
the bulk as it will be used during the TensErLEED calculations. You can
check whether the bulk unit cell was recognized correctly by opening
|POSCAR_bulk| in :program:`VESTA` and editing the boundary such that
multiple cells are shown in all three directions. For the same purpose,
the :file:`POSCAR_bulk_appended` is the original |POSCAR| file with 1–3
bulk units (depending on the bulk thickness) appended at the bottom. It
is meant to help checking whether the bulk cell is aligned correctly with
the slab.


.. _poscar_mincell:

POSCAR_mincell
--------------

If the :ref:`SYMMETRY_CELL_TRANSFORM` parameter is set, or if a
smaller-area unit cell is found during the symmetry search, an
additional |POSCAR_mincell| file will be written, containing the
atoms in the reduced unit cell.


.. _poscar_vacuum_corrected:

POSCAR_vacuum_corrected
-----------------------

A |POSCAR_vacuum| file is provided (in folder |SUPP|) if the original |POSCAR|
file does not have a suitable vacuum gap (> 5 Å) between its topmost and (a
|c|-periodic replica of) its bottommost atom. The following |POSCAR| input
files will be considered unsuitable for ViPErLEED:

- The vacuum gap is somewhere in the middle. This means that there are multiple
  atoms above a large (> 5 Å) vacuum gap. A typical example is a 'symmetric'
  slab centered at :math:`c=0`;
- There are atoms very close (:math:`< 1\times10^{-4}` in fractional
  coordinates) to both :math:`c=0` and :math:`c=1`.

In these cases, the |POSCAR_vacuum| file may be used as a starting point to
produce an acceptable input |POSCAR| for a subsequent run.

.. note::
    When preparing a new set of input files from |POSCAR_vacuum|, be careful to
    adapt any :ref:`PARAMETERS` that are defined as fractions of the unit-cell
    |c| vector (e.g., :ref:`LAYER_CUTS`, :ref:`BULK_LIKE_BELOW`,
    :ref:`BULK_REPEAT`).

A |POSCAR| file with a gap smaller than 5 Å will not cause ViPErLEED to stop,
but a |POSCAR_vacuum| file is nonetheless generated. It can be used, e.g., to
judge the appropriate identification of which atoms are at the top and which
ones belong to the bulk-like portion of the input |POSCAR|.
