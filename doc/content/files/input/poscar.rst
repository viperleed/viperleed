.. include:: /substitutions.rst

.. _poscar:

POSCAR
======

A POSCAR file describes the structure: The unit cell, the atom types (and the
number of atoms for each type), and their coordinates. Note that the atom names
in a POSCAR file ("POSCAR elements") need not be actual chemical elements (as
given in the periodic table).
See :ref:`element name collision<ElementNameCollision>`  for the distinction
between POSCAR elements and chemical elements, and the
:ref:`ELEMENT_MIX<ELSPLIT>`  and :ref:`ELEMENT_RENAME<ELDEF>`
parameters for mapping between POSCAR element names and chemical elements.

The POSCAR file format is the same as used in VASP input, and can be exported
directly from VESTA :cite:p:`mommaVESTAThreedimensionalVisualization2011`
(File -> Export Data -> file type 'VASP (POSCAR;*.vasp)' -> Fractional
coordinates). The atoms in the POSCAR file are ordered by element, but
otherwise, their order does not matter. In the data lines containing
coordinates, only the first three columns are read by ViPErLEED, so
anything to the right will be ignored. POSCAR files with "Cartesian"
coordinates will be accepted, but ViPErLEED uses fractional ("Direct")
coordinates in all output files, so it is recommended to follow this
convention in the input as well.

**See also:** `POSCAR in the VASP wiki <https://www.vasp.at/wiki/index.php/POSCAR>`__

ViPErLEED has some requirements concerning the orientation of the structure
in the POSCAR file. The first two unit cell vectors a and b must be parallel
to the surface plane, i.e. their z component (perpendicular to the surface)
must be zero. The third unit cell vector c must have a non-zero component
in the z direction, but does not necessarily have to be perpendicular to the
surface. Slabs must be asymmetric, with +z pointing away from the surface and
the lowest-lying layers (i.e. smallest z coordinates) bulk-like.

After the |calc| :ref:`initialization<initialization>` is run for the first
time, the original POSCAR file will be copied to POSCAR_user. The following
changes are then made to POSCAR:

-  The unit cell will be simplified to a higher-symmetry / lower-circumference
   form, if possible. If this happens, a warning will appear in the log
   together with the transformation matrix.
-  Rhombic, oblique and hexagonal cells will be transformed from acute to
   obtuse
-  The origin will be shifted to a high-symmetry point, that is, the
   highest-order rotation axis, or if no rotation axis is found, a
   mirror or glide plane.
-  For atoms that were recognized as symmetry-equivalent within
   :ref:`SYMMETRY_EPS<sym_eps>`, the atomic positions will be
   averaged to fully correspond to the system's symmetry (using
   either an automatically determined plane group, or the one
   defined in :ref:`SYMMETRY_FIX<ISYM>`). This behavior can be
   altered with the :ref:`SYMMETRIZE_INPUT<SYMMETRY_NOMOVE>`
   parameter.
-  Atoms that lie within :ref:`SYMMETRY_EPS<sym_eps>`  of a rotation
   axis or mirror plane will be moved onto that axis or plane to fully
   correspond to the system's symmetry. This can also be prevented via
   :ref:`SYMMETRIZE_INPUT<SYMMETRY_NOMOVE>`.
-  If a symmetry reduction that requires rotation of the unit cell has
   been set in the :ref:`SYMMETRY_FIX<ISYM>`  parameter, the unit cell
   will be rotated in the POSCAR.
-  Comments will be added in the POSCAR file, which predict the behaviour
   of the system during the subsequent calculations (see below).

Comment lines
-------------

The POSCAR file contains the following **comment lines** after initialization:

-  **Plane group**: planar symmetry group of the slab.
   If :ref:`SYMMETRY_FIX<ISYM>` is used to select a certain group, this will be
   indicated, while the full symmetry of the slab is mentioned in brackets. In
   the POSCAR_oricell file, this is prepended by stars, since the original cell
   might not display the correct symmetry. See also the
   :ref:`list of planegroups <planegroups>`.

The atoms are then listed one per line, grouped by element.
For each atom the following information is given:

-  **N**: Consecutive numbering of the atoms. Same as atom number in VESTA
   :cite:p:`mommaVESTAThreedimensionalVisualization2011`. Atom numbering is
   conserved from the original POSCAR. This numbering convention is applied
   everywhere in |calc|.
-  **SiteLabel**: ``element_sitetype``, as determined from
   :ref:`SITE_DEF<SITEDEF>`.
-  **Layer**: The layer that the atom is in, as determined from
   :ref:`LAYER_CUTS<layer_cuts>`.
-  **Linking**: Progressive label that indicates which atoms are related to one
   another by the symmetry **Group**. When one of the atoms from an equivalence
   group is moved via the :ref:`DISPLACEMENTS<DISPLACEMENTS>`, its equivalent
   ones will be also moved such that the symmetry is conserved (see the
   :ref:`DISPLACEMENTS<DISPLACEMENTS>`  file for further details).
-  **FreeDir**: Allowed in-plane movement direction for the atom during LEED
   optimization. Will be ``locked`` if the atom is on a rotation axis, and
   ``[i j]`` if the atom is on a mirror plane, where the allowed direction is
   ``ia + jb``. This column is not displayed in the POSCAR_oricell file, since
   the cell (and therefore the unit vectors) might be different. Bulk atoms
   will be labelled ``bulk`` in this column, since they cannot be moved during
   optimization.

.. _poscar_oricell:

POSCAR_oricell
--------------

A separate **POSCAR_oricell** file is created (see SUPP folder), which contains
comments and corrections of atomic positions, but with the same orientation and
position of the unit cell as in the original POSCAR.
This can be used for direct comparison (e.g., in VESTA
:cite:p:`mommaVESTAThreedimensionalVisualization2011`) with the original file,
and can be useful to judge whether the :ref:`SYMMETRY_EPS<sym_eps>`  value
chosen is appropriate.

.. _poscar_bulk:

POSCAR_bulk
-----------

In addition, a **POSCAR_bulk** file is created (see SUPP folder) based on
the :ref:`LAYER_CUTS<layer_cuts>`, :ref:`N_BULK_LAYERS<n_bulk_layers>`,
:ref:`BULK_REPEAT<BULK_REPEAT>` and :ref:`SUPERLATTICE<SUPERLATTICE>`
parameters. The structure in this file represents the bulk as it will be
used during the TensErLEED calculations. You can check whether the bulk
unit cell was recognized correctly by opening POSCAR_bulk in VESTA and
editing the boundary such that multiple cells are shown in all three
directions. For the same purpose, the **POSCAR_bulk_appended** is the
original POSCAR file with 1–3 bulk units (depending on the bulk thickness)
appended at the bottom, meant to check whether the bulk cell is aligned
correctly with the slab.

.. _poscar_mincell:

POSCAR_mincell
--------------

If the :ref:`SYMMETRY_CELL_TRANSFORM<SYMMETRY_CELL_TRANSFORM>`  parameter
is set, or if a smaller-area unit cell is found during the symmetry search,
an additional **POSCAR_mincell** file will be written, containing the atoms
in the reduced unit cell.


.. _poscar_out:

POSCAR_OUT
----------

After executing a search, a POSCAR_OUT file will be produced in the OUT folder.
This takes the same format as the POSCAR file after initialization, and the new
positions are those of the best-fit structure found during the search (i.e.,
corresponding to the lowest |R factor|).


.. _poscar_vacuum_corrected:

POSCAR_vacuum_corrected
-----------------------

A **POSCAR_vacuum_corrected** file is provided (in folder SUPP) if the original
POSCAR file does not have a suitable vacuum gap (> 5 Å) between its topmost and
(a **c**-periodic replica of its) bottommost atom. The following POSCAR input
files will be considered unsuitable for ViPErLEED:

* The vacuum gap is somewhere in the middle. This means that there are multiple
  atoms above a large (> 5 Å) vacuum gap. A typical example is a 'symmetric'
  slab centred at :math:`c=0`;
* There are atoms very close (:math:`< 1\times10^{-4}` in fractional
  coordinates) to both :math:`c=0` and :math:`c=1`.

In these cases, the POSCAR_vacuum_corrected file may be used as a starting
point to produce an acceptable input POSCAR for a subsequent run.

.. note::
    When preparing a new set of input files from POSCAR_vacuum_corrected, be
    careful to adapt any :ref:`PARAMETERS<parameters>` that are defined as
    fractions of the unit-cell c vector (e.g., :ref:`LAYER_CUTS<layer_cuts>`,
    :ref:`BULK_LIKE_BELOW<BULK_LIKE_BELOW>`, :ref:`BULK_REPEAT<BULK_REPEAT>`).

A POSCAR file with a gap smaller than 5 Å will not cause ViPErLEED to stop, but
a POSCAR_vacuum_corrected file is nonetheless generated. It can be used, e.g.,
to judge the appropriate identification of which atoms are at the top and which
ones belong to the bulk-like portion of the input POSCAR.
