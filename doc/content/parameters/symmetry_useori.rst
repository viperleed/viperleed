.. _symmetry_useori:

SYMMETRY_FIND_ORI
=================

SYMMETRY_FIND_ORI = False is used to declare that the origin of the input :ref:`POSCAR<POSCAR>`  file is already one of the highest-symmetry points in the cell. The symmetry search is then restrained to the origin. Only the origin will be tested as a rotation axis, and only planes passing through the origin will be tested as the main mirror or glide planes.

**Default**: SYMMETRY_FIND_ORI = True if no symmetry group is defined in the :ref:`POSCAR<POSCAR>`  header, False otherwise

**Syntax**:

::

   SYMMETRY_FIND_ORI = False

**Acceptable values**: True, False, true, false, T, F, t, f

SYMMETRY_FIND_ORI = False is intended for user input files only in very large cells that do **not** have at least one low-occupancy sublayer. In such cells, the symmetry search can sometimes time out due to the very high number of candidate symmetries that need to be tested. If the user is aware of the highest-symmetry point in the cell, moving the origin there beforehand and setting SYMMETRY_FIND_ORI = False is recommended.

For POSCAR files optimized by ViPErLEED, the origin has been moved to a point of highest symmetry. These files are recognized by the POSCAR header, which includes the symmetry group.
