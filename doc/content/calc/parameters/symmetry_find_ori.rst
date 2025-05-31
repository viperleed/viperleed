.. include:: /substitutions.rst

.. _symmetry_find_ori:

SYMMETRY_FIND_ORI
=================

SYMMETRY_FIND_ORI = False is used to declare that the origin of the input
:ref:`POSCAR` file is already one of the highest-symmetry points in the cell.
The symmetry search is then restrained to the origin. Only the origin will be
tested as a rotation axis, and only planes passing through the origin will be
tested as the main mirror or glide planes.

**Default**: SYMMETRY_FIND_ORI = True

**Syntax**:

::

   SYMMETRY_FIND_ORI = False

**Acceptable values**: True, False, true, false, T, F, t, f


.. versionchanged:: 0.8.0
   Up to and including |calc| version 0.7.2 the default was only True if no
   symmetry group was defined in the :ref:`POSCAR` header and False otherwise.
   This was changed in newer versions as it could introduce errors if the
   POSCAR was edited, but the symmetry group comment was left unchanged.


``SYMMETRY_FIND_ORI = False`` is intended for user input files only in very
large cells that do **not** have at least one low-occupancy sublayer. In such
cells, the symmetry search can sometimes time out due to the very high number
of candidate symmetries that need to be tested. If the user is aware of the
highest-symmetry point in the cell, moving the origin there beforehand and
setting ``SYMMETRY_FIND_ORI = False`` is recommended.

For POSCAR files optimized by ViPErLEED, the origin has been moved to a point
of highest symmetry. These files can be recognized by the POSCAR header, which
includes the symmetry group.