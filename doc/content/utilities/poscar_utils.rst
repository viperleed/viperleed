.. _poscar_utils:

POSCAR utilities
################

All POSCAR utilities support the following options:

- ``-h``, ``--help``: print help message with the available options
- ``-v``, ``--verbose``: increase output verbosity

.. _poscar_utils_attach_bulk:

attach_bulk
===========

.. _poscar_utils_delete_above:

delete_above
============

.. _poscar_utils_delete_below:

delete_below
============

.. _poscar_utils_delete_between:

delete_between
==============

.. _poscar_utils_enforce_symmetry:

enforce_symmetry
================

.. _poscar_utils_find_symmetry:

find_symmetry
=============

.. _poscar_utils_get_bulk_repeat:

get_bulk_repeat
===============

.. _poscar_utils_merge:

merge
=====

Merges two or more POSCAR files into one.
All files must have the same unit cell dimensions (within a tolerance defined by the ``--eps`` option).

The resultant POSCAR file will contain all atoms from all input files.
This can be used to stich together superstructures and a bulk cell, for example.

The utility raises an error if any two atoms are closer than ``--eps-collision``.
This can also be used to check if atoms from different slabs (with the same unit cell dimensions) are in the same positions.

**Usage**

.. code-block:: console

    $ viperleed poscar merge POSCAR1 POSCAR2 ... >POSCAR_OUT

**Options**

- ``--eps-cell``: tolerance for checking that the unit cell dimensions of the input files are the same (default: 1e-6)
- ``--eps-collision``: tolerance for checking that no two atoms are closer than this distance (default: 1e-6)

.. _poscar_utils_modify_vacuum:

modify_vacuum
=============

.. _poscar_utils_project_c_to_z:

project_c_to_z
==============

.. _poscar_utils_reorder_elements:

reorder_elements
================

.. _poscar_utils_rescale_cell:

rescale_cell
============


.. _poscar_utils_sort_by_z:

sort_by_z
=========

Sorts the atoms in the file by their z-coordinate within each element block.
To reorder the element blocks themselves, use :ref:`poscar_utils_reorder_elements`.

**Usage**

.. code-block:: console

    $ viperleed poscar sort_by_z <POSCAR_IN >POSCAR_OUT

**Options**

None



.. _poscar_utils_strip_comments:

strip_comments
==============

.. _poscar_utils_vasp_relax:

vasp_relax
==========

