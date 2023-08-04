.. _poscar_utils:

POSCAR utilities
################

All POSCAR utilities support the following options:

- ``-h``, ``--help``: print help message with the available options
- ``-v``, ``--verbose``: increase output verbosity

.. hint::

    Most ViPErLEED POSCAR utilities read from ``stdin`` and write to ``stdout``.
    This means that you can easily use them in a pipeline, e.g.:

    .. code-block:: console

        $ cat POSCAR | viperleed poscar delete_above C | viperleed poscar enforce_symmetry >POSCAR_OUT


.. tip::

    If you find these utilities useful, consider adding an alias to your ``.bashrc`` so you don't have to type ``viperleed poscar`` every time.

.. _poscar_utils_attach_bulk:

attach_bulk
===========

.. _poscar_utils_delete_above:

delete_above
============

Deletes all atoms in the POSCAR file above the specified fraction of the :math:`\vec{c}` vector.

With the ``--verbose`` option, the utility prints the number of atoms deleted.
This can also be useful to quickly check the number of atoms above a certain height.

**Usage**

.. code-block:: console

    $ viperleed poscar delete_above C <POSCAR_IN >POSCAR_OUT

**Additional Options**

- ``c`` (required): the fraction of the :math:`\vec{c}` above which to delete atoms


.. _poscar_utils_delete_below:

delete_below
============

Same as :ref:`poscar_utils_delete_above`, but deletes all atoms below the specified fraction of the :math:`\vec{c}` vector.

**Usage**

.. code-block:: console

    $ viperleed poscar delete_below C <POSCAR_IN >POSCAR_OUT

**Options**

- ``c`` (required): the fraction of the :math:`\vec{c}` below which to delete atoms

.. _poscar_utils_delete_between:

delete_between
==============

Same as :ref:`poscar_utils_delete_above` and :ref:`poscar_utils_delete_below`, but deletes all atoms between the specified fractions of the :math:`\vec{c}` vector.

**Usage**

.. code-block:: console

    $ viperleed poscar delete_between C1 C2 <POSCAR_IN >POSCAR_OUT

**Additional Options**

- ``c1`` (required): delete atoms with :math:`c_1 < c < c_2`
- ``c2`` (required): see above


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

**Additional Options**

- ``--eps-cell``: tolerance for checking that the unit cell dimensions of the input files are the same (default: 1e-6)
- ``--eps-collision``: tolerance for checking that no two atoms are closer than this distance (default: 1e-3)

.. _poscar_utils_modify_vacuum:

modify_vacuum
=============

.. _poscar_utils_project_c_to_z:

project_c_to_z
==============

.. _poscar_utils_reorder_elements:

reorder_elements
================

Reorders the element blocks in the POSCAR file.

By default, the element blocks are reordered by ascending atomic number.
Use options listed below to change the order.

**Usage**

.. code-block:: console

    $ viperleed poscar reorder_elements <POSCAR_IN >POSCAR_OUT               # ascending atomic number
    $ viperleed poscar reorder_elements --custom=O,Fe <POSCAR_IN >POSCAR_OUT # custom order

**Additional Options**

- ``--alphabetical``: sort elements by alphabetical order of the element symbols
- ``--descending``: sort elements by descending atomic number
- ``--custom``: sort elements by a custom order (comma-separated list of element symbols)

.. _poscar_utils_rescale_cell:

rescale_cell
============

Rescales the unit cell dimensions of the POSCAR file by the specified factor.

**Usage**

.. code-block:: console

    $ viperleed poscar rescale_cell 1.01 <POSCAR_IN >POSCAR_OUT           # stretch isotropically by 1%
    $ viperleed poscar rescale_cell 1.01 1.02 0.99 <POSCAR_IN >POSCAR_OUT # stretch anisotropically

**Additional Options**

- ``scaling``: (required) One or three scaling factors for the unit cell.
  If three values are given, the scaling factors are applied to the a, b, and c vector, respectively.
  If only one value is given, an isotropic scaling is applied.

.. _poscar_utils_sort_by_z:

sort_by_z
=========

Sorts the atoms in the file by their z-coordinate within each element block.
To reorder the element blocks themselves, use :ref:`poscar_utils_reorder_elements`.

**Usage**

.. code-block:: console

    $ viperleed poscar sort_by_z <POSCAR_IN >POSCAR_OUT

**Additional Options**

None


.. _poscar_utils_strip_comments:

strip_comments
==============

Strips all comments from the POSCAR file (e.g. :ref:`SITE_DEFs<sitedef>` information added by ViPErLEED).

**Usage**

.. code-block:: console

    $ viperleed poscar strip_comments <POSCAR_IN >POSCAR_OUT

**Additional Options**

None

.. _poscar_utils_vasp_relax:

vasp_relax
==========

