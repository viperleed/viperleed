.. _poscar_utils:

.. include:: /substitutions.rst

POSCAR utilities
################

ViPErLEED comes with a set of utilities for easy manipulation of POSCAR files.
These utilities can be helpful for creating and modifying input structures for
ViPErLEED calculations.

All POSCAR utilities support the following options:

- ``-h, --help``: show a help message with the available options
- ``-v, --verbose``: increase output verbosity

.. hint::

    Most ViPErLEED POSCAR utilities read from ``stdin`` and write to
    ``stdout``. This means that you can easily use them in a pipeline,
    e.g.,

    .. code-block:: console

        $ cat POSCAR | viperleed poscar delete_above C | viperleed poscar enforce_symmetry >POSCAR_OUT


.. tip::

    If you find these utilities useful, consider adding an alias to your
    ``.bashrc`` so you don't have to type ``viperleed poscar`` every time.

.. _poscar_utils_attach_bulk:

attach_bulk
===========

Interactive script that takes a slab POSCAR and adds a bulk POSCAR on the
bottom, expanding the the unit cell in the :math:`\vec{c}` direction
accordingly.

**Usage**

.. code-block:: console

    $ viperleed poscar attach_bulk

**Additional Options**

None.


.. _poscar_utils_delete_above:

delete_above
============

Deletes all atoms in the POSCAR file above the specified fraction of the
:math:`\vec{c}` vector.

With the ``--verbose`` option, the utility prints the number of atoms deleted.
This can also be useful to quickly check the number of atoms above a certain
height.

**Usage**

.. code-block:: console

    $ viperleed poscar delete_above C <POSCAR_IN >POSCAR_OUT

**Additional Options**

- ``c`` (required): the fraction (floating point number) of :math:`\vec{c}`
  above which to delete atoms.


.. _poscar_utils_delete_below:

delete_below
============

Same as :ref:`poscar_utils_delete_above`, but deletes all atoms below the
specified fraction of the :math:`\vec{c}` vector.

**Usage**

.. code-block:: console

    $ viperleed poscar delete_below C <POSCAR_IN >POSCAR_OUT

**Options**

- ``c`` (required): the fraction (floating point number) of :math:`\vec{c}`
  below which to delete atoms.

.. _poscar_utils_delete_between:

delete_between
==============

Same as :ref:`poscar_utils_delete_above` and :ref:`poscar_utils_delete_below`,
but deletes all atoms between the specified fractions of the :math:`\vec{c}`
vector.

**Usage**

.. code-block:: console

    $ viperleed poscar delete_between C1 C2 <POSCAR_IN >POSCAR_OUT

**Additional Options**

- ``c1`` (required): delete atoms with :math:`c_1 < c < c_2`
- ``c2`` (required): see above


.. _poscar_utils_enforce_symmetry:

enforce_symmetry
================

Finds the planegroup of the POSCAR file and enforces it by moving atoms to
symmetric positions.

Symmetry detection works the same as the
:ref:`find_symmetry<poscar_utils_find_symmetry>`
utility but here a symmetrized POSCAR file is returned.

**Usage**

.. code-block:: console

    $ viperleed poscar enforce_symmetry <POSCAR_IN >POSCAR_OUT

**Additional Options**

- ``-e, --symmetry-eps``: Epsilon for in-plane symmetry detection in Å. Behaves
  like :ref:`sym_eps` in the :ref:`PARAMETERS<parameters>` file. Default: 0.1Å
- ``--symmetry-eps-z``: Epsilon for out-of-plane symmetry detection in Å.
  Behaves like the second argument of :ref:`sym_eps` in the
  :ref:`PARAMETERS<parameters>` file. If not provided, the value of
  ``--symmetry-eps`` is used.
- ``--planegroup``: Planegroup to enforce.
  Default: detected automatically from the slab.
  Use this option to override the automatic detection and
  manually lower the symmetry.

**Example**

.. code-block:: console

    $ viperleed poscar enforce_symmetry <POSCAR_IN >POSCAR_OUT --symmetry-eps 0.01

.. _poscar_utils_find_symmetry:

find_symmetry
=============

Finds the planegroup of the POSCAR file and prints it to ``stdout``.
This utility uses the same algorithm for symmetry detection as is used
in ViPErLEED calculations.

**Usage**

.. code-block:: console

    $ viperleed poscar find_symmetry <POSCAR_IN

**Additional Options**

- ``-e, --symmetry-eps``: Epsilon for in-plane symmetry detection in Å.
  Behaves like :ref:`sym_eps` in the :ref:`PARAMETERS<parameters>` file.
  Default: 0.1Å
- ``--symmetry-eps-z``: Epsilon for out-of-plane symmetry detection in Å. .
  Behaves like the second argument of :ref:`sym_eps` in the
  :ref:`PARAMETERS<parameters>` file. If not provided, the value of
  ``--symmetry-eps`` is used.


.. _poscar_utils_get_bulk_repeat:

get_bulk_repeat
===============

Interactive script that reads a POSCAR file, asks at what c value the bulk
starts, then automatically reduces the size of the POSCAR to non-redundant
bulk layers only, and outputs the appropriate
:ref:`N_BULK_LAYERS<n_bulk_layers>` and :ref:`BULK_REPEAT` values.
Additionally, the files ``POSCAR_bulk`` containing the bulk unit-cell and
a file ``POSCAR_min`` containing the minimal surface slab will be written.

**Usage**

.. code-block:: console

    $ viperleed poscar get_bulk_repeat

**Additional Options**

None.

.. _poscar_utils_merge:

merge
=====

Merges two or more POSCAR files into one.
All files must have the same unit cell dimensions (within a tolerance defined
by the ``--eps`` option).

The resultant POSCAR file will contain all atoms from all input files. This
can be used to stitch together superstructures and a bulk cell, for example.

The utility raises an error if any two atoms are closer than
``--eps-collision``. This can also be used to check if atoms from different
slabs (with the same unit cell dimensions) are in the same positions.

**Usage**

.. code-block:: console

    $ viperleed poscar merge POSCAR1 POSCAR2 ... >POSCAR_OUT

**Additional Options**

- ``--eps-cell``: tolerance for checking that the unit cell dimensions of the
  input files are the same (default: 1e-1)
- ``--eps-collision``: tolerance for checking that no two atoms are closer than
  this distance (default: 0.1)

.. _poscar_utils_modify_vacuum:

modify_vacuum
=============

Modifies the vacuum spacing of a POSCAR file.

While most :term:`DFT` codes use periodic boundary conditions along the z
direction, in LEED calculations the symmetry has to be broken in order to
simulate a surface. This utility allows to modify the vacuum spacing of a
POSCAR file by adding or removing vacuum around the slab.

**Usage**

.. code-block:: console

    $ viperleed poscar modify_vacuum 10 <POSCAR_IN >POSCAR_OUT # add 10 Å of vacuum

**Additional Options**

- ``vacuum`` (required): Add or remove this amount of vacuum in Å.
  If the flag ``--absolute`` is set, the total vacuum spacing (measured
  from topmost to bottommost atom) will be set to this value.
- ``--absolute``: see above.

.. _poscar_utils_project_c_to_z:

project_c_to_z
==============

Projects the :math:`\vec{c}` vector of the POSCAR file onto the
:math:`\vec{z}` axis. Note this does not alter atomic coordinates,
only the orientation of the lattice vectors. The bulk-stacking
direction is assumed to be along the :math:`\vec{z}` vector.
See also the :ref:`page on used conventions<conventions>`.

**Usage**

.. code-block:: console

    $ viperleed poscar project_c_to_z <POSCAR_IN >POSCAR_OUT

**Additional Options**

None.

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

- ``--alphabetical``: sort elements by alphabetical order of the element
  symbols
- ``--descending``: sort elements by descending atomic number
- ``--custom``: sort elements by a custom order (comma-separated list of
  element symbols)

.. _poscar_utils_rescale_cell:

rescale_cell
============

Rescales the unit cell dimensions of the POSCAR file by the specified factor.
Irrespective of how many values are give, this utility will directly alter the
unit cell basis vectors in the POSCAR file, not the scaling factor (line 2).

**Usage**

.. code-block:: console

    $ viperleed poscar rescale_cell 1.01 <POSCAR_IN >POSCAR_OUT           # stretch isotropically by 1%
    $ viperleed poscar rescale_cell 1.01 1.02 0.99 <POSCAR_IN >POSCAR_OUT # stretch anisotropically

**Additional Options**

- ``scaling``: (required) One or three scaling factors for the unit cell.
  If three values are given, the scaling factors are applied to the
  :math:`\vec{a}`, :math:`\vec{b}`, and :math:`\vec{c}` vector, respectively.
  If only one value is given, an isotropic scaling is applied.

.. _poscar_utils_sort_by_z:

sort_by_z
=========

Sorts the atoms in the file by their z-coordinate within each
element block. To reorder the element blocks themselves, use
:ref:`poscar_utils_reorder_elements`.

**Usage**

.. code-block:: console

    $ viperleed poscar sort_by_z <POSCAR_IN >POSCAR_OUT

**Additional Options**

- ``--reversed``: sort elements bottom to top (default: top to bottom)


.. _poscar_utils_strip_comments:

strip_comments
==============

Strips all comments from the POSCAR file (e.g. :ref:`SITE_DEF<sitedef>`
information added by ViPErLEED). This can also be used to strip ion
velocities from a VASP POSCAR file.

**Usage**

.. code-block:: console

    $ viperleed poscar strip_comments <POSCAR_IN >POSCAR_OUT

**Additional Options**

None

.. _poscar_utils_vasp_relax:

vasp_relax
==========

Formats the POSCAR file for use with :term:`VASP`.

It can often be useful to "pre-relax" a surface structure with :term:`DFT`
calculations before performing a |LEED-IV| analysis. This utilities facilitates
this by formatting the POSCAR file for relaxation with :term:`VASP`.
The vasp_relax utility adds the following information to the POSCAR file:

- the tag ``Selective dynamics``, which indicates to VASP that selected ion
  positions are allowed to move
- three boolean flags (`T`, `F`) for each atom indicating whether the atom
  is  allowed to move along the :math:`\vec{a}`, :math:`\vec{b}`, and
  :math:`\vec{c}` unit cell vectors, respectively

In general, it can be useful to optimize the positions of the topmost layers
of atoms, while keeping the positions of the atoms in the bulk fixed.
The ``above_c`` value should be chosen such that bulk atoms are not allowed
to move to prevent the bulk lattice parameters from changing.

**Usage**

.. code-block:: console

    $ viperleed poscar vasp_relax 0.20 <POSCAR_IN >POSCAR_OUT
    $ viperleed poscar vasp_relax 0.35 --all_directions <POSCAR_IN >POSCAR_OUT

**Additional Options**

- ``above_c``: (required) the fraction of the :math:`\vec{c}` vector above
  which to allow atoms to move
- ``--all_directions``: allow all atoms to move along all three unit cell
  vectors (default: only allow movement along :math:`\vec{c}`)
