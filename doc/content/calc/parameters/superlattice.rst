.. _superlattice:

SUPERLATTICE
============

SUPERLATTICE defines the relationship between the real-space unit vectors
of the bulk and those of the surface.

**Default:** Detect automatically from bulk. See below for details.

**Syntax examples:**

::

   SUPERLATTICE = c(2x2)
   SUPERLATTICE = p(2x1)
   SUPERLATTICE = (2x1)                    # same as p(2x1)
   SUPERLATTICE = (sqrt(2) x sqrt(2))R45   # same as M = 1 1, -1 1
   SUPERLATTICE = (2*sqrt(2) x sqrt(2))R45
   SUPERLATTICE M = 1 1, -1 1              # same as (sqrt(2) x sqrt(2))R45

If SUPERLATTICE is not defined, ViPErLEED will perform a search for possible
unit cell vectors for the bulk (as defined by :ref:`LAYER_CUTS` and
:ref:`N_BULK_LAYERS`), and choose two vectors that minimize the unit
cell area and circumference. The resulting SUPERLATTICE matrix will
then be further optimized to bring the bulk unit cell to its highest symmetry
form (preserving the area). Note that for many systems, multiple choices of the
bulk unit cell are possible, but result in different SUPERLATTICE definitions;
this choice will be made arbitrarily if SUPERLATTICE is not actively defined.

Even when SUPERLATTICE is defined, an automatic search for the best bulk unit
cell will run. If a unit cell with a smaller area is found, the program will
warn and stop. If required, halting can be suppressed using the
:ref:`HALTING` parameter.

.. note::

   -  For Wood notation, use only positive float values. The separator
      in the brackets should be the letter 'x', as in the examples.
   -  In matrix notation use the convention ``'M = m11 m12, m21 m22``', with
      **a**\ :sub:`surf` = m\ :sub:`11`\ ·\ **a**\ :sub:`bulk`
      + m\ :sub:`12`\ ·\ **b**\ :sub:`bulk`,
      and **b**\ :sub:`surf` = m\ :sub:`21`\ ·\ **a**\ :sub:`bulk`
      + m\ :sub:`22`\ ·\ **b**\ :sub:`bulk`, i.e., the bulk row vectors
      **a**\ :sub:`bulk` and **b**\ :sub:`bulk` get multiplied with **M**
      from the left to obtain the surface row vectors **a**\ :sub:`surf`
      and **b**\ :sub:`surf` respectively. In matrix notation all elements
      should be signed **integers** (floats will also be accepted, but this
      is not recommended).

.. note::
   Make sure to get the relation between the structure in the :ref:`POSCAR`
   file and the superlattice right. For example, in c(4x8)  "4" must refer
   to the **a** and "8" to the **b** lattice vector of the POSCAR, not reverse,
   unless :ref:`SYMMETRY_CELL_TRANSFORM` specifies that the superstructure cell
   differs from the cell given in the POSCAR file.

.. todo::
   Is it correct like this? Then one should mention that the POSCAR must
   contain exactly one superstructure cell, unless SYMMETRY_CELL_TRANSFORM
   defines it otherwise (also in the POSCAR description)! **TODO** which
   POSCAR? The user-supplied one or the one after symmetry detection??? -ms


When using Wood notation, the regular expression interpreting the input is:

::

   (?P<type>[PCpc]*)\s*\(\s*(?P<g1>.+)\s*[xX]\s*(?P<g2>.+)\s*\)\s*[rR]*\s*(?P<alpha>[\d.]*)

The expressions ``g1`` and ``g2`` in (g1 x g2) are interpreted afterwards;
here, only sqrt(value) and multiplicative expressions are parsed, so
``'2*sqrt(2)``' will be read correctly, but e.g. ``'(1/3)*sqrt(3)``'
will yield an error. For these cases, using a matrix is usually better.
