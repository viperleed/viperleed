.. _occdelta:

================================================
Chemical substitution for structure optimization
================================================

Variations of chemical concentrations at given sites during the search
must be specified in a block starting with the

::

   = OCC_DELTA

header flag, followed by a list of *all* elements and the range of their
concentrations.
If the total site occupation does not sum to one at any point, this is
implicitly interpreted as a vacancy.
The command follows a slightly different syntax than those for the
:ref:`Geometrical displacements<GEODELTA>` or
:ref:`Vibrational amplitudes<VIBDELTA>`.


Syntax
======

:ref:`See also the example below<occdelta_example>`


::

   POSCAREL(_site) number(s) = chem1 start end step, chem2 start end step, chem3 start end step, chem4 start end step, chem5 start end step
   POSCAREL(_site) number(s) = fix1 chem1 + fix2 chem2 start end step, chem3 start end step, chem4 start end step, chem5 start end step

where

-  ``POSCAREL`` is the :ref:`POSCAR` element and **not** the chemical species 
   defined via :ref:`ELEMENT_MIX`.
   See also :ref:`element name collision<elementnamecollision>`.
-  ``site`` is optional, and has the same functionality as the one in
   the :ref:`Geometrical displacements<GEODELTA>`.
-  ``number(s)`` also behave as in the
   :ref:`Geometrical displacements<GEODELTA>`.

On the right-hand side of the ``=`` sign, blocks are included in a
comma-separated list. Each list element can have the two possible forms:

::

   chem start end step
   fix1 chem1 + fix2 chem2 (+ ...) start end step

where ``chem*`` is one of the **chemical** elements that you have defined
via :ref:`ELEMENT_MIX` or the special flag ``Vac`` (not case sensitive) for
vacancies. Notice that a **maximum of five distinct chemical elements**
(including vacancies) can be used at each atomic position.

``start``, ``end``, ``step``, and ``fix*`` are floating point numbers
between 0 and 1 defining the fractional occupations in percent.

When using the first syntax, the program will search all independent
combinations of ``chem*`` in the range specified.
The second syntax allows you to fix the relative concentration of (at
least two, at maximum five) elements via the numbers ``fix*``, creating
a "combined atom", whose concentration with respect to the other
chemical species specified in the following blocks will be varied
according to the ``start``, ``end`` and ``step`` values defined
afterwards.
When you combine chemical species with fixed concentration,
the ``fix*`` numbers will give the ratio between those chemical species,
while the overall concentration of the "combined atom" is still given
by ``start end step``.
The total number of chemical species on each site
must be **at most five**.


.. _occdelta_example:

Example
=======

::

   = OCC_DELTA
   ! Concentration of oxygen atom 1 (and symmetry-equivalent atoms) will
   ! be varied from 80% to 100% with 5% steps (rest: vacancies)
   O 1 = O 0.8 1.0 0.05

   ! Occupation of M_top sites (M is a POSCAR element) will be varied
   ! between 40% iron + 60% nickel and 60% iron + 40% nickel, with 5% steps.
   M_top = Fe 0.4 0.6 0.05, Ni 0.6 0.4 0.05

   ! As above, varying Fe 30-50%, Ni 60-40%, while at the same time
   ! keeping Ti constant at 10%
   M_top = Fe 0.3 0.5 0.05, Ni 0.6 0.4 0.05, Ti 0.1

How atoms are addressed on the left is described on the main
:ref:`DISPLACEMENTS` page. Note that if you want to vary the
concentrations of multiple elements, the *number of steps* in the ranges must
be the same, and none of the steps must have a total occupation greater than 1.
To keep an element at a constant concentration while varying others, you can
assign that concentration using only one number instead of three, as in the
example for Ti above. Alternatively, the same input format with of
``start end step`` with ``start == end`` and arbitrary ``step`` is
also interpreted as a constant concentration.

In the OCC_DELTA block, the element on the left *must* be the element as
defined in the :ref:`POSCAR` file, and the elements on the right
*must* be chemical elements, defined either by :ref:`ELEMENT_MIX`
or :ref:`ELEMENT_RENAME` in the :ref:`PARAMETERS` file.

Note that a **maximum of five distinct chemical elements** (including
vacancies) can be used at each atomic position.

For some applications, it can be useful to apply a static offset, without
re-doing the reference calculation. For this purpose, the OCC_DELTA block
also accepts single-value input (per element) on the right:

::

   = OCC_DELTA
   O 1 = O 0.8                         ! Concentration of oxygen atom 1 (and symmetry-equivalent atoms) will be fixed to 80% (rest: vacancies)
   M_top = Fe 0.6, Ni 0.4              ! Occupation of M_top sites (M is a POSCAR element) will be fixed to 60% iron + 40% nickel.


.. note::
   -  Due to the Fortran format currently used, ``start``, ``end``, and
      ``step`` will be truncated at the *second decimal digit* by
      rounding (i.e., 85.263 -> 85.26, while 85.265 -> 85.27).
   -  As for the :ref:`Geometrical displacements<GEODELTA>` and for the
      :ref:`Vibrational amplitudes<VIBDELTA>`, the concentration steps
      above will be applied to all symmetry-equivalent atoms, unless
      you turn off symmetry via :ref:`ISYM` or :ref:`SYMDELTA`.
   -  The **minimum** number of blocks is **one**.
      You can use this to specify a *fixed* chemical substitution on the
      atomic site, which can differ from the one you specified in the
      :ref:`VIBROCC`  file. This is generally not recommended.
