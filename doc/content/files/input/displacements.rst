.. include:: /substitutions.rst

.. _displacements:

=============
DISPLACEMENTS
=============

The DISPLACEMENTS file defines the variations of geometry, vibrational
amplitudes and element concentrations that should be considered in the
search. In other words, the DISPLACEMENTS file defines the parameter
space for the search.

The file is split into three main blocks: Geometry, Vibrations, and
Occupations. The blocks are delimited by lines starting with an equals
sign "``=``", as in the following example

..  code-block:: none

   = GEO_DELTA
   O 1 z = -0.05 0.05 0.005      ! Oxygen atom 1 (and symmetry-equivalent atoms) will be displaced in z direction over the range [-0.05, 0.05] with step 0.005
   Ir 2-11 z = -0.05 0.05 0.01   ! Iridium atoms 2-11 (and symmetry-equivalent atoms) will be displaced in z direction over the range [-0.05, 0.05] with step 0.01
   O L(2-4) z = -0.05 0.05 0.01  ! Oxygen atoms in layers 2-4 will be displaced in z direction over the range [-0.05, 0.05] with step 0.01

   = VIB_DELTA
   O 1 = -0.05 0.05 0.02         ! Vibrational amplitude of oxygen atom 1 (and symmetry-equivalent atoms) will be varied over the range [-0.05, 0.05] with step 0.02
   Ir_top = -0.05 0.05 0.01      ! Vibrational amplitude of Iridium atoms in Ir_top sites will be varied over the range [-0.05, 0.05] with step 0.01

   = OCC_DELTA
   O 1 = O 0.8 1.0 0.05          ! Concentration of oxygen atom 1 (and symmetry-equivalent atoms) will be varied from 80% to 100% with 5% steps (rest: vacancies)

Indentation is allowed, but does not affect the function.

In each of the lines in the blocks, the left side of the ``=`` sign defines
which atoms are being assigned displacements, while the right side defines
the range of displacements to apply. Generally, atoms can be addressed by
:ref:`POSCAR<POSCAR>`  element, chemical element, or site type. If there is
an :ref:`element name collision<ElementNameCollision>`, the element will be
interpreted as the POSCAR element, so the assignment will be made for *all*
elements in :ref:`ELEMENT_MIX<ELSPLIT>`. Supplying a list of atom numbers
(**N** in :ref:`POSCAR<POSCAR>`, or by layer number using ``L(x)`` for layer
``x``) to further limit which atoms are being addressed is optional. If no
numbers are given, the displacements on the right are applied to all atoms
of the given element/site.

The exact syntax of the three blocks differs slightly, and is explained in
detail here:

-  :ref:`Geometrical displacements<GEODELTA>`
-  :ref:`Vibrational amplitudes<VIBDELTA>`
-  :ref:`Chemical substitution<OCCDELTA>`

Generally, any displacement applied to one atom will also be applied to
all symmetry-equivalent atoms (see Linking in :ref:`POSCAR<POSCAR>`),
such that the symmetry is preserved during the search
(eg. in-plane geometrical displacements will be mirrored for atoms linked
by a mirror symmetry plane). If multiple assignments are made for the same
atom, the assignments will be ignored if they are consistent, but the user
will be warned and the program may stop if there is a contradiction. For
example, if some of the iridium atoms 2-11 in the code above were linked
by a mirror plane, assigning a displacement in z direction to all of them
would be accepted, but assigning the very same in-plane displacement
(not parallel to the mirror plane) to all of them would lead to a
contradiction, as this would break the symmetry. In that case, it would
be preferable to assign the displacement explicitly to only *one* of the
symmetry-equivalent atoms.

.. note::
    See the :ref:`Domain calculations<domain_calculation>` page for information
    on how to format the DISPLACEMENTS file when optimizing multiple structures
    simultaneously.

.. toctree::
   :hidden:

   displacements/geodelta
   displacements/vibdelta
   displacements/occdelta


Advanced functionality
======================

.. note:: Indentation is allowed, but does not affect the function.

**Further constraints: Using CONSTRAIN blocks**

If some displacements should be constrained or linked in ways that go beyond
simple symmetry conservation, arbitrary displacements can be linked by using
:ref:`CONSTRAIN blocks<SEARCHCONSTRAINTS>`.

..  code-block:: none

   = VIB_DELTA
   Ir_top = -0.05 0.05 0.02     ! see above; for example below

   = CONSTRAIN
   geo O_top, Ir_top = linked   ! keep geometrical displacement index the same for O_top and Ir_top atoms. Requires the displacement ranges to have the same number of steps.
   vib Ir_def = linked          ! keep vibrational displacement index the same for all Ir_def atoms
   vib Ir_top = -0.03           ! although a displacement range is defined for Ir_top, fix its value to -0.03 instead
   vib Ir_top = ind(2)          ! same as the line before: Fix index to 2, i.e. the second entry in the displacement range

.. toctree::
   :hidden:

   displacements/searchconstraints
   displacements/symdelta

Running multiple searches
-------------------------

If you want to optimize multiple parameters not simultaneously, but end-to-end
(necessary e.g. for geometrical optimization), you can use multiple blocks in
the DISPLACEMENTS file to express this. After finishing one set of delta
calculations and search, the program will then loop back to execute delta
calculations and search again, starting from the optimized results of the
previous search.

For example, the following DISPLACEMENTS file would first optimize z position
and vibrational amplitudes simultaneously for the given set of atoms,
then run another search from the optimized z and vibrational amplitudes,
this time optimizing the x coordinate:

..  code-block:: none

   == SEARCH z

     = GEO_DELTA
     Ir L(1-6) z = -0.05 0.05 0.01

     = VIB_DELTA
     Ir = -0.005 0.005 0.0005

   == SEARCH x

     = GEO_DELTA
     Ir L(1-6) xy[1 0] = -0.03 0.03 0.01

The successive search blocks are each introduced by a line starting with
``== SEARCH``. Any text after the ``== SEARCH`` tag in the same line will
be treated as a 'name' for the block (in the above example ``z`` and ``x``).
The name will be referenced in the log files, but does not have any influence
on the behaviour of the program.


In-plane optimization shorthand
-------------------------------

If you want to run one search to optimize positions in one in-plane
direction (e.g., x), then another search for the other direction,
you can either write out two search blocks, or abbreviate by entering
just ``xy`` or ``ab`` as the direction, without the
:ref:`brackets to indicate a specific direction<GEODELTA>`. This will
effectively expand the block into two subsequent blocks, using the ``[1 0]``
direction in the first and the ``[0 1]`` direction in the second one.
For example

..  code-block:: none

   == SEARCH xy

     = GEO_DELTA
     Ir L(1-6) xy = -0.03 0.03 0.01

     = VIB_DELTA
     Ir = -0.005 0.005 0.0005

is equivalent to

..  code-block:: none

   == SEARCH xy[1 0]

     = GEO_DELTA
     Ir L(1-6) xy[1 0] = -0.03 0.03 0.01

     = VIB_DELTA
     Ir = -0.005 0.005 0.0005

   == SEARCH xy[0 1]

     = GEO_DELTA
     Ir L(1-6) xy[0 1] = -0.03 0.03 0.01

     = VIB_DELTA
     Ir = -0.005 0.005 0.0005


Looping searches
----------------

It is possible to have a set of search blocks running in a loop using the
tags ``<loop>`` and ``</loop>``. Note that these flags are expected to
always be directly before a ``== SEARCH`` statement, or before the end of
the file. Loops will end when their latest iteration does not yield a better
|R factor| than the previous iteration (note that this means each looped block
will be executed at least twice). For example, to optimize the z coordinate,
then search x based on the optimized z, and then loop back to searching z 
with the optimized x, the example above could be modified like this:

..  code-block:: none

   <loop>
   == SEARCH z

     = GEO_DELTA
     Ir L(1-6) z = -0.05 0.05 0.01

     = VIB_DELTA
     Ir = -0.005 0.005 0.0005

   == SEARCH x

     = GEO_DELTA
     Ir L(1-6) xy[1 0] = -0.03 0.03 0.01

   </loop>

You can also nest loops; for example, to optimize z, then loop in-plane 
optimization to get optimal x and y for that z, then start again with z, 
you could do:

..  code-block:: none

   <loop>

   == SEARCH z

     = GEO_DELTA
     Ir L(1-6) z = -0.05 0.05 0.01

     = VIB_DELTA
     Ir = -0.005 0.005 0.0005

   <loop>

   ! This block will automatically be split into separate blocks for optimizing x and y
   == SEARCH xy

     = GEO_DELTA
     Ir L(1-6) xy = -0.03 0.03 0.01

   </loop>
   </loop>
