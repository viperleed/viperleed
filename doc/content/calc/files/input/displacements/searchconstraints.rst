.. _searchconstraints:

Constraints for structure optimization
======================================

This block of the :ref:`DISPLACEMENTS` file allows you to define specific 
constraints other than symmetry for the structure-optimization routine. The 
commands are introduced by a header line of the form

..  code-block:: none

   = CONSTRAIN

The general syntax for a line in the CONSTRAIN block is (examples given below):

..  code-block:: none

   flag atoms[, more atoms] = value

-  ``flag`` is either ``geo``, ``vib`` or ``occ``, depending what type of
   parameter should be addressed (compare to the other blocks in the
   :ref:`DISPLACEMENTS` file).
-  The syntax for addessing ``atoms`` is the same as in the other blocks, i.e.,
   by POSCAR element, chemical element, or site type, with optional further
   restriction by atom or layer number. Listing multiple sets of atoms in one
   line is also allowed using commas to separate them.
-  The ``value`` on the right is either ``linking`` to link parameters to
   each other, a floating point number to lock in a specific value from the
   displacements range, or ``ind(N)``, where N is an integer (counting from 1)
   giving an index in the displacement range that should be locked. See below
   for examples for each of these cases.

Linking parameters
------------------

The main function of CONSTRAIN blocks is to **link** parameters together, i.e.,
if one of them is at the *N*\ th index of its displacement range, all of them
will be. A simple example would be to keep all vibration amplitudes for a
given site type the same:

::

   = VIB_DELTA
   Ir_top = -0.05 0.05 0.01       ! vary vibration amplitudes for all atoms in Ir_top sites over the range [-0.05, 0.05] with step 0.01
   Ir_def = -0.03 0.03 0.01       ! vary vibration amplitudes for all atoms in Ir_def sites over the range [-0.03, 0.03] with step 0.01

   = CONSTRAIN
   vib Ir_top = linked            ! keep the vibration amplitudes of all Ir_top atoms the same at all times
   vib Ir_def = linked            ! keep the vibration amplitudes of all Ir_def atoms the same at all times

Note that for linking to work, the linked parameters have to have the same
number of steps. In the above example, the following would **not** work,
because the displacement ranges have different sizes:

::

   = CONSTRAIN
   vib Ir = linked                ! keep the parameter index for the vibration amplitudes of all Ir atoms the same at all times

However, it can be made to work if the displacement
ranges have the same size, i.e.,

::

   = VIB_DELTA
   Ir_top = -0.05 0.05 0.01     ! vary vibration amplitudes for all atoms in Ir_top sites over the range [-0.05, 0.05] with step 0.01
   Ir_def = -0.025 0.025 0.005  ! vary vibration amplitudes for all atoms in Ir_def sites over the range [-0.025, 0.025] with step 0.005

   = CONSTRAIN
   vib Ir = linked              ! keep the vibration amplitude index of all Ir_top atoms the same at all times

Note that while in this example, the *indices* for vibration variation
of Ir_top and Ir_def are tied together, that does in no way mean that the
*values* have to be the same. Both the *base values* (defined in the
:ref:`VIBROCC` file) and the actual values within the displacement
ranges can differ.

You can also force atoms from different sites and/or elements to move in a
concerted fashion by linking the respective parameters:

::

   = GEO_DELTA
   O L(1-2) z = -0.05 0.05 0.005     ! move oxygen atoms in layers 1 and 2 along z over range [-0.05, 0.05] with step 0.005
   Ir L(1) z = 0.03 -0.03 0.003      ! move iridium atoms in layer 1 along z over range [0.03, -0.03] with step 0.003

   = CONSTRAIN
   geo O L(1-2), Ir L(1) = linked    ! link the indices for the displacement ranges above

In this example, all oxygen atoms in layers 1 and 2 would move in a concerted
fashion, as would the iridium atoms in layer 1. Furthermore, by linking them
*in the same line*, the movements of oxygen and iridium are coupled as given
by the two ranges, so when oxygen moves up, iridium moves down (e.g., when O
moves to +0.01 A, Ir moves to -0.006 A).

Freezing parameters
-------------------

Apart from linking parameters, CONSTRAIN blocks can also be used to freeze a
parameter to a specific value. This can be used to enforce an offset for a
given value (within the limitations of tensor LEED) without having to re-do
the reference calculation. It is also equivalent to giving a displacement
range with only a single value.

Freezing parameters can be achieved by specifying a value from the
displacements range, or by giving a specific index:

::

   = VIB_DELTA
   Ir_top = -0.05 0.05 0.02     ! vary vibration amplitudes for Ir_top atoms over the range [-0.05, 0.05] with step 0.02

   = CONSTRAIN
   vib Ir_top = -0.03           ! although a displacement range is defined for Ir_top, fix it's value to -0.03 instead
   vib Ir_top = ind(2)          ! same as the line before: Fix index to 2, i.e. the second entry in the displacement range

Note that in the ind(N) function, indices are counted starting at 1,
not at 0, to keep them consistent with values in the SD.TL file.
