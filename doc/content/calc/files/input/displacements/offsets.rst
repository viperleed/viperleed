.. _offsets:

OFFSETS block
-------------

The offsets block is an optional part of the :ref:`DISPLACEMENTS` file.
If used, it has to be placed at the top of the DISPLACEMENTS file. The 
block is introduced by a header line of the form

..  code-block:: none

   = OFFSETS

This is followed by a list of offset declarations for atoms that should be
displaced away from their original
position in the reference calculation.

Note that offsets applied to parameters that are not varied during the search
will fix these parameters to the offset value.
The general syntax for a line in the OFFSETS block is (examples given below):

..  code-block:: none

   flag target[, target] = value

-  ``flag`` is either ``geo``, ``vib`` or ``occ``, depending what type of
   parameter should be addressed (compare to the other blocks in the
   :ref:`DISPLACEMENTS` file).
-  The syntax for addressing ``atoms`` is the same as in the other blocks, i.e.,
   by POSCAR element, chemical element, or site type, with optional further
   restriction by atom or layer number. Listing multiple sets of atoms in one
   line is also allowed using commas to separate them.
-  The ``value`` on the right is the offset value to be applied to the
   specified parameter(s). For geometric and vibrational displacements, this
   is a floating point number in Ångström. For occupation displacements, this
   is a chemical element followed by a floating point number between 0 and 1
   giving the occupation of the specified site(s).

.. todo:: Add example for occupation offsets
