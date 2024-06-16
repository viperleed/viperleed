.. _symdelta:

=================
The SYM_DELTA tag
=================

.. warning::
  SYM_DELTA is not well tested yet, and has at least one known bug.
  Use of :ref:`ISYM` instead of SYM_DELTA is encouraged where possible.

SYM_DELTA is a tag in the :ref:`DISPLACEMENTS` file that allows you to
temporarily change the symmetry of your slab, such that displacements
then assigned to atoms will be applied to fewer or no "symmetry-equivalent"
atoms. If you require lowering of symmetry, this will apply only to those
input lines between two SYM_DELTA lines (e.g., start a block with one line
turning off symmetry or changing the symmetry group, and end it with one
line turning the full symmetry back on or changing the symmetry to another
group).

Example
-------

..  code-block:: none

   = GEO_DELTA
   O 1 z = -0.05 0.05 0.005      ! Oxygen atom 1 (and symmetry-equivalent atoms) will be displaced in z direction over the range [-0.05, 0.05] with step 0.005
   SYM_DELTA = False
   Ir 1 3-5 z = -0.05 0.05 0.01  ! Iridium atoms 1 and 3-5 (but NOT their symmetry-equivalent atoms) will be displaced in z direction over the range [-0.05, 0.05] with step 0.01
   SYM_DELTA = True

   = VIB_DELTA
   Ir 1-6 = -0.05 0.05 0.02      ! Vibrational amplitude of iridium atoms 1-6 (and symmetry-equivalent atoms) will be varied over the range [-0.05, 0.05] with step 0.02


Acceptable values
-----------------

-  ``T``, ``True``, ``F``, ``False`` (not case sensitive): ``False`` turns off
   symmetry linking entirely, which is equivalent to setting the symmetry group
   to p1. ``True`` turns symmetry linking back on.
-  One can also directly specify a symmetry group with ``SYM_DELTA = group``.
   That group will then be used to restrict geometrical displacements or link
   symmetry-equivalent atoms. This is functionally equivalent to changing the
   value of :ref:`ISYM`, but only for the operations that follow. Note that
   only symmetry *reduction* from the overall slab symmetry is allowed. See
   :ref:`ISYM`  for a more detailed explanation of allowed symmetry changes.

The use of SYM_DELTA should be reserved for **highly specific** cases and is
**generally discouraged**. If the desired effect can be achieved by lowering
the *overall* symmetry, then using the :ref:`ISYM` tag of the :ref:`PARAMETERS`
file is always preferable.
