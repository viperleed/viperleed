.. include:: /substitutions.rst

.. _geodelta:

Geometric displacements
=======================

Geometric displacements to be used during the search must be specified in a
block starting with the

..  code-block:: none

   = GEO_DELTA

header flag, followed by a list of displacements for each of the atoms that
one wants to move during the search.

.. note::
    Geometric displacements for any atom can only be applied along one
    :ref:`direction<geodelta_direction>`` (e.g. z, x, y, along arcs etc.) at a
    time. To optimize positions in multiple directions, multiple subsequent 
    search blocks are necessary.

**Example**:

..  code-block:: none

   = GEO_DELTA
   O 1 z = -0.05 0.05 0.005      ! Oxygen atom 1 (and symmetry-equivalent atoms) will be displaced in z direction over the range [-0.05, 0.05] with step 0.005
   Ir 2-11 z = -0.05 0.05 0.01   ! Iridium atoms 2-11 (and symmetry-equivalent atoms) will be displaced in z direction over the range [-0.05, 0.05] with step 0.01

   O 1 [1 0] = -0.05 0.05 0.005  ! Oxygen atom 1 will be displaced along a over the range [-0.05, 0.05] with step 0.005
   O 1 [0 1] offset = 0.02       ! Oxygen atom 1 displacements will additionally be offset by 0.02 A along b

The values on the right are interpreted as range of displacements in Ångström.
How atoms are addressed on the left is described on the main
:ref:`DISPLACEMENTS` page. However, The GEO_DELTA block contains 
one additional value on the left, which defines the *direction*
of the displacement, as described below.

When multiple searches are executed consecutively or in a loop, the
displacement ranges are centered around the optimized position from
previous searches per default. If you define an offset for an atom,
the optimized position from previous searches is *discarded* for this
atom, and the offset is applied to its original position. If you want
to center the displacement range around the original position of the
atom, you can also clear the offset without specifying a direction:

..  code-block:: none

   O 1 offset = 0         ! center around original atom position, instead of the position resulting from previous searches
   O 1 offset = clear     ! equivalent
   O 1 offset = original  ! equivalent

.. _geodelta_direction:

Direction
---------

Possible directions for displacements are:

.. Using a definition list here. Could also be a bullet list, not sure...

``z``
   Displacements along the direction orthogonal to the surface.
   Positive *z* values correspond to movements of the atoms away from the bulk.

   ..  code-block:: none

      O 1 z = -0.05 0.05 0.005      ! Oxygen atom 1 (and symmetry-equivalent atoms) will be displaced in z direction over the range [-0.05, 0.05] with step 0.005

``ab[n1 n2]`` (or just ``[n1 n2]``)
   In-plane displacements along the direction identified by the lattice
   vector with the whitespace-separated indices    ``n1`` (:math:`n_1`)
   and ``n2`` (:math:`n_2`). The first integer refers to the first vector
   |a| in the :ref:`POSCAR` file.
   The displacement direction will be positive in the direction of the
   vector :math:`\mathbf{v} = n_1 \mathbf{a} + n_2 \mathbf{b}`.
   ``n1`` and ``n2`` accept floating point values.
   Notice that the direction vector will be normalized, so
   ``[1 3]`` and ``[9 27]`` correspond to the same direction.

   ..  code-block:: none

      O 1 [1 0] = -0.05 0.05 0.005      ! Oxygen atom 1 will be displaced along a over the range [-0.05, 0.05] with step 0.005
      O 1 ab[1 -1] = -0.05 0.05 0.005   ! Oxygen atom 1 will be displaced diagonally along (a-b) over the range [-0.05, 0.05] with step 0.005
                                        ! Symmetry-equivalent atoms will be displaced such that the symmetry is preserved.


``xy[m1 m2]``

   In-plane displacements along the vector
   :math:`\mathbf{v}` = :math:`m1 \mathbf{x} + m2 \mathbf{y}`, where
   :math:`\mathbf{x}` and :math:`\mathbf{y}` are the unit vectors along
   the Cartesian axes. The direction vectors will be normalized.

   ..  code-block:: none

      O 1 xy[0 1] = -0.05 0.05 0.005    ! Oxygen atom 1 will be displaced along y over the range [-0.05, 0.05] with step 0.005
                                        ! Symmetry-equivalent atoms will be displaced such that the symmetry is preserved.

``azi(ab[c1 c2])``
   In-plane displacement around a *circular* trajectory centered at a
   specified point :math:`C`. The same convention as in the previous
   commands is used to specify the center:

   -  ``azi(ab[c1 c2])`` or just ``azi([c1 c2])`` means
      :math:`C = c_1 \mathbf{a} + c_2 \mathbf{b}` and
   -  ``azi(xy[c3 c4])`` means
      :math:`C = c_3 \mathbf{x} + c_4 \mathbf{y}`.

   The range on the right again defines a range of displacements in
   Ångström, in this case measured along the defined circular arc.
   Positive translations will translate to counterclockwise rotation as
   seen from vacuum. Zero displacement is the original position of the
   atom.
   Note that, since the displacement is given along the circular arc,
   the absolute displacement from the original position can be
   significantly smaller than for a linear displacement when the circle
   is small.

   ..  code-block:: none

      O 1 azi([0 0]) = -0.05 0.05 0.005   ! Oxygen atom 1 will be displaced along a circle centered on the origin by ±0.05 Å following the circular arc, with step 0.005
                                          ! Symmetry-equivalent atoms will be displaced such that the symmetry is preserved.

``r(ab[c1 c2])``
   In-plane displacement **r**\ elative (radial) to a specified point **C**.
   The same convention as in the previous commands is used to specify the
   point of reference, i.e.,

   -  ``r(ab[c1 c2])`` or ``r([c1 c2])`` means
      :math:`C = c_1 \mathbf{a} + c_2 \mathbf{b}` and
   -  ``r(xy[c3 c4])`` means :math:`C = c_3 \mathbf{x} + c_4 \mathbf{y}`.

   Positive values are interpreted as moving the atom *away* from point
   C, negative values move the atoms *towards* point C.

   ..  code-block:: none

      O 1 r([0 0]) = -0.05 0.05 0.005  ! Oxygen atom 1 will be displaced away from the origin over the range [-0.05, 0.05] with step 0.005
                                       ! Symmetry-equivalent atoms will be displaced such that the symmetry is preserved.

Offset
------

In addition to displacement along a specific direction, an offset along
a different direction can be defined. That offset will be added to the
"neutral" position of the atom, i.e. apply to **all** points in the
displacement range.

..  code-block:: none

   O 1 [0 1] offset = 0.02           ! Oxygen atom 1 displacements will be offset by 0.02 A along b
   O 1 ab[0 1] offset = 0.02         ! same as above
   O 1 xy[0 1] offset = 0.02         ! Oxygen atom 1 displacements will be offset by 0.02 A along y

Unlike the displacement ranges themselves, the offset flag allows
multiple assignment, as long as one of the assignments is in-plane
and the other one is out-of-plane:

..  code-block:: none

   ! in-plane and out-of-plane offsets can be combined:
   O 1 [0 1] offset = 0.02           ! Oxygen atom 1 displacements will be offset by 0.02 A along b ...
   O 1 z offset = 0.03               !       ... and by 0.03 A along z


.. note::
   -  If your input ``start``, ``stop``, and ``step`` values do not lead to
      an odd integer number of steps, the extremes of the interval will be
      extended in a symmetric fashion around the midpoint
      [= (``start``\ +\ ``stop``)/2] (i.e., ``step`` has precedence).
   -  Displacements of atoms will be **cross-checked for symmetry **
      **conservation** (unless you have turned off symmetry via
      :ref:`ISYM` and/or :ref:`SYMDELTA`),
      and the program will throw an **error** if inconsistencies arise.
      In general: atoms at *n*-fold rotational axes cannot be displaced;
      atoms on mirror planes can be moved only along the planes. You can
      find which displacement directions conserve the symmetry of your
      structure input in the comments added to the :ref:`POSCAR`
      file. Refer to the relation between plane groups in the
      :ref:`ISYM`  page in case you required a lowering of the 
      symmetry of your slab via :ref:`ISYM` or :ref:`SYMDELTA`.
   -  During one optimization run, an atom can only be displaced along
      **one** axis (so, for example, **not** sampling all in-plane directions
      at once). This is due to the way that the TensErLEED search is currently
      designed, with geometric displacements being optimized along a 1D array
      of points only. Since LEED is much more sensitive to variations of the
      out-of-plane geometry of your sample (small :math:`k_{\textrm{par}}`),
      it is a good idea to *first* run a few optimization runs on the *z*
      positions only, and treat in-plane displacements later as a refinement
      (unless your :ref:`POSCAR` model is *very far off* from the real 
      structure).
