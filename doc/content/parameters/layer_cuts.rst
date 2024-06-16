.. include:: /substitutions.rst

.. _layer_cuts:

LAYER_CUTS
==========

LAYER_CUTS determines how the slab provided in the :ref:`POSCAR` file is 
divided in layers. Full dynamic multiple scattering is calculated inside 
each layer. The electron beams are then propagated as attenuated plane 
waves in between adjacent layers.

.. note::

   The :ref:`BULK_LIKE_BELOW` parameter offers an easy way to detect the bulk 
   repeat unit automatically, which will set both :ref:`N_BULK_LAYERS` and 
   LAYER_CUTS. However, be careful combining :ref:`BULK_LIKE_BELOW` with 
   LAYER_CUTS because the former will overwrite the latter during 
   initialization producing a fixed list of cut positions. Particular 
   care should be taken in adapting ``dc(float)`` and ``dz(float)``
   cut types after initialization.

**Allowed values:** ``list`` of floats (0 < ``x`` < 1),
``dz(float)``, or ``dc(float)``

**Default:** dz(1.2)

**Syntax:**

::

   LAYER_CUTS = 0.09 0.17 0.22 0.39 0.65    ! explicit list of c values
   LAYER_CUTS = dz(1.0)                     ! cut wherever distance along z is greater than 1.0 Å
   LAYER_CUTS = dc(1.2)                     ! cut wherever distance along c is greater than 1.2 Å
   LAYER_CUTS = 0.3 < dz(1.2) < 0.45 0.52   ! cut at c = 0.3, c = 0.45 and c = 0.52, and wherever distance along z is greater than 1.2 Å in range [0.3, 0.45]
   LAYER_CUTS = 0.15 0.3 < dz(1.2)          ! cut at c = 0.15, c = 0.3, and wherever distance along z is greater than 1.2 Å above c = 0.3

-  The ``list`` input explicitly contains the positions at which the layers
   should be cut (order does not matter). These positions are expressed as a
   fraction of the height of the cell as defined by the **c** vector in the
   POSCAR file. The list should not contain the edges of the cell (i.e., 0
   and 1). The topmost layer goes from *c* = the largest value in the list
   to *c* = 1.0.
-  ``dc(<float>)`` will automatically generate a list of cut positions, such
   that the distance along the **c** direction between the highest atom of
   each layer and the lowest atom of the layer above is larger than ``float``
   (in Ångström).
-  ``dz(<float>)`` works similarly to ``dc(<float>)``, but the distances are
   expressed along the *z* cartesian coordinate (in Ångström), i.e.,
   perpendicular to the plane of the **a** and **b** unit vectors.
-  Both ``dz`` and ``dc`` can be limited in scope by combining them with
   ``<`` or ``>`` characters, to specify that automatic cuts should be
   applied only above or below the values to the left and/or right. An
   example where this may be useful would be to specify the cuts in the
   bulk (which must not contain more than two layers) explicitly, but
   still get automatic cuts above.

.. note::

   -  When specifying a list or limiting the ``c`` range of ``dz``, take
      special care when a reference calculation is run from a previous
      optimized slab. Atoms may have moved substantially during the
      structural optimization (search), and can end up in the wrong
      layer if the LAYER_CUTS is given as a list of cut positions.
   -  The scattering calculations in the underlying TensErLEED package
      can become unstable the distances between layers are less than
      than approximately 1.0 Å. Thus, the program will output a *warning*
      in case the smallest interlayer distance resulting from LAYER_CUTS
      is smaller than 1.0 Å. To avoid this issue, it is advisable to use
      the ``dz(float)`` version of LAYER_CUTS, with ``float`` :math:`\geq` 1.
   -  If your system has only interlayer distances significantly smaller
      than 1.0 Å, the only viable option is treating the whole slab as a
      *single layer*. Notice that one should be **very careful** in handling
      this case, and the computing time will be high. In order to treat this,
      the :ref:`POSCAR` input should contain a very thick slab (at
      least three times the largest inelastic mean free path for the energy
      range used) such that all beams have negligible amplitude when 'exiting'
      the bottom surface. This should be constructed by manually adding layers
      with the bulk structure at the bottom of the slab. In addition, one
      should introduce at least one 'fake' bulk layer with an appropriate
      cut position (use the ``list`` version of LAYER_CUTS), since the code
      requires the presence of a 'bulk' section of the slab. If the thick
      slab is thick enough to prevent transmission of electrons towards this
      fake bulk layer, its presence should not affect the results. It is
      advisable to iteratively **test** the minimum thickness of the slab:
      run a reference calculation with a certain number *N* of true bulk
      layers + one fake bulk repeat; then a second reference calculation
      with *N*\ +1 true bulk layers + the fake bulk repeat. Then run an
      |R-factor| calculation between the theoretical beams of the two
      reference calculations to confirm that the difference in |R factor|
      is small enough for your purposes.
