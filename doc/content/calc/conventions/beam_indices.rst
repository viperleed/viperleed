
.. _conventions_beam_indices:

Beam indices
------------

Diffraction-beam orders are defined by the indices :math:`h` and
:math:`k`, relating the incident (\ :math:`\mathbf{k}`) and outgoing
wave vectors (\ :math:`\mathbf{k'}`) via

.. math::
    \mathbf{k'} = \mathbf{k} + (h, k) \begin{pmatrix}\mathbf{a}^*_\mathrm{bulk} \\ \mathbf{b}^*_\mathrm{bulk} \end{pmatrix} ,

where :math:`\mathbf{a}^*_\mathrm{bulk}` and :math:`\mathbf{b}^*_\mathrm{bulk}`
are the reciprocal-lattice vectors of the (1 × 1) cell of the bulk.
Notice that beam indices are always defined by the reciprocal unit cell of
the bulk, irrespective of the periodicity of the superstructure. Beams
originating exclusively from the superstructure have fractional :math:`h`
or :math:`k` indices. It is common practice to label beams as :math:`(h|k)`,
:math:`(h, k)`, or :math:`(hk)`. ViPErLEED uses the ``(h | k)`` syntax for
machine readability and to avoid ambiguity. Spaces are also accepted as
separators in the inputs.
