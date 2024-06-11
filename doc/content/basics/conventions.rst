.. include:: /substitutions.rst

.. _conventions:

Conventions
===========

Coordinate system
-----------------

The lattice vectors describing the unit cell are named |a|, |b|, and |c|.
|a| and |b| must be in the :math:`(x, y)` plane, and |c| must have a positive
:math:`z` component. The :math:`+z` direction is considered to be pointing
outwards from the surface (i.e., towards vacuum). For (1 × 1) cells, it
is common (and highly recommended) to choose the |a| and |b| vectors such that
:math:`|\mathbf{a}| \leq |\mathbf{b}|`.


Units
-----

ViPErLEED uses units of ångström (Å) for all distances and vibrational
amplitudes in inputs and outputs. Some plots may display picometres
(1 pm = 0.01 Å), but this will **always** be clearly labeled.

Energies for input and output are in units of electronvolts (eV). Parts
of TensErLEED and certain phase-shift-generation programs use units of
`hartree <https://en.wikipedia.org/wiki/Hartree>`__ internally and in
the raw outputs. ViPErLEED takes care of the conversion automatically.


Case sensitivity
----------------

Names of :ref:`input files<list_input_files>` are **case sensitive**,
but :ref:`parameter<parameters>` keywords and
:ref:`element names<elementnamecollision>` are not.


Beam Indices
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


Comments
--------

ViPErLEED input files allow comments starting with ``#`` and ``!``.
Anything on a line after either of these symbols will be ignored.
Indentation is allowed, but will be ignored.