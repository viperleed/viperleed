.. include:: /substitutions.rst

.. _conventions:

Conventions
===========

Coordinate system
-----------------

The lattice vectors describing the unit cell are named |a|, |b|, and |c|.
|a| and |b| must be in the :math:`xy` plane, and |c| must have a positive
:math:`z` component. The :math:`+z` direction is considered to be pointing
outwards from the surface (i.e., towards the vacuum). For :math:`(1 \times 1)`
cells, it is common (and highly recommended) to choose the |a| and |b| vectors
such that :math:`\mid \mathbf{a} \mid \leq \mid \mathbf{b} \mid`

Units
-----

ViPErLEED uses units of Ångström (Å) for all distances and vibrational
amplitudes in inputs and outputs. Some plots may display picometers
(1 pm = 0.01 Å), but this will **always** be clearly labeled.

Energies for input and output are in units of electronvolts (eV).
Parts of TensErLEED and certain phaseshift-generation programs may use
units of `Hartree <https://en.wikipedia.org/wiki/Hartree>`__ internally
(and in raw outputs).
ViPErLEED takes care of the conversion automatically.


Case sensitivity
----------------

Names of :ref:`input files<list_input_files>` are **case sensitive**, but
:ref:`parameter<parameters>` keywords and element names (see below) are not.

Beam Indices
------------

Diffraction beam orders are defined by the indices :math:`h` and
:math:`k`, relating the incident (\ :math:`\mathbf{k}`) and outgoing
wave vectors (\ :math:`\mathbf{k'}`) via

.. math::
    \mathbf{k'} = \mathbf{k} + \begin{pmatrix}h \\ k \end{pmatrix} \begin{pmatrix}\mathbf{a^*} \\ \mathbf{b^*}\end{pmatrix}

where :math:`\mathbf{a^*}` and :math:`\mathbf{b^*}` are the reciprocal-lattice
vectors. Beams are commonly labeled as :math:`(h|k)`, :math:`(h,k)` or
:math:`(hk)`. ViPErLEED generally uses the syntax ``(h | k)`` for machine
readability and to avoid ambiguity.

Comments
--------

ViPErLEED input files allow comments starting with``#`` and ``!``.
Anything on a line after either of these symbols will be ignored.
Indentation is allowed, but will be ignored.