.. _conventions:

Conventions
===========

Coordinate system
-----------------

The lattice vectors describing the unit cell are named :math:`\vec{a}`, 
:math:`\vec{b}`, and :math:`\vec{c}`.
:math:`\vec{a}` and :math:`\vec{b}` must be in the :math:`xy` plane, 
and :math:`\vec{c}` must have a positive :math:`z` component.
The :math:`+z` direction is considered to be pointing outwards from the surface.

Units
-----

ViPErLEED uses units of Ångström (Å) for all distances and vibrational amplitudes in inputs, and outputs. 
Some plots may display picometers (pm), but this will **always** be clearly labeled.

Energies are in- and output in units of electron volts (eV).
Note that parts of TensErLEED and certain phaseshift generation scripts my use units of `Hartree <https://en.wikipedia.org/wiki/Hartree>`__ internally (and in raw outputs).
ViPErLEED takes care of the conversion automatically.

Case sensitivity
----------------

Names of :ref:`input files<list_input_files>` are **case sensitive**, but :ref:`parameter<parameters>` keywords are not.
