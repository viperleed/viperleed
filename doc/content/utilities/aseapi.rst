.. _aseapi:

Atomic Simulation Environment API
=================================

The tleedm package of ViPErLEED supports running directly from surface 
structures from the `Atomic Simulation Environment 
<https://wiki.fysik.dtu.dk/ase/>`__ (ASE).
ASE is a widely used python 
library for setting up and controlling atomistic simulations such as 
density functional theory calculations and molecular dynamics s
imulations.
Via the ASE interface provided in tleedm, an ASE atoms object can be 
used as surface model for a LEED I(V) calculation.

This can be used to enable automated runs and sample surface structures 
far beyond what is possible in the Tensor LEED approach [1].
Surface structures passed to tleedm via the application programming 
interface (API) need to fulfill the same conventions (e.g. vacuum side o
f surface towards +\ **z** direction) as applicable for POSCAR files.

**TODO: should we provide an example file and then link to it?**

The API provides a number of python function that allow calling and 
starting LEED I(V) calculations.

**TODO: copy docstring for run_from_ase?**

| [1] Rous, P. J.; Pendry, J. B. *The Theory of Tensor LEED*. `Surface Science 1989, 219 (3), 355–372 <https://doi.org/10.1016/0039-6028(89)90513-X>`__.
