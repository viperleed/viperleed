.. _aseapi:

Atomic Simulation Environment Interface
=======================================

The tleedm package of ViPErLEED supports running directly from surface structures from the `Atomic Simulation Environment <https://wiki.fysik.dtu.dk/ase/>`__ (ASE).
ASE is a widely used python library for setting up and controlling atomistic simulations such as density functional theory calculations and molecular dynamics simulations.
Via the ASE interface provided in tleedm, an ASE atoms object can be used as surface model for a LEED I(V) calculation.

This can be used to enable automated runs and sample surface structures 
far beyond what is possible in the Tensor LEED approach :cite:p:`rousTheoryTensorLEED1989,rousTensorLEEDTechnique1986`.
Surface structures passed to tleedm via the application programming 
interface (API) need to fulfill the same conventions (e.g. vacuum side of surface towards +\ **z** direction) as applicable for POSCAR files.

The API provides a number of python function that allow calling and 
starting LEED I(V) calculations.
For an example of how to use the ASE interface, see :ref:`this tutorial<example_Ni(110)>`.

.. note:: 
    The ViPErLEED ASE interface is an area of active development and may change in the future.
    If you have any suggestions for features or feedback concerning the ASE interface, feel free to contact the ViPErLEED team (preferably by opening an issue on the ViPErLEED Github page).

**TODO: Docstring from viperleed_from_ase can be put here after package reorganization**

.. _aseapi_auto_sites:

Automatic site definitions
--------------------------

To facilitate batch processing of structures, the ASE interface allows for a simplistic automated assignment of site definitions (usually done via the :ref:`parameter SITE_DEF<sitedef>`).

If the SITE_DEF parameter is not defined in a ViPErLEED run using the ASE interface, ViPErLEED will try to assign "surface" sites on its own.
For this calculation, every atom is considered as a solid sphere with a radius proportional to the elements' covalent radius.
Then, going from highest to lowest atom (:math:`z` position), ViPErLEED checks if an atom is "visible" from vacuum, or if the line of sight is blocked by higher-up atoms.
Atoms that are "visible" are declared as surface atoms (site ``surf``), while all other atoms will be given the default site type (``def``).
Note, this will create, at maximum, two site types per species.
