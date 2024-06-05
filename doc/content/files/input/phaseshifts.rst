.. _phaseshifts:

PHASESHIFTS
===========

The ``PHASESHIFTS`` file is generated automatically with the eeasisss script
during initialization if needed. If a ``PHASESHIFTS`` file is already present,
a consistency check will be performed for the following features:

-   The energy range in the PHASESHIFTS file must be at least as large
    as the energy range to be calculated
    (:ref:`THEO_ENERGIES<theo_energies>`)
-   The number of blocks per energy must either match the number of
    chemical elements (i.e. elements in :ref:`POSCAR<POSCAR>` plus
    potentially elements added by :ref:`ELEMENT_MIX<ELSPLIT>`), or the
    number of distinct sites (see :ref:`SITE_DEF<SITEDEF>`) times the
    number of elements that can occupy any given site (i.e. the format
    generated automatically, :ref:`see below<phaseshift_format>`).
-   If the real part of the inner potentially is not defined explicitly,
    then the first line in the PHASESHIFTS file should contain the
    parameters defining that potential (see :ref:`V0_REAL<MUFTIN>`)

.. _phaseshift_format:

Format
______

The format of the ``PHASESHIFTS`` file is as follows:

**The first line** contains the number of blocks per energy, and four
parameters for the real part of the inner potential, and an optional timestamp.

The data below is listed by **energy**, where the energy is given as a
single floating-point value before the phase shifts for that energy.
Energies are given in Hartree.

For each energy, there is **one block per site and element occupying
that site**, where sites are defined by :ref:`SITE_DEF<SITEDEF>` and
elements are at least the elements from the :ref:`POSCAR file<POSCAR>`,
with additional entries if elements were added via the
:ref:`ELEMENT_MIX<ELSPLIT>`  parameter.

In each such block, there is **one floating point value per angular
momentum number** :math:`L`, representing the phase shift that an
electron with that :math:`L` experiences when it scatters at the given
site, occupied by the given element, at the given energy.

Example
_______

Below, you find an example of what a ``PHASESHIFTS`` file may look like. The
shown file is generated as part of the
:ref:`Ir(100)-O example<example_Ir(100)-O>`. It has of three block for
:ref:`elements<elementnamecollision>` (``O_surf``, ``Ir_surf``, ``Ir_def``),
each with phase shifts for angular momentum quantum number 0 through 17.
Note, that the phase shifts for a single element may be split across two
lines as is the case here.

::

    3   -11.89   -0.15  -88.82   17.74      Ir(100)-O_example 221130-104323
    0.5512
    1.4793-1.0892 0.0178 0.0006 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    -0.7672-0.3853-0.9502 0.0147 0.0007 0.0000 0.0000 0.0000 0.0000 0.0000
    0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    -0.7517-0.3790-0.9060 0.0152 0.0007 0.0000 0.0000 0.0000 0.0000 0.0000
    0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    0.6615
    1.3761-1.0791 0.0266 0.0010 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    -0.8554-0.4496-0.8502 0.0263 0.0014 0.0001 0.0000 0.0000 0.0000 0.0000
    0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    -0.8396-0.4421-0.8243 0.0271 0.0014 0.0001 0.0000 0.0000 0.0000 0.0000
    0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    0.7717
    1.2832-1.0943 0.0360 0.0015 0.0001 0.0000 0.0000 0.0000 0.0000 0.0000
    0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    -0.9590-0.5219-0.8529 0.0406 0.0024 0.0001 0.0000 0.0000 0.0000 0.0000
    0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    -0.9430-0.5132-0.8329 0.0421 0.0025 0.0001 0.0000 0.0000 0.0000 0.0000
    0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    0.8820
    ...
