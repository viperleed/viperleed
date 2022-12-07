.. _phaseshifts:

The PHASESHIFTS file
====================

The PHASESHIFTS file is generated automatically with the EEASiSSS.x script during initialization if needed. If a PHASESHIFTS file is already present, a consistency check will be performed for the following features:

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

The format of the automatically generated PHASESHIFT file is as follows:

**The first line** contains the number of blocks per energy, and four 
parameters for the real part of the inner potential.

The data below is listed by **energy**, where the energy is given as a 
single floating-point value before the phaseshifts for that energy. 
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

**TODO** Should we give a sample file like we do in :ref:`Expbeams` ?