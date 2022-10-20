.. _beamlist:

The BEAMLIST file
=================

The BEAMLIST file is always created during initialization using the beamgen3 script.

The file contains a list of all beams used for the calculation, which is defined by the energy range (see :ref:`THEO_ENERGIES<REFENERGIES>`). Any beam not listed in BEAMLIST cannot be used in calculations, so the beams in the :ref:`IVBEAMS file<IVBEAMS>`  should always be a subset of the beams in BEAMLIST.
