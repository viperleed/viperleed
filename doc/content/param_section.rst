.. _paramsection:

List of parameters by section
=============================

This page lists parameters for ViPErLEED by which section of ViPErLEED/TensErLEED they are (most) relevant for. Parameters may occur multiple times, and not all parameters are listed.
See the :ref:`main PARAMETERS page<PARAMETERS>` for a complete alphabetical list and other groupings.

.. note::
    While all parameters have a default value, parameters marked with a 
    "**→**" usually require user input for simple but non-trivial 
    systems.

Initialization
--------------

.. table::
  :width: 100%
  :widths: 25 75

  +--------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | Parameter                                  | Function                                                                                                     |
  +============================================+==============================================================================================================+
  | :ref:`BEAM_INCIDENCE<BEAMINCIDENCE>`       | Incidence angle and direction of the electron beam in experiment                                             |
  +--------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`BULK_LIKE_BELOW<BULK_LIKE_BELOW>`    | Onset of unrelaxed slab, for automatic detection of minimal bulk and bulk repeat                             |
  +--------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`BULK_REPEAT<BULK_REPEAT>`            | Thickness of the bulk repeat unit, or a bulk repeat vector                                                   |
  +--------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`ELEMENT_MIX<ELSPLIT>`                | Declare that sites in the :ref:`POSCAR file<POSCAR>`  can be occupied by different chemical elements         |
  +--------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`LAYER_CUTS<layer_cuts>`              | How to separate the :ref:`POSCAR file<POSCAR>`  into layers                                                  |
  +--------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | **→** :ref:`N_BULK_LAYERS<n_bulk_layers>`  | Define how many layers in the :ref:`POSCAR file<POSCAR>`  represent the bulk                                 |
  +--------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | **→** :ref:`SITE_DEF<SITEDEF>`             | Define which sites in the :ref:`POSCAR file<POSCAR>`  are special, i.e. have different vibrational amplitude |
  +--------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | **→** :ref:`SUPERLATTICE<SUPERLATTICE>`    | The relationship between the surface and bulk unit cells                                                     |
  +--------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`SYMMETRIZE_INPUT<SYMMETRY_NOMOVE>`   | Whether to move atoms in the :ref:`POSCAR file<POSCAR>`  to perfectly match the symmetry                     |
  +--------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`SYMMETRY_EPS<sym_eps>`               | Error tolerance during symmetry search                                                                       |
  +--------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`SYMMETRY_FIX<ISYM>`                  | Manually set a symmetry, or turn symmetry off                                                                |
  +--------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`SYMMETRY_FIND_ORI<SYMMETRY_FIND_ORI>`| Whether the symmetry search should look for the highest-symmetry point.                                      |
  +--------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`THEO_ENERGIES<theo_energies>`        | What energy range to calculate                                                                               |
  +--------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`T_DEBYE<T_DEBYE>`                    | Debye temperature of the system (only for automatically generating :ref:`VIBROCC<VIBOCCIN>`)                 |
  +--------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`T_EXPERIMENT<T_EXPERIMENT>`          | Measurement temperature in experiment (only for automatically generating :ref:`VIBROCC<VIBOCCIN>`)           |
  +--------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`VIBR_AMP_SCALE<VIBR_AMP_SCALE>`      | Scaling factor, only for automatically generating :ref:`VIBROCC<VIBOCCIN>`                                   |
  +--------------------------------------------+--------------------------------------------------------------------------------------------------------------+

.. note::
  Parameters setting the symmetry strongly affect all sections, but are 
  not listed again below.

Reference calculation
---------------------

.. table::
  :width: 100%
  :widths: 25 75

  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | Parameter                                   | Function                                                                                                     |
  +=============================================+==============================================================================================================+
  | :ref:`ATTENUATION_EPS<attenuation_eps>`     | Cutoff for beam propagation                                                                                  |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`BEAM_INCIDENCE<BEAMINCIDENCE>`        | Incidence angle and direction of the electron beam in experiment                                             |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`BULKDOUBLING_EPS<BULKDOUBLEEPS>`      | Convergence criterion for bulk thickness in the TensErLEED calculation                                       |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`BULKDOUBLING_MAX<BULKDOUBLEITER>`     | Maximum bulk thickness in TensErLEED calculation                                                             |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`BULK_REPEAT<BULK_REPEAT>`             | Thickness of the bulk repeat unit, or a bulk repeat vector                                                   |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`ELEMENT_MIX<ELSPLIT>`                 | Declare that sites in the :ref:`POSCAR file<POSCAR>`  can be occupied by different chemical elements         |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`FILAMENT_WF<FILWF>`                   | The LEED filament work function                                                                              |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`KEEP_REFCALC_DIRS<keep_refcalc_dirs>` | Toggle to keep the reference calculating execution directories                                               |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`LAYER_CUTS<layer_cuts>`               | How to separate the :ref:`POSCAR file<POSCAR>`  into layers                                                  |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`LMAX<LMAX>`                           | Maximum angular momentum number; usually determined via :ref:`PHASESHIFT_EPS<PHASESHIFTMIN>`                 |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | **→** :ref:`N_BULK_LAYERS<n_bulk_layers>`   | Define how many layers in the :ref:`POSCAR file<POSCAR>`  represent the bulk                                 |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | **→** :ref:`N_CORES<NCORES>`                | The number of CPUs to use                                                                                    |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`PHASESHIFT_EPS<PHASESHIFTMIN>`        | Cutoff in phaseshifts magnitudes to determine :ref:`LMAX<LMAX>`                                              |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`SCREEN_APERTURE<SCREEN_APERTURE>`     | The aperture of the acceptance cone of the LEED screen                                                       |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | **→** :ref:`SITE_DEF<SITEDEF>`              | Define which sites in the :ref:`POSCAR file<POSCAR>`  are special, i.e. have different vibrational amplitude |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | **→** :ref:`SUPERLATTICE<SUPERLATTICE>`     | The relationship between the surface and bulk unit cells                                                     |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`TENSOR_OUTPUT<TOUTPUT>`               | Disable Tensor output for some or all layers                                                                 |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`THEO_ENERGIES<theo_energies>`         | What energy range to calculate                                                                               |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`V0_IMAG<v0_imag>`                     | Imaginary part of the inner potential                                                                        |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`V0_REAL<MUFTIN>`                      | Real part of the inner potential                                                                             |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`V0_Z_ONSET<INPOTZ>`                   | How far from the topmost atom the inner potential begins                                                     |
  +---------------------------------------------+--------------------------------------------------------------------------------------------------------------+

R-factor calculation
--------------------

.. table::
  :width: 100%
  :widths: 25 75

  +----------------------------------------+--------------------------------------------------------------------------------------------------+
  | Parameter                              | Function                                                                                         |
  +========================================+==================================================================================================+
  | :ref:`BEAM_INCIDENCE<BEAMINCIDENCE>`   | Incidence angle and direction of the electron beam in experiment                                 |
  +----------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`IV_SHIFT_RANGE<IVSHIFTRANGE>`    | Range and step size for shifting experimental and theoretical curves during R-factor calculation |
  +----------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`PLOT_IV<PLOT_COLORS_RFACTOR>`    | Change appearance of the R-factor plot files                                                     |
  +----------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`R_FACTOR_TYPE<RFACTORTYPE>`      | Which definition of the R-factor to use                                                          |
  +----------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`R_FACTOR_SMOOTH<RFACTORSMOOTH>`  | How strongly experimental beams are smoothed                                                     |
  +----------------------------------------+--------------------------------------------------------------------------------------------------+
  | **→** :ref:`SUPERLATTICE<SUPERLATTICE>`| The relationship between the surface and bulk unit cells                                         |
  +----------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`THEO_ENERGIES<theo_energies>`    | What energy range to calculate                                                                   |
  +----------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`V0_IMAG<v0_imag>`                | Imaginary part of the inner potential                                                            |
  +----------------------------------------+--------------------------------------------------------------------------------------------------+

Delta-amplitudes calculation
----------------------------

Behaviour is mainly governed by the :ref:`DISPLACEMENTS file<DISPLACEMENTS>`. Some relevant parameters are:

.. table::
  :width: 100%
  :widths: 25 75

  +----------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | Parameter                              | Function                                                                                                     |
  +========================================+==============================================================================================================+
  | :ref:`ELEMENT_MIX<ELSPLIT>`            | Declare that sites in the :ref:`POSCAR file<POSCAR>`  can be occupied by different chemical elements         |
  +----------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`LMAX<LMAX>`                      | Maximum angular momentum number; usually determined via :ref:`PHASESHIFT_EPS<PHASESHIFTMIN>`                 |
  +----------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | **→** :ref:`N_CORES<NCORES>`           | The number of CPUs to use                                                                                    |
  +----------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`PHASESHIFT_EPS<PHASESHIFTMIN>`   | Cutoff in phaseshifts magnitudes to determine :ref:`LMAX<LMAX>`                                              |
  +----------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | **→** :ref:`SITE_DEF<SITEDEF>`         | Define which sites in the :ref:`POSCAR file<POSCAR>`  are special, i.e. have different vibrational amplitude |
  +----------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | **→** :ref:`SUPERLATTICE<SUPERLATTICE>`| The relationship between the surface and bulk unit cells                                                     |
  +----------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`TENSOR_INDEX<TENSOR_INDEX>`      | Which Tensor files to use for the delta calculation and search                                               |
  +----------------------------------------+--------------------------------------------------------------------------------------------------------------+
  | :ref:`THEO_ENERGIES<theo_energies>`    | What energy range to calculate                                                                               |
  +----------------------------------------+--------------------------------------------------------------------------------------------------------------+

Search
------

Behaviour is also governed by the :ref:`DISPLACEMENTS file<DISPLACEMENTS>`. The most relevant parameters are:

.. table::
  :width: 100%
  :widths: 25 75

  +------------------------------------------------+--------------------------------------------------------------------------------------------------+
  | Parameter                                      | Function                                                                                         |
  +================================================+==================================================================================================+
  | :ref:`BEAM_INCIDENCE<BEAMINCIDENCE>`           | Incidence angle and direction of the electron beam in experiment                                 |
  +------------------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`IV_SHIFT_RANGE<IVSHIFTRANGE>`            | Range and step size for shifting experimental and theoretical curves during R-factor calculation |
  +------------------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`LOG_SEARCH<LOG_SEARCH>`                  | Output the search log file (may be very large, mostly for debugging)                             |
  +------------------------------------------------+--------------------------------------------------------------------------------------------------+
  | **→** :ref:`N_CORES<NCORES>`                   | The number of CPUs to use                                                                        |
  +------------------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`R_FACTOR_TYPE<RFACTORTYPE>`              | Which definition of the R-factor to use                                                          |
  +------------------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`R_FACTOR_SMOOTH<RFACTORSMOOTH>`          | How strongly experimental beams are smoothed                                                     |
  +------------------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`SEARCH_BEAMS<SEARCHBEAMS>`               | Whether to use R-factor of integer, fractional, or all beams for the search                      |
  +------------------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`SEARCH_CONVERGENCE<SEARCH_CONVERGENCE>`  | Convergence criteria for the search, and convergence-dependent parameter control                 |
  +------------------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`SEARCH_CULL<SEARCH_CULL>`                | Controls regular culling of worst-performing structures, and what to replace them with           |
  +------------------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`SEARCH_MAX_GEN<SEARCHGENMAX>`            | Maximum total number of generations that the search should run for                               |
  +------------------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`SEARCH_POPULATION<SEARCHPOP>`            | Number of trial structures used in the search                                                    |
  +------------------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`SEARCH_START<SEARCHSTART>`               | How to initialize the search population                                                          |
  +------------------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`STOP<STOP>`                              | Stop execution of ViPErLEED at next opportunity                                                  |
  +------------------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`TENSOR_INDEX<TENSOR_INDEX>`              | Which Tensor files to use for the delta calculation and search                                   |
  +------------------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`THEO_ENERGIES<theo_energies>`            | What energy range to calculate                                                                   |
  +------------------------------------------------+--------------------------------------------------------------------------------------------------+
  | :ref:`V0_IMAG<v0_imag>`                        | Imaginary part of the inner potential                                                            |
  +------------------------------------------------+--------------------------------------------------------------------------------------------------+

Domain search
-------------

As the :ref:`domain search<domain_calculation>`  may involve all of the segments above, the parameters listed there are relevant. The following additional parameters affect domains specifically:

.. table::
  :width: 100%
  :widths: 25 75
  
  +----------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------+
  | Parameter                                                | Function                                                                                                                                 |
  +==========================================================+==========================================================================================================================================+
  | :ref:`DOMAIN<DOMAIN>`                                    | Define a domain for :ref:`calculations involving multiple coexisting structural domains<domain_calculation>`                             |
  +----------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------+
  | :ref:`DOMAIN_STEP<DOMAIN_STEP>`                          | Step width for structural domain coverage during search                                                                                  |
  +----------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------+
  | :ref:`SYMMETRY_CELL_TRANSFORM<SYMMETRY_CELL_TRANSFORM>`  | Relationship between a supercell and the primitive surface unit cell (only relevant for :ref:`domain calculations<domain_calculation>`)  |
  +----------------------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------+
