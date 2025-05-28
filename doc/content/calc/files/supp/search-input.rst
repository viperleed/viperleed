.. include:: /substitutions.rst

.. _search-input:

Search input files
==================

-   search.steu: Lists some search parameters (from :ref:`SEARCH_CONVERGENCE`
    and :ref:`SEARCHGENMAX`), which |R factor| to use (\ :ref:`RFACTORTYPE`),
    the atom variations (as determined from :ref:`DISPLACEMENTS`), and the
    starting configuration (see :ref:`SEARCHSTART`).
-   rf.info: Contains the information for calculating the |R factor| during the
    search: The energy range (combined from :ref:`THEO_ENERGIES`
    and :ref:`EXPBEAMS.csv<EXPBEAMS>`), :ref:`V0_IMAG`, :ref:`IVSHIFTRANGE`,
    :ref:`RFACTORSMOOTH`, and the relationship of experimental and theoretical
    beams (including how to average them). Also contains the entire
    :ref:`AUXEXPBEAMS` file.
-   restrict.f: Fortran source code file. Contains additional parameter
    constraints, based on the :ref:`DISPLACEMENTS` file.
    Compiled with the rest of the search fortran files at runtime.
-   search-PARAM: Defines static fortran array sizes. Mostly generated
    automatically from all other input, but also defines the search population 
    size :ref:`SEARCHPOP`. Compiled with the search fortran files at runtime.
