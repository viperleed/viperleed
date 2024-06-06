.. include:: /substitutions.rst

.. _search-input:

Search input files
==================

-   search.steu: Lists some search parameters (from
    :ref:`SEARCH_CONVERGENCE<SEARCH_CONVERGENCE>` and
    :ref:`SEARCH_MAX_GEN<SEARCHGENMAX>`), which |R factor| to use
    (\ :ref:`R_FACTOR_TYPE<RFACTORTYPE>`), the atom variations
    (as determined from :ref:`DISPLACEMENTS<DISPLACEMENTS>`),
    and the starting configuration (see :ref:`SEARCH_START<SEARCHSTART>`).
-   rf.info: Contains the information for calculating the |R factor| during the
    search: The energy range (combined from :ref:`THEO_ENERGIES<theo_energies>`
    and :ref:`EXPBEAMS.csv<EXPBEAMS>`), :ref:`V0_IMAG<v0_imag>`,
    :ref:`IV_SHIFT_RANGE<IVSHIFTRANGE>`, :ref:`R_FACTOR_SMOOTH<RFACTORSMOOTH>`,
    and the relationship of experimental and theoretical beams (including how
    to average them). Also contains the entire :ref:`AUXEXPBEAMS<AUXEXPBEAMS>`
    file.
-   restrict.f: Fortran source code file. Contains additional parameter
    constraints, based on the :ref:`DISPLACEMENTS file<DISPLACEMENTS>`.
    Compiled with the rest of the search fortran files at runtime.
-   search-PARAM: Defines static fortran array sizes. Mostly generated
    automatically from all other input, but also defines the search
    population size :ref:`SEARCH_POPULATION<SEARCHPOP>`. Compiled with
    the search fortran files at runtime.
