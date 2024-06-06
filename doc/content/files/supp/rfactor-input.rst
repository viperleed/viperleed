.. include:: /substitutions.rst

.. _rfactor-input:

|R-factor| input files
======================

-   rfactor-WEXPEL: Defines how to calculate the |R factor|: The energy range
    (from either :ref:`THEO_ENERGIES<theo_energies>` or
    :ref:`EXPBEAMS.csv<EXPBEAMS>`), :ref:`V0_IMAG<v0_imag>`,
    :ref:`IV_SHIFT_RANGE<IVSHIFTRANGE>`, :ref:`R_FACTOR_TYPE<RFACTORTYPE>`,
    :ref:`R_FACTOR_SMOOTH<RFACTORSMOOTH>`, and the relationship of experimental
    and theoretical beams (including how to average them). Also contains the
    entire :ref:`AUXEXPBEAMS<AUXEXPBEAMS>`  file.
-   rfactor-PARAM: Defines static fortran array sizes. Generated automatically
    from other input, compiled with the rfactor fortran files at runtime.
