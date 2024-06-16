.. include:: /substitutions.rst

.. _rfactor-input:

|R-factor| input files
======================

-   rfactor-WEXPEL: Defines how to calculate the |R factor|: The energy range
    (from either :ref:`THEO_ENERGIES` or :ref:`EXPBEAMS.csv<EXPBEAMS>`), 
    :ref:`V0_IMAG`, :ref:`IVSHIFTRANGE`, :ref:`RFACTORTYPE`,
    :ref:`RFACTORSMOOTH`, and the relationship of experimental
    and theoretical beams (including how to average them). Also 
    contains the entire :ref:`AUXEXPBEAMS`  file.
-   rfactor-PARAM: Defines static fortran array sizes. Generated automatically
    from other input, compiled with the rfactor fortran files at runtime.
