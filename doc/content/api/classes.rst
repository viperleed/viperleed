.. include:: /substitutions.rst

.. _api_calc_classes:

Classes
=======

The |calc| package uses a number of classes to represent the input
data, parameters and the results of |LEED-IV| calculations.

.. currentmodule:: viperleed.calc.classes

.. autosummary::
    :toctree: classes
    :recursive:

    atom.Atom
    atom_containers.AtomContainer
    atom_containers.AtomList
    beam.Beam
    layer.Layer
    r_error.R_Error
    rparams.DomainParameters
    rparams.Rparams
    searchpar.SearchPar
    sitetype.Sitetype
    sitetype.Atom_type
    slab.base_slab.BaseSlab
    slab.bulk_slab.BulkSlab
    slab.surface_slab.SurfaceSlab
    sym_entity.SymPlane
