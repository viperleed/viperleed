.. include:: /substitutions.rst

.. _api_index:

====================
ViPErLEED Python API
====================

A core part of the :ref:`ViPErLEED project<index>` is the ``viperleed``
Python package. (See :ref:`here<installation>` for installation instructions.)
While ``viperleed`` can be used as a standalone tool to, e.g., perform |LEED-IV|
calculations via |calc|, it also contains a number of useful library functions
related  to surface science and electron diffraction. Through the :term:`API`,
users can access all of the functionality of ``viperleed`` from within their
own Python programs.

To use the API, you can simply import the ``viperleed`` package and call
the desired functions. Here is a simple example of how to use the API to 
find the :ref:`plane group<planegroups>` of a slab from a :ref:`POSCAR` file:

.. code-block:: python

    from viperleed.calc import symmetry
    from viperleed.calc.classes.rparams import Rparams
    from viperleed.calc.files import poscar

    # Load a slab from a POSCAR file
    example_slab = poscar.read('POSCAR')

    # Create an empty parameters object and update the slab
    run_params = Rparams()
    example_slab.full_update(run_params)

    # Detect the plane group of the slab and print it
    planegroup = symmetry.findSymmetry(example_slab, run_params)
    print(planegroup)


.. toctree::
    :maxdepth: 1

    ASE interface<api/aseapi>
    api/classes
    api/files
    api/libraries
    api/sections
    api/errors
