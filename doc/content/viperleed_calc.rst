.. include:: /substitutions.rst

.. _viperleed_calc:

==============
viperleed.calc
==============

|calc| is the Python package that makes up the core of the computational part
of the ViPErLEED project. It provides a set of tools for calculating and
analyzing |LEED-IV| curves for arbitrary input structures.

It is available from the Python Package Index (PyPI) and can easily be
installed using ``pip``:

.. code-block:: bash

    pip install viperleed

The ViPErLEED Python code can also be obtained from the
`viperleed <https://github.com/viperleed/viperleed>`__ GitHub repository.

The |calc| package is written in Python, but many core parts require
the :term:`TensErLEED` source code to be available. It can be downloaded
from the
`viperleed-tensorleed <https://github.com/viperleed/viperleed-tensorleed>`__
GitHub repository. See the :ref:`installation page<installation>` for further
information.

Stand-alone calculations with |calc| can be executed using the
command-line interface (see :ref:`here<how_to_run>`).

An interface to the Atomic Simulation Environment (:term:`ASE`) is
:ref:`also available<aseapi>`.


.. todo:
    Mention you can import the package and use it in your own code.
