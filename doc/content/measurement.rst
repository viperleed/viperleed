.. include:: /substitutions.rst

.. _measurement:

==========================
ViPErLEED data acquisition
==========================

|measure| is the Python package that contains all functionality revolving
around data acquisition. It provides a graphical user interface through which
various calibration tasks to improve data quality and |LEED IV| measurements
can be performed.

It is available from the Python Package Index (PyPI) and can easily be
installed using ``pip``:

.. code-block:: bash

    pip install viperleed

The ViPErLEED Python code can also be obtained from the
`viperleed <https://github.com/viperleed/viperleed>`__ GitHub repository.

.. toctree::
    :maxdepth: 1

    measurement/setup_measurements
    measurement/download_firmware
    measurement/download_imaging_source_drivers
    measurement/energy_calibration
    measurement/time_resolved
    measurement/best_practice
    measurement/quantities
