.. include:: /substitutions.rst

.. _download_imaging_source_drivers:

===============================
Download Imaging Source Drivers
===============================

To use the ViPErLEED Imaging Source camera implementation, additional driver
components must be installed. These are not shipped as part of ViPErLEED due to
vendor licensing and redistribution restrictions. First, download and install
the
`GigE driver <https://www.theimagingsource.com/en-us/product/software/driver/>`__
provided by *The Imaging Source*.

In addition, ViPErLEED requires a set of Imaging Source specific ``.dll``
files. These can be obtained from our
`repository <https://github.com/viperleed/imagingsource-core>`__.

Download the repository, extract the DLL files, and place them inside the
**Drivers** directory that you configured in the
:ref:`System settings <setup_measurements>`.

Once these drivers are in place, Imaging Source cameras can be used within the
ViPErLEED measurement GUI.
