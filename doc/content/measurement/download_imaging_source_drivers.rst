.. include:: /substitutions.rst

.. _download_imaging_source_drivers:

===============================
Download Imaging Source drivers
===============================

To use the ViPErLEED Imaging Source camera implementation, additional driver
components must be installed. These are not shipped as part of ViPErLEED due to
vendor licensing and redistribution restrictions. First, download and install
the
`GigE driver <https://www.theimagingsource.com/en-us/product/software/driver/>`__
provided by *The Imaging Source*. In addition, ViPErLEED requires a set of
Imaging Source specific ``.dll`` files. These can be obtained from our
`repository <https://github.com/viperleed/imagingsource-core>`__. Download the
repository, extract the DLL files, and place them inside the **Drivers**
directory that you configured in the
:ref:`System settings <setup_measurements>`.

Once these drivers are in place, Imaging Source cameras can be used within the
ViPErLEED measurement GUI. If you are using a GigE camera and you have
difficulty detecting the device, try setting the IP configuration mode to DHCP
in the GigECam IP Config tool. If you are connecting the cameras via ethernet,
ensure that the network interface is configured to use jumbo frames and set the
MTU value to 9014 bytes.
