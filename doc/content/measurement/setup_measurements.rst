.. include:: /substitutions.rst

.. _setup_measurements:

=====
Setup
=====

After ViPErLEED has been installed successfully, the graphical user interface
(GUI) can be started by executing ``viperleed gui`` in a command prompt or
terminal. This opens the main GUI window, which provides two options:
**Simulate LEED Pattern** and **Measure LEED-IV**. To proceed to the data
acquisition interface, select **Measure LEED-IV**.


.. _fig_gui_selection:
.. figure:: /_static/gui/gui_selection.jpg
    :width: 30%
    :align: center

    Press *Measure LEED-IV* to enter the GUI for data acquisition.

System settings
===============

Before measurements can be carried out, ViPErLEED must know where to store and
retrieve essential configuration files, measurement data, and (if applicable)
hardware-related drivers and components. If you are entering the measurement
GUI for the first time, a warning will appear indicating that mandatory system
settings need to be defined.

.. _fig_error_sys_settings:
.. figure:: /_static/gui/error_sys_settings.jpg
    :width: 60%
    :align: center

    A warning that mandatory system settings need to be selected.

Select *OK* to continue. This opens the **System settings** dialogue, where
paths to relevant folders or directories can be specified. The hints displayed
next to each label provide short explanations of how each directory is used
within ViPErLEED.

.. _fig_sys_settings:
.. figure:: /_static/gui/sys_settings.jpg
    :width: 55%
    :align: center

    The system settings dialogue.

The **System Settings** section allows the definition of the following paths:

 * **Configuration** – Location in which device and GUI settings are stored and from which they are loaded.
 * **Measurements** – Directory where measurement data are saved.
 * **Arduino CLI** – Path to the installation of the `Arduino CLI <https://docs.arduino.cc/arduino-cli/>`__ used for firmware uploading.
 * **Drivers** – Directory containing additional device drivers.
 * **Firmware** – Directory containing the ViPErLEED hardware controller firmware files.

Only the **Configuration** and **Measurements** paths are mandatory. However,
if you plan to use the ViPErLEED hardware controller in your measurements,
you should also define the **Arduino CLI** and **Firmware** paths at this
stage. If your setup uses our Imaging Source camera implementation, it is
recommended to specify the **Drivers** directory at this point as well.

Once all relevant paths have been defined, the main measurement GUI will
appear.

.. _fig_measurement_gui:
.. figure:: /_static/gui/measurement_gui.jpg
    :width: 40%
    :align: center

    The data acquisition interface.

The toolbar at the top provides additional functionality not directly related
to data acquisition:

 * **File** – Load existing measurements for review.
 * **Devices** – Edit camera and controller settings. Selecting a camera opens a live preview window.
 * **Tools** – Hardware-related utilities.
 * **View** – Reopen the data plot if it has been closed.
 * **Settings** – Reopen the system settings dialogue.
 * **About** – Display software and version information.

Once the system settings have been configured, ViPErLEED is ready for firmware
upload, bad pixel calibration and data acquisiotion.

Uploading Firmware
==================

Before a ViPErLEED hardware controller can be used for data acquisition, the
firmware must be uploaded to the controller. ViPErLEED uses the Arduino Micro
as its hardware controller platform, and new controllers must be flashed once
before their first use. With the controller connected to the PC, open the
**Tools** menu and select *Upload/upgrade Firmware*. Make sure that your
**Firmware** directory contains a copy of the firmware for the controller.

.. _fig_measurement_gui_select_firmware:
.. figure:: /_static/gui/measurement_gui_select_firmware.jpg
    :width: 40%
    :align: center

    Select *Upload/upgrade Firmware*.

If no `Arduino CLI <https://docs.arduino.cc/arduino-cli/>`__ installation is
available in either the system paths or the **Arduino CLI** directory,
ViPErLEED will request permission to perform a local installation of the
Arduino CLI.

.. _fig_install_arduino_cli:
.. figure:: /_static/gui/install_arduino_cli.jpg
    :width: 68%
    :align: center

    Press *Agree and install Arduino CLI*.

If your operating system does not support OpenSSL 1, you will be prompted to
install the Python ``requests`` package. If your operating system supports
OpenSSL 1, or after you have installed the ``requests`` package, the Arduino
CLI will be installed automatically and you can start setting up your
ViPErLEED hardware controller.

.. _fig_firmware_gui:
.. figure:: /_static/gui/firmware_gui.jpg
    :width: 58%
    :align: center

    The firmware upload tool.

Press *Refresh* to detect connected controllers. New Arduino Micro boards will
be shown with the COM port they are detected on. If a controller already
contains ViPErLEED firmware, this will be indicated by the name and the
installed version will be displayed. You can select the firmware you wish to
upload in the dropdown window next to the *Upload firmware* button. ViPErLEED
will automatically search your configured **Firmware** directory for matching
firmware files, and will automatically select the newest available version.
After you have selected your desired firmware, press *Upload firmware* to start
the upload. Once the progress bar reaches 100%, the controller is ready for
use.

Bad pixel calibration
=====================

This calibration step is optional, but strongly recommended before performing
|LEED-IV| measurements. Camera sensors typically contain a small number of
defective or unstable pixels. To improve image quality during data acquisition,
these pixels can be detected and replaced during image processing. To perform a
bad pixel calibration, open the *Find bad pixels* tool under the **Tools**
menu.

.. _fig_measurement_gui_select_pixels:
.. figure:: /_static/gui/measurement_gui_select_pixels.jpg
    :width: 40%
    :align: center

    The bad pixel calibration tool.

After connecting the camera and allowing it to fully initialise, it will be
detected and can be selected in the dropdown menu at the top of the tool.

.. _fig_find_bad_pixels:
.. figure:: /_static/gui/find_bad_pixels.jpg
    :width: 33%
    :align: center

    The bad pixel finder.

When selecting a camera for the first time, you will be prompted to choose a
directory in which the resulting bad pixel files will be stored.

.. _fig_find_bad_pixels_with_camera:
.. figure:: /_static/gui/find_bad_pixels_with_camera.jpg
    :width: 45%
    :align: center

    If the bad pixel calibration has been performed previously, the stored
    results will be shown.

Press *Find* to start the acquisition process. The tool will guide you through
the required steps interactively.

.. _fig_find_bad_pixels_started:
.. figure:: /_static/gui/find_bad_pixels_started.jpg
    :width: 50%
    :align: center

    Instructions will appear as each step is performed.

Once the calibration is completed, the stored calibration will automatically be
used for subsequent measurements.
