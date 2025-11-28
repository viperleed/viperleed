.. include:: /substitutions.rst

.. _setup_measurements:

=====
Setup
=====

After ViPErLEED has been installed successfully, the graphical user interface
(GUI) can be started by executing ``viperleed gui`` in a command prompt or
terminal. This opens the main GUI window, which provides two options:
**Simulate LEED Pattern** and **Measure LEED-IV**. To proceed to the data
acquisition interface, select **Measure LEED-IV**. If your GUI is rendered
differently, set the DPI settings of python.exe to override by application.

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

Select *OK* to continue. This opens the **System settings** dialog, where
paths to relevant folders or directories can be specified. The hints displayed
next to each label provide short explanations of how each directory is used
within ViPErLEED.

.. _fig_sys_settings:
.. figure:: /_static/gui/sys_settings.jpg
    :width: 55%
    :align: center

    The system settings dialog.

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

 * **File** – Load existing measurements for review. Loaded measurements are displayed in the data plot.
 * **Devices** – Edit camera and controller settings. Selecting a camera opens a live preview window. Settings for new devices can be created by selecting the corresponding device from this menu.
 * **Tools** – Hardware-related utilities, such as the firmware upload tool for the ViPErLEED hardware controller and a bad-pixel calibration tool.
 * **View** – Reopen the data plot if it has been closed.
 * **Settings** – Reopen the system settings dialog.
 * **About** – Display software information.

Once the system settings have been configured, ViPErLEED is ready for firmware
upload, bad pixel calibration, and data acquisition.

Uploading firmware
==================

Before a ViPErLEED hardware controller can be used for data acquisition, the
firmware must be uploaded to the controller. ViPErLEED uses the Arduino Micro
as its hardware controller platform, and new controllers must be flashed once
before their first use. With the controller connected to the PC, open the
**Tools** menu and select *Upload/upgrade Firmware*. Make sure that your
**Firmware** directory contains a copy of the controller firmware.

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
be listed with the COM port on which they were detected. If a controller
already contains ViPErLEED firmware, this will be indicated by its name and the
installed version will be displayed. You can select the firmware you wish to
upload in the dropdown window next to the *Upload firmware* button. ViPErLEED
will automatically search your configured **Firmware** directory for matching
firmware files, and will automatically select the newest available version.
After you have selected your desired firmware, press *Upload firmware* to start
the upload. Once the progress bar reaches 100%, the firmware has been
installed on the controller.

.. _fig_upload_finished:
.. figure:: /_static/gui/upload_finished.jpg
    :width: 58%
    :align: center

    The upload progress has to reach 100%. If it fails and reports an error for
    any reason, the upload has to be repeated.

ViPErLEED hardware controller setup
===================================

To operate the controller during a measurement, it must be identifiable. To
assign your controller a serial number and to create settings for it, open the
**Devices** menu and select the controller on which the firmware was just
installed.

.. _fig_select_controller:
.. figure:: /_static/gui/select_controller.jpg
    :width: 53%
    :align: center

    If the controller has not been used before, it will not yet have a serial
    number.

You should receive a notification that no settings were found for this
controller. Press *Create a new settings file* to generate new configuration
files for the controller in your **Configuration** directory.

.. _fig_create_new_settings_no_serial_nr:
.. figure:: /_static/gui/create_new_settings_no_serial_nr.jpg
    :width: 65%
    :align: center

    Create new settings for your controller.

After creating new settings, assign a serial number to the controller. Either
generate a random serial number by pressing *Generate randomly* or type your
desired serial number into the line edit. Press *Set* to write it to the
controller. This serial number will act as the identification marker for the
controller, but it can be changed later if needed.

.. _fig_new_settings_no_serial_nr:
.. figure:: /_static/gui/new_settings_no_serial_nr.jpg
    :width: 45%
    :align: center

    Set a serial number on your controller.

After setting the serial number, delete the old, temporary settings file and
close the settings dialog. Ensure the controller is powered, then reopen the
controller settings dialog. The **Devices** menu should now display the
controller with the serial number you assigned it.

.. _fig_delete_old_settings:
.. figure:: /_static/gui/delete_old_settings.jpg
    :width: 65%
    :align: center

    The old settings file can be deleted; it corresponds to a controller
    without a serial number.

.. _fig_new_settings_with_serial_nr:
.. figure:: /_static/gui/new_settings_with_serial_nr.jpg
    :width: 48%
    :align: center

    Once the controller is powered, the hardware configuration will be
    displayed.

Click on the *Temperature (°C) - ?? thermocouple* dropdown menu and click the
option with the same content in the dropdown, as seen in
:numref:`_fig_dropdown_thermocouple`.

.. _fig_dropdown_thermocouple:
.. figure:: /_static/gui/dropdown_thermocouple.jpg
    :width: 30%
    :align: center

    Select the thermocouple option.

In the window that opens, select the type of thermocouple installed on your
controller.

.. _fig_select_thermocouple:
.. figure:: /_static/gui/select_thermocouple.jpg
    :width: 45%
    :align: center

    Select the installed thermocouple type.

When finished, press *OK* to save your settings.

Before acquiring data with your ViPErLEED controller, verify that it is
configured for the correct |I0| input range for your LEED unit. If the wrong
mode is selected, you must adjust the jumpers inside the ViPErLEED hardware
controller box.

Camera setup
============

To set up your camera, first ensure that it is powered and connected to your
PC. Give the camera enough time to boot. Go to the **Devices** menu and select
your camera to open a live view.

.. _fig_select_camera:
.. figure:: /_static/gui/select_camera.jpg
    :width: 48%
    :align: center

    Select the camera.

You should receive a notification that no settings were found for this camera.
Press *Create a new settings file* to generate new configuration files for the
camera in your **Configuration** directory.

.. _fig_create_new_settings_camera:
.. figure:: /_static/gui/create_new_settings_camera.jpg
    :width: 65%
    :align: center

    Create new settings for your camera.

Once the settings have been created, a live view of the camera will open.
Right-click inside the live view.

.. _fig_camera_view:
.. figure:: /_static/gui/camera_view.jpg
    :width: 100%
    :align: center

    Right-click in the live view to open its menu.

Select *Properties* from the menu to open the camera settings. Set the binning
value such that your images have between 400×400 and 500×500 pixels after
binning. If your sensor is rectangular, set the binning so that the shorter
side is between 400 and 500 pixels. Later, we will set the ROI so that the
resulting images are square.

.. _fig_new_camera_settings:
.. figure:: /_static/gui/new_camera_settings.jpg
    :width: 55%
    :align: center

    Set the binning value to reduce image size and improve the signal-to-noise
    ratio.

Bad pixel calibration
=====================

This calibration step is optional but strongly recommended before performing
|LEED-IV| measurements. Camera sensors typically contain a small number of
defective or unstable pixels. To improve image quality during data acquisition,
these pixels can be detected and replaced during image processing. To perform a
bad pixel calibration, open the bad pixel calibration tool under the **Tools**
menu.

.. _fig_measurement_gui_select_pixels:
.. figure:: /_static/gui/measurement_gui_select_pixels.jpg
    :width: 40%
    :align: center

    Select *Find bad pixels*.

The camera can be selected in the dropdown menu at the top of the tool.

.. _fig_find_bad_pixels:
.. figure:: /_static/gui/find_bad_pixels.jpg
    :width: 33%
    :align: center

    The dropdown menu is populated automatically. If you just connected your
    camera to the PC, give it some time to initialise and then check the menu
    again.

When selecting a camera for the first time, you will need to choose a directory
in which the resulting bad pixel files will be stored.

.. _fig_find_bad_pixels_with_camera:
.. figure:: /_static/gui/find_bad_pixels_with_camera.jpg
    :width: 45%
    :align: center

    If the bad pixel calibration has been performed previously, the stored
    results will be shown.

Press *Find* to start the acquisition process. The tool will guide you through
the required steps interactively. First you will have to acquire dark images.
For this you should have the cap of the camera lens ready, to ensure that you
can cover it entirely.

.. _fig_bad_pixels_instruction_1:
.. figure:: /_static/gui/bad_pixels_instruction_1.jpg
    :width: 50%
    :align: center

    Instructions will appear as each step is performed.

There are two dark image steps in which different pixel errors are detected.
After this you will have to take off the cap of the camera lens again and you
need to provide uniform illumination. Usually just turning the lights on and
putting a sheet of white paper a few centimeters infront of the lense is good
enough. Make sure that there are no shadows cast onto the sheet.

.. _fig_bad_pixels_instruction_2:
.. figure:: /_static/gui/bad_pixels_instruction_2.jpg
    :width: 50%
    :align: center

    Take off the cap and put a white sheet of paper infront of the camera.

Once the calibration is completed, the calibration will be stored in your
chosen directory and is automatically used in subsequent measurements.

.. _fig_bad_pixels_finished
.. figure:: /_static/gui/bad_pixels_finished.jpg
    :width: 50%
    :align: center

    After completion you can close the bad-pixel finder dialog.



Camera placement
================

After pixel calibration, the camera can be moved into its final position, as it
will no longer be necessary to adjust it. The camera aperture should always be
fully opened. To proceed, open a live view of your camera. For optimal
settings, position the camera so that the LEED screen fills out the entire
sensor. (Align the edges of the screen with the edges of the sensor.) After
adjusting the camera position, right-click the camera view and select
“Allow setting of ROI”. Set the region of interest such that the edges of the
ROI align with the edges of the screen.
