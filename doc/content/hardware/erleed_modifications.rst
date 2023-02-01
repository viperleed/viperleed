.. _erleed_modification:

#########################
ErLEED Modification Guide
#########################

The :term:`ErLEED` LEED electronics produced by :term:`SPECS` are commonly used in many surface science laboratories.
The ViPErLEED electronics are designed (and tested) to work with the ErLEED electronics.
However, while the ErLEED electronics enable external control over beam energy, they do not allow direct read-out of the applied voltage.

When available, ViPErLEED can directly measure the beam voltage and use this value to automatically calibrate the electron energy.
This is not strictly required, but **highly recommend** to prevent distortions of the :math:`I(V)` curves.

It is possible to make the beam voltage accessible for read out in the ErLEED electronics by making a small modification to the control unit.

.. warning::
    The modifications described below **will** void any warranty unless you explicitly get a permission by the supplier.
    The ViPErLEED developers take no responsibility for any malfunction of that may occur as a result of modifications to the LEED electronics.


**TODO Michele, Michael, Alex: Details on which versions numbers are supported; why we perform the modifications, warnings etc.**

Required Components
===================

Before you start with the modifications, make sure you have all required components and tools.


For the beam HV pin:
    - 2 330 :math:`\Omega` resistors,
    - **TODO: name of pin component**
    - a short shrink tube (~1.5 cm),

For the beam HV port:
    - **TODO: cable, components**


Additionally, you will need the following tools to perform the modifications:
    - a suitable set of Phillips and flat head screw drivers,
    - a soldering iron and solder,
    - tweezers,
    - cable clippers,
    - **TODO: name of hole punch machine**
    - a cable tie cutter,
    - a heat gun for shrink tubes.

.. _fig_resistors_1:
.. figure:: /_static/hardware/ErLEED_modification/resistors/resistors_1.png
    :width: 20%
    :align: center

    Components for the beam HV pin.

**TODO Michele, Michael: names/numbers of HV connector components & hole-punch machine**

Opening up the electronics
==========================

.. important::
    Before starting the modifications described below, turn off and completely disconnect all plugs from the control unit.
    **Disconnect the power cable** and wait 10 minutes before proceeding with opening the unit to allow all capacitors to fully discharge.


To start, we need to open up the electronics.
First, fully disconnect all plugs from the unit and place it on a suitable anti-static electronics workbench.
Start to open up the electronics unit by removing all screws that hold the cover plate in place, then remove the cover plate.
:numref:`fig_cover_plate_removed` shows the inside of the ErLEED control unit with the cover plate removed.


.. _fig_cover_plate_removed:
.. figure:: /_static/hardware/ErLEED_modification/electronics_overview.svg
    :width: 75%
    :align: center

    ErLEED control unit with cover plate removed.

Removing the back plate
=======================

To make the required modifications, it is also necessary to partially take off the back plate of the control unit.
There is no need to completely disconnect the back plate from the rest of the electronics, but fashioning the new port, as described below will likely require tiling the plate horizontally.

For best accessibility, you most likely want to remove the power plug and screen connection from the back plate (compare :numref:`fig_new_port_location`).
Additionally, you may need to remove a number of cable ties and unplug various connections to the motherboard (e.g. connections for "ANODE", "L2", "L1/3", "FILAMENT", etc.).
We highly recommend labeling each cable before removal.

Ultimately, you should be able to remove and tilt back the back plate of the unit as shown in :numref:`fig_taking_off_backplate`.


.. _fig_taking_off_backplate:
.. figure:: /_static/hardware/ErLEED_modification/taking_off_backplate.svg
    :width: 75%
    :align: center

    Removing the back plate of the control unit.


Beam HV pin
===========

In the following steps, you will build and place a simple pin connector which allows reading out the beam potential.
A suitable pin can easily be improvised from two 330 :math:`\Omega` resistors, a (?) and a short shrink tube as shown in :numref:`fig_resistors_1`.


First, bend the wires of the resistors around the (?) as shown in :numref:`fig_resistors_1`.
Then, place the (?) in the loop and solder it in place and repeat the procedure with the second resistor (see :numref:`fig_resistors_3` and :numref:`fig_resistors_4`).
Make sure the soldered spot is stable and the connection across both resistors and to the pin is intact.
Finally, cut off protruding wires (see :numref:`fig_resistors_4`) from the resistors and place a shrink tube over the resistors and soldered spot (see :numref:`fig_resistors_5`).


.. list-table::
    :align: center
    :width: 100%

    * - .. _fig_resistors_2:
  
        .. figure:: /_static/hardware/ErLEED_modification/resistors/resistors_2.png

            Bending resistor wire around the pin.

      - .. _fig_resistors_3:

        .. figure:: /_static/hardware/ErLEED_modification/resistors/resistors_3.png

            Two resistors soldered to the pin.

      - .. _fig_resistors_4:

        .. figure:: /_static/hardware/ErLEED_modification/resistors/resistors_4.png

            Protruding wires removed.

      - .. _fig_resistors_5:

        .. figure:: /_static/hardware/ErLEED_modification/resistors/resistors_5.png

            Finished connector with shrink tube.

To read out the beam voltage during operation, we need access to the new pin from the outside.
For this, there needs to be a connection from the pin on the motherboard to the backside of the chassis.
To make the connection cable, take the **TODO cable specifications**, strip off ~0.5 cm of the insulation on both side, and solder the **TODO name of female pin part** to the cable.
You should be able to firmly connect the cable to the pin as shown in :numref:`fig_resistors_cable_attached`.

.. _fig_resistors_cable_attached:
.. figure:: /_static/hardware/ErLEED_modification/resistors_cable_attached.jpeg
    :width: 25%
    :align: center

    New beam HV pin with connection cable attached.


New port
========


Next, we will need to fashion a new port on the control unit cassis.
Fortunately, there is ample space on the backplate, next to the existing connectors.
We recommend placing the "Beam HV" port below the ":math:`I0_{\text{MON}}`" port, as shown in :numref:`fig_new_port_location`.

Carefully punch a hole in the backplate of the chassis in the desired location.
This may require removing the power plug and screen connector from the backplate.
Then, place the new high-voltage coaxial connector, solder the connection wire to it and place a shrink tube over he solder spot.
Secure the coaxial connector in place by tightening the nut that came with it on the inside of the backplate.

Finally, re-attach the power and screen plugs if you had to remove them.
At this point, we highly recommend labeling the newly fashioned port appropriately (e.g. "Beam HV").


**TODO Michele: How is this port & machine called?**


.. _fig_new_port_location:
.. figure:: /_static/hardware/ErLEED_modification/new_port_location.svg
    :width: 75%
    :align: center

    Location of the new beam HV port.



Placing the beam HV pin
=======================

Next, you need to place the custom pin on the motherboard.
The pin is intended to read out the potential applied to the electron beam in the LEED setup.
To do this, we can measure the voltage at the filament where the electrons originate.
By placing a voltage divider with two equal resistors parallel to the filament, we can get a reference potential.
Since the filament is essentially a short circuit, a negligible current will pass through the parallel resistors.

To directly access the filament, we can conveniently place the new pin right next to the high-voltage diodes next to the filament port.
The exact location is shown in :numref:`fig_pin_location_zoomed_out` and :numref:`fig_pin_location_zoomed_in`.

Using a long needle, form two hooks with "feet" of the beam HV pin made in the last step.
The hooks should tightly loop around the connections of the diodes to the circuit board.
Once securely placed, solder both wires to the side of the diodes, while being careful not to heat the diodes directly.


.. list-table::
    :align: center
    :width: 100%

    * - .. _fig_pin_location_zoomed_out:
  
        .. figure:: /_static/hardware/ErLEED_modification/pin_location/location_medium.svg

            Filament output and highlighted location of new pin.

      - .. _fig_pin_location_zoomed_in:

        .. figure:: /_static/hardware/ErLEED_modification/pin_location/location_large.svg

            Zoomed in filament output and marked solder spots.


When finished, the pin should look as shown in :numref:`fig_pin_soldered`.
Test the connections using a multimeter.
Finally, connect the new pin to the cable leading to the new port.


.. _fig_pin_soldered:
.. figure:: /_static/hardware/ErLEED_modification/pin_location/connector_soldered.svg
    :width: 50%
    :align: center

    Beam HV pin soldered to the board.

**TODO Michael, Michele: Quick test to see if working as intended??**

Reassembly
==========

At this point, the modification is complete and you can reassemble the control unit.
Plug in all connectors on the motherboard and make sure they are securely connected, including the newly placed beam HV pin.
We also highly recommend you replace all cable ties that you cut during disassembly.
Finally, screw the back plate and the cover plate back on.
