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
    The ViPErLEED developers take no responsibility for any malfunctions that may occur as a result of modifications to the LEED electronics.

In order to make the beam voltage accessible, we need to measure the average potential at the filament.
We can do this by adding a voltage divider with two equal resistors in parallel to the filament.
The output of the voltage divider is then connected to a new port on the back of the control unit.
A schematic of the circuit is shown in :numref:`fig_ebeam_circuit`.

.. _fig_ebeam_circuit:
.. figure:: /_static/hardware/ErLEED_modification/upgrade_ebeam_circuit.svg
    :align: center

    Schematic of the modifications to the ErLEED electronics to enable beam voltage read-out.

:math:`F+` and :math:`F-` are the potentials at either end of the filament :math:`R_{\mathrm{fil}}`.
They relate to the real electron energy :math:`E` and the filament voltage :math:`\Delta V_{\mathrm{fil}}` as

.. math::
    F+ = E + \frac{\Delta V_{\mathrm{fil}}}{2}, \quad F- = E - \frac{\Delta V_{\mathrm{fil}}}{2}.

.. important::

    Always disconnect the BEAM HV plug when measuring the beam current :math:`I_0`.
    Otherwise, the current drawn by the voltage measurement will distort the current measurement.
    The current from the ViPErLEED interface electronics :math:`I_{\mathrm{ViPErLEED}}` is :math:`\approx6.25\,\mu\mathrm{A}/100` eV.

Required Components
===================

Before you start with the modifications, make sure you have all required components and tools.


For the beam HV pin:
    - 2 equal resistors (between 200 and 500 :math:`\Omega`) with ~1 W power rating, tolerance :math:`\leq 1\%`,
    - male and female plug pins (e.g., Vogt part no. 1365a.61 and 1361.61),
    - a short shrink tube (~1.5 cm),

For the beam HV port:
    - a SHV plug (e.g., RS stock no. 212-7444),
    - **TODO: cable length, diameter?, other components**


Additionally, you will need the following tools to perform the modifications:
    - a suitable anti-static electronics workbench,
    - a suitable set of Phillips, Pozi, and flat-head screwdrivers,
    - a soldering iron and solder,
    - tweezers,
    - cable clippers,
    - a cable tie cutter,
    - a heat gun for shrink tubes,
    - a sheet metal hole punch machine or a drill with a suitable drill bit.

.. _fig_resistors_1:
.. figure:: /_static/hardware/ErLEED_modification/resistors/resistors_1.png
    :width: 20%
    :align: center

    Components for the beam HV pin.

The two resistors should be of equal value and have a large power rating but they do not need to be precision resistors.
The measured beam voltage relates to the real electron energy :math:`E` as

.. math::
    \text{BEAM HV} = E \left(1 + \frac{R}{2R_{\mathrm{ViPErLEED}}}\right) + \frac{\epsilon}{4}\Delta V_{\mathrm{fil}}.

Here, :math:`\epsilon` is the relative difference between the two :math:`R` resistors in the voltage divider (see :numref:`fig_ebeam_circuit`).
Thus a 1% difference in the resistors at 20 V filament voltage will result in a 50 meV offset in the energy calibration.

Furthermore, the smaller :math:`R` is relative to :math:`R_{\mathrm{ViPErLEED}}`, the smaller the gain error of the energy calibration will be.
With 330 :math:`\Omega` resistors, the gain error is about :math:`1\times 10^{-5}`, or ~0.1 eV at 1000 eV.
However, the resistors should not be chosen too small as the additional current :math:`I_{\mathrm{extra}}` (see :numref:`fig_ebeam_circuit`) drawn by the voltage divider is

.. math::
    I_{\mathrm{extra}} = \frac{\Delta V_{\mathrm{fil}}}{2R} = I_{\mathrm{fil}} \frac{R_{\mathrm{fil}}}{2R}

when the BEAM HV plug is disconnected.
This equates to about 30 mA with 20 V filament voltage and a 330 :math:`\Omega` resistor, i.e. ~0.6 W.

Opening up the electronics
==========================

.. important::
    Before starting the modifications described below, turn off and completely disconnect all plugs from the control unit.
    **Disconnect the mains (power) cable** and wait 5 min before proceeding with opening the unit to allow all capacitors to fully discharge.


To start, we need to open up the electronics.
First, fully disconnect all plugs from the unit and place it on a suitable anti-static electronics workbench.
Start to open up the electronics unit by removing all screws that hold the cover plate in place, then remove the cover plate. The cover may have a protective-earth connector to be disconnected.
:numref:`fig_cover_plate_removed` shows the inside of the ErLEED control unit with the cover plate removed.


.. _fig_cover_plate_removed:
.. figure:: /_static/hardware/ErLEED_modification/electronics_overview.svg
    :width: 75%
    :align: center

    ErLEED control unit with cover plate removed.

Removing the back plate
=======================

To make the required modifications, it is also necessary to partially take off the back plate of the control unit.
There is no need to completely disconnect the back plate from the rest of the electronics, but fashioning the new port, as described below will likely require tilting the plate horizontally.

For best accessibility, you most likely want to remove the mains plug and screen connection from the back plate (compare :numref:`fig_new_port_location`).
Additionally, you may need to remove a number of cable ties and unplug various connections to the motherboard (e.g., connections for "ANODE", "L2", "L1/3", "FILAMENT", etc.).
We highly recommend labeling each cable before removal.

Ultimately, you should be able to remove and tilt back the back plate of the unit as shown in :numref:`fig_taking_off_backplate`.


.. _fig_taking_off_backplate:
.. figure:: /_static/hardware/ErLEED_modification/taking_off_backplate.svg
    :width: 75%
    :align: center

    Removing the back plate of the control unit.


.. _section_beam_hv_pin:

Beam HV pin
===========

In the following steps, you will build and place a simple pin connector which allows reading out the beam potential.
A suitable pin can easily be improvised from two 330 :math:`\Omega` resistors, a male plug pin and a short shrink tube as shown in :numref:`fig_resistors_1`.


First, bend the wires of the resistors around the male plug pin as shown in :numref:`fig_resistors_1`.
Then, place the pin in the loop and solder it in place. Repeat the procedure with the second resistor (see :numref:`fig_pin_from_resistors` b and c).
Make sure the soldered spot is stable and the connection across both resistors and to the pin is intact.
Finally, cut off protruding wires (see :numref:`fig_pin_from_resistors` c) from the resistors and place a shrink tube over the resistors and soldered spot (see :numref:`fig_pin_from_resistors` d).
It is a good idea to check at this stage that the two resistors still have approximately the same resistance as before. It is sometimes possible to overheat the resistors while soldering or applying a shrinking tube.

.. _fig_pin_from_resistors:
.. figure:: /_static/hardware/ErLEED_modification/resistors/pin_from_resistors.svg
    :width: 100%
    :align: center

    Steps to fashion the new pin. (a) Bending resistor wire around the pin. (b) Two resistors soldered to the pin. (c) Protruding wires removed. (d) Finished connector with shrink tube.


To read out the beam voltage during operation, we need access to the new pin from the outside.
For this, there needs to be a connection from the pin on the motherboard to the backside of the chassis.
To make the connection, take the new cable, strip off ~0.5 cm of the insulation on both side, and solder the female pin plug to the cable.
You should be able to firmly connect the cable to the pin as shown in :numref:`fig_resistors_cable_attached`.

.. _fig_resistors_cable_attached:
.. figure:: /_static/hardware/ErLEED_modification/resistors_cable_attached.jpeg
    :width: 25%
    :align: center

    New beam HV pin with connection cable attached.


New port
========


Next, we will need to fashion a new port on the chassis of the control unit.
Fortunately, there is ample space on the backplate, next to the existing connectors.
We recommend placing the "Beam HV" port below the ":math:`I0_{\text{MON}}`" port, as shown in :numref:`fig_new_port_location`.

Carefully punch (or drill) a hole in the backplate of the chassis in the desired location.
This may require removing the mains plug and screen connector from the backplate.
Then, place the new SHV connector, solder the connection wire to it and place a shrink tube over the solder spot.
Secure the coaxial connector in place by tightening the nut that came with it on the inside of the backplate.

Finally, re-attach the mains and screen plugs if you had to remove them.
At this point, we highly recommend labeling the newly fashioned port appropriately (e.g., "Beam HV").

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
By placing a voltage divider with two equal resistors parallel to the filament, we can get a reference potential (see the schematic circuit diagram in :numref:`fig_ebeam_circuit`).
Since the filament is essentially a short circuit, a negligible current (~20 mA) will pass through the parallel resistors.

To directly access the filament, we can conveniently place the new pin right next to the high-voltage diodes next to the filament port.
The exact location is shown in :numref:`fig_pin_location_zoomed_out` and :numref:`fig_pin_location_zoomed_in`.

Using a long needle, form two hooks with the remaining wires of the resistors used to make the beam HV pin in step :ref:`section_beam_hv_pin`.
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
Make sure to connect back the protective-earth connector to the chassis.
We also highly recommend you replace all cable ties that you cut during disassembly.
Finally, screw the back plate and the cover plate back on.
