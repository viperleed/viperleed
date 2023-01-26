.. _erleed_modification:

#########################
ErLEED Modification Guide
#########################
The :term:`ErLEED` LEED electronics produced by :term:`SPECS` are commonly used in many surface science laboratories.
The ViPErLEED electronics are designed (and tested) to work with the ErLEED electronics.
**TODO Michele, Michael, Alex: Details on which versions numbers are supported; why we perform the modifications, warnings etc.**


Opening up the electronics
=======================


.. _fig_connector_soldered:
.. figure:: /_static/hardware/ErLEED_modification/electronics_overview.svg
    :width: 75%
    :align: center

    LEED control unit with cover plate removed.

Removing the back plate
=======================

- cable ties

# image of cable ties

New port
========

**TODO Michele, Michael: names/numbers of HV connector components & hole-punch machine**

To read out the beam voltage during operation, we need access to the new pin from the outside.
For this, we will need to fashion a new port on the control unit cassis.
Fortunately, there is ample space on the backplate, next to the existing connectors.
We recommend placing the "Beam HV" port below the ":math:`I0_{\text{MON}}`" port, as shown in :numref:`fig_new_port_location`.


.. _fig_new_port_location:
.. figure:: /_static/hardware/ErLEED_modification/new_port_location.svg
    :width: 75%
    :align: center

    Components ...


Making the beam HV connector
============================

In the following steps, you will build and place a simple pin connector which allows reading out the beam potential.

You will need two 330 :math:`\Omega` resistors, a (?) and a short shrink tube as shown in :numref:`fig_resistors_1`.

.. _fig_resistors_1:
.. figure:: /_static/hardware/ErLEED_modification/resistors/resistors_1.png
    :width: 20%
    :align: center

    Components ...

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

.. _fig_resistors_cable_attached:
.. figure:: /_static/hardware/ErLEED_modification/resistors_cable_attached.jpeg
    :width: 25%
    :align: center

    Components ...



Placing the beam current pin
============================

.. list-table::
    :align: center
    :width: 100%

    * - .. _fig_resistors_2:
  
        .. figure:: /_static/hardware/ErLEED_modification/pin_location/location_medium.svg

            Bending resistor wire around the pin.

      - .. _fig_resistors_3:

        .. figure:: /_static/hardware/ErLEED_modification/pin_location/location_large.svg

            Two resistors soldered to the pin.



.. _fig_connector_soldered:
.. figure:: /_static/hardware/ErLEED_modification/pin_location/connector_soldered.svg
    :width: 50%
    :align: center

    Components ...

Reassambeling
=============

