.. _erleed_modification:

#########################
ErLEED Modification Guide
#########################

The :term:`ErLEED` LEED electronics produced by :term:`SPECS` are commonly used in many surface science laboratories.
The ViPErLEED electronics are designed (and tested) to work with the ErLEED electronics.
** TODO Michele, Michael, Alex: Details on which version numbers are supported; why we perform the modifications, etc. **


Opening up the electronics
=======================

# image of electronics inside

Removing the back plate
=======================

- cable ties

# image of cable ties

New port
========


First, we will make a new plug reading the beam current.


Beam current port
=================

You will need two 330 :math:`\Omega` resistors, a (?) and a short shrink tube as shown in :numref:`fig_resistors_1`.

.. _fig_resistors_1:
.. figure:: /_static/hardware/ErLEED_modification/resistors/resistors_1.jpeg
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
  
        .. figure:: /_static/hardware/ErLEED_modification/resistors/resistors_2.jpeg

            Bending resistor wire around the pin.

      - .. _fig_resistors_3:

        .. figure:: /_static/hardware/ErLEED_modification/resistors/resistors_3.jpeg

            Two resistors soldered to the pin.

      - .. _fig_resistors_4:

        .. figure:: /_static/hardware/ErLEED_modification/resistors/resistors_4.jpeg

            Protruding wires removed.

      - .. _fig_resistors_5:

        .. figure:: /_static/hardware/ErLEED_modification/resistors/resistors_5.jpeg

            Finished connector with shrink tube.

.. _fig_resistors_cable_attached:
.. figure:: /_static/hardware/ErLEED_modification/resistors_cable_attached.jpeg
    :width: 25%
    :align: center

    Components ...

.. _fig_connector_soldered:
.. figure:: /_static/hardware/ErLEED_modification/connector_soldered.jpeg
    :width: 25%
    :align: center

    Components ...

Placing beam current pin
========================

# zoomed in images


Reassambeling
=============

