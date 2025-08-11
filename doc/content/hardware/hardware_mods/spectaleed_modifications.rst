.. _spectaleed_modification:

#############################################
Omicron SpectaLEED/NG LEED Modification Guide
#############################################


The :term:`SpectaLEED<SpectaLEED>` LEED/Auger optics produced by :term:`Omicron<Scienta Omicron>` is among the most common optics used in surface science laboratories. It is controlled by the NG LEED (or NG LEED S) electronics. The ViPErLEED electronics are designed (and tested) to work with the NG LEED (S) electronics.
At variance with the :term:`ErLEED` electronics, the NG LEED (S) requires deeper modifications of the circuitry to allow reasonable LEED-:math:`I(V)` measurements. In the following sections, five main categories of modifications are described:

#. Modifications that are necessary for a :ref:`reliable measurement<ngleed_i0_mod_mandatory>` of :math:`I_0`
#. Modifications aimed at :ref:`reducing impact of noise<ngleed_improve_noise>` on :math:`I_0`
#. Modifications aimed at :ref:`improving the time response<ngleed_faster_beam>` of the beam energy, and, as a consequence, of :math:`I_0`
#. Modification of the :ref:`range<ngleed_suppressor_range>` of the ``SUPPRESSOR`` voltage
#. Modifications needed to :ref:`access the beam energy at high voltage<ngleed_beamhv_port>` (only needed for older units)

In addition, you will also need to :ref:`prepare a control cable in order to allow ViPErLEED to set the LEED energy<ngleed_control_cable>`.

.. note::
    The modifications described in :secref:`ngleed_i0_mod_mandatory` are **necessary** in order to be able to measure the beam current :math:`I_0` with a Omicron NG LEED (S) electronics. Without the modifications in :secref:`ngleed_i0_mod_mandatory` the :math:`I_0` signal that can be measured at the rear panel is **meaningless**.

The modifications described in :secref:`ngleed_improve_noise` and :secref:`ngleed_faster_beam` are optional, but **strongly recommended** for users with a microchannel-plate (**MCP**) LEED optics, as the small beam currents are strongly affected by noise (:secref:`ngleed_improve_noise`) and the slow response times whenever the beam energy is ramped (:secref:`ngleed_faster_beam`).

The modification described in :secref:`ngleed_suppressor_range` are **strongly suggested** for all users.

.. warning::
    The modifications described below **will** void any warranty unless you explicitly get a permission by the supplier.
    The ViPErLEED developers take no responsibility for any malfunctions that may occur as a result of modifications to the LEED electronics.

.. warning::
    Some of the modifications described below require overriding safety features of the NG LEED unit. This means that the unit will be connected to the mains voltage and will be allowed to produce high-voltages. In some cases, the high voltages should be directly measured. The ViPErLEED developers take no responsibility for any damage to people and/or equipment that may result from inappropriate application of safety measures.


Overview of NG LEED unit
========================

.. admonition:: TL;DR

    To open the LEED unit and access the electronics, follow these steps:

        * Make sure to disconnect the unit from the mains voltage and wait 5 min for the capacitors to discharge.
        * Open the chassis as shown in :numref:`fig_ngleed_new_overview`.
        * Remove the desired module(s) as shown in :numref:`fig_ngleed_remove_hv_cables`.
        * If necessary, remove the rear panel as shown in :numref:`fig_ngleed_rear_panel_tilted`.

:numref:`fig_ngleed_old_overview` and :numref:`fig_ngleed_new_overview` show an overview of the **older** and **newer** versions of the NG LEED electronics, respectively. The two versions can be told apart from the different appearance of the rear \[:numref:`fig_ngleed_old_overview`\ (a), :numref:`fig_ngleed_new_overview`\ (a)\] and front \[:numref:`fig_ngleed_old_overview`\ (c), :numref:`fig_ngleed_new_overview`\ (c)\] panels.

.. _fig_ngleed_old_overview:
.. figure:: /_static/hardware/SpectaLEED_modification/overview_ngleed_old.svg
    :align: center

    Overview of the older version of the NG LEED electronics with key instructions on how to open the chassis and how to remove modules. (a) View of the rear panel. **TODO** (b) Removal of side covers. (c) Removal of top and bottom cover grids. (d, e) Top and bottom views of the inside of the NG LEED electronics.

.. _fig_ngleed_new_overview:
.. figure:: /_static/hardware/SpectaLEED_modification/overview_ngleed_new.svg
    :align: center

    Overview of the newer version of the NG LEED electronics with key instructions on how to open the chassis and how to remove modules. (a) Partial view of the rear panel. (b) Removal of side covers. (c) Removal of top and bottom cover grids. (d, e) Top and bottom views of the inside of the NG LEED electronics.

The rear panel \[:numref:`fig_ngleed_old_overview`\ (a), :numref:`fig_ngleed_new_overview`\ (a)\] houses:

* The high-voltage plug providing the voltages to the ``FILAMENT``, ``ANODE``, suppressor (sometimes referred to as ``RETARD`` inside the units), ``WEHNELT``, and lenses (``L1/3`` and ``L2``). All these voltages depend on the ``BEAM`` energy.
* The SHV port for the ``SCREEN`` voltage and, only in the newer version, the one for monitoring the beam energy at high voltages. See :secref:`ngleed_beamhv_port` for instructions on how to add this **important** port also on the older NG LEED units.
* The programming input switch and port; switching to "external" --- by pulling on the lever, then flipping it --- deactivates the "beam energy" potentiometer knob on the front panel. The energy is set externally via an analog input at the ``BEAM ENERGY PROGRAM INPUT`` connector (the 5-way circular connector, **not to be confused** with the ``PROGRAMMING INPUT`` BNC that is used only for Auger). This is how ViPErLEED can control the beam energy. See :secref:`ngleed_control_cable` for instructions on how to prepare a suitable control cable.
* The :math:`I_0` BNC. See :secref:`ngleed_i0_mod_mandatory` for the relevant modifications.
* Other connectors not relevant for LEED-:math:`I(V)` measurements. Notice that the :math:`E_0` monitor is **not** the beam energy at high voltage, but a version downscaled to 0--10 V. This port should **not** be used instead of the high-voltage version, as it can (and will) have non-zero offset and non-unity gain.

.. important::
    Before starting the modifications described below, turn off and completely disconnect all plugs from the control unit.
    **Disconnect the mains (power) cable** and wait 5 min before proceeding with opening the unit to allow all capacitors to fully discharge.

.. _ngleed_opening:

Opening the chassis
-------------------

All the modifications described in the following sections require direct access to (various) electronics boards in the unit. First, fully disconnect all plugs from the unit and place it on a suitable anti-static electronics workbench.

To open the chassis, remove the side covers by inserting a screwdriver under it and forcing the two velcro-like pads to come apart (a little force may be necessary) --- see :numref:`fig_ngleed_old_overview`\ (b) or :numref:`fig_ngleed_new_overview`\ (b). You can now access the eight screws holding in place the top and bottom chassis covers. Both covers need to be removed to gain access to the module boards. Remove the screws, then, with the help of a flat screwdriver, lift the grids up. You may experience a little resistance: they are held in place by six metal clips (three on the front, three on the back) that will snap out when enough force is applied --- see :numref:`fig_ngleed_old_overview`\ (c) or :numref:`fig_ngleed_new_overview`\ (c).

.. warning::
    The clips holding the cover grids are very sharp. Do not use your fingers to lift the grid. Also, the clips can be deformed a bit in the process. They can be bent back using a plier. **Do not** attempt to bend back the clips with your bare hands.

.. _ngleed_extract_hv_board:

Extracting one HV module
------------------------

The high-voltage modules are mounted on the motherboard of the electronics. See :numref:`fig_ngleed_old_overview`\ (d) or :numref:`fig_ngleed_new_overview`\ (d). To remove any of the high-voltage modules you will need to completely remove all the aluminium profile bars on the top side of the chassis (except for the one closest to the rear panel, on which the cables carrying the high voltages are secured). The bars can be removed by unscrewing the four bolts holding them in place from the sides. We suggest to mark the positions and directions of the bars before removing them for easier reassembly at the end.

.. note::
    The threads on the bars may wear out quite quickly. Consider cutting them again with a suitable thread cutter before assembling the bars back.

The modules are secured to the motherboard with two (grounded) bolts that can be removed from the bottom side of the unit. See :numref:`fig_ngleed_old_overview`\ (e) or :numref:`fig_ngleed_new_overview`\ (e). To extract the HV modules, it is also necessary to remove the high-voltage cables that carry their output voltage(s) toward the rear panel. :numref:`fig_ngleed_remove_hv_cables` shows an overview of the HV cables to be removed (left). It is **highly recommended** to label each cable before removing it.

The cables are plugged via spade connectors onto the boards. The mating connectors soldered on the boards are easily damaged. To minimize the risk, it is advisable to use an angled tweezer, as shown in the right panel of :numref:`fig_ngleed_remove_hv_cables`. Very lightly grab the spade on the board from the sides, and fit the the tips of the tweezers in between the board and the connector of the cable. Then, use the tweezers as a lever, by rotating them around the corners of the pincers, that are lying on the board.

.. _fig_ngleed_remove_hv_cables:
.. figure:: /_static/hardware/SpectaLEED_modification/removing_hv_cables.svg
    :align: center

    Cables carrying the high voltages generated in each of the HV modules (left), and hint on how to unplug them without damaging the mating spade connector that is soldered to the board (right).

Now that the holding screws, the bars, and the HV cables have been removed, the module(s) can be unplugged from the motherboard by pulling them out. In the process, it may happen that the long hex nut of the holding screw gets stuck on the adjacent board (marked with a circle in the left panel of :numref:`fig_ngleed_remove_hv_cables`). There is no need to unplug all boards: they can be extracted by wiggling the board to be removed and lightly tilting the neighbouring one.

.. _ngleed_remove_rear:

Removing the rear panel
-----------------------

Some of the modifications described below require access to the high- and low-voltage connectors mounted on the rear panel. It is sufficient to tilt the rear panel into a (close-to) horizontal position. There is no need to fully remove it.

The rear panel can be simply removed by loosening the six screws holding it on the chassis. In order to tilt the panel, one has to also loosen the tension on the high-voltage-carrying cables. For this, it is usually enough to remove the screws holding the cables attached to the aluminium profile bar closest to the rear panel (see :numref:`fig_ngleed_remove_hv_cables` and :numref:`fig_ngleed_rear_panel_tilted`). Possibly, also the cable tie holding the ribbon cable in place needs to be removed. For some of the modifications --- especially those on the HV plug --- it may be useful to remove the mains fuse box as well.

.. _fig_ngleed_rear_panel_tilted:
.. figure:: /_static/hardware/SpectaLEED_modification/rear_panel_tilted.svg
    :align: center

    The newer NG LEED electronics with its rear panel free to be tilted back. This gives access to the HV plug as well as the :math:`I_0` output.

.. _ngleed_i0_mod_mandatory:

Making :math:`I_0` measurements possible
========================================


.. admonition:: TL;DR

    The :math:`I_0` output provided by the unmodified LEED electronics is next to meaningless. Some modifications are necessary to make it usable:.

      * On the BEAM HV module (see :numref:`fig_ngleed_i0_beam_module_mod`):

        * Cut the trace and insert a 1 kΩ resistor between ``HV_GND`` and the inverting input of U4.
        * Remove trimmer R43.
        * Replace U4 with OPA627.
        * Add a 1 nF (or 10 nF) capacitor on the feedback of U3.

      * On the E0 BUFFER board (see :numref:`fig_ngleed_i0_control_and_buffer`):

        * Add two 10 nF capacitors in parallel to resistors R7 and R8.
        * Consider replacing R1 and R2 as described in :secref:`ngleed_i0_modify_e0_buffer`.

      * On the BNC output, ensure that the 1 µF capacitor is not connected (see :numref:`fig_ngleed_i0_rear`).

.. _i0_instability:

Overview – Why are modifications necessary
------------------------------------------

The NG LEED electronics provides an analog 0--10 V output BNC on the rear panel \[see :numref:`fig_ngleed_old_overview`\ (a) or :numref:`fig_ngleed_new_overview`\ (a)\] that is supposed to provide measurements for the :math:`I_0` current. 1 V on the output corresponds to 1 µA. The :math:`I_0` current is measured internally with circuitry present on the ``BEAM`` HV module. The relevant section of the circuit is shown in :numref:`fig_ngleed_i0_circuit_beam_module_original`.

.. _fig_ngleed_i0_circuit_beam_module_original:
.. figure:: /_static/hardware/SpectaLEED_modification/i0_circuit_original.svg
    :align: center

    The portion of the ``BEAM`` HV module circuit responsible of producing the measurement of :math:`I_0`.

The measurement of :math:`I_0` is performed in the following manner: all HV modules that generate voltages for the electron gun (i.e., ``BEAM``, ``FILAMENT``, ``ANODE``, ``WEHNELT``, ``L1/3``, and ``L2``) are floating. Their reference potential is ``HV_GND`` (at times also referred to as ``HV_GUARD`). This means that the net current drawn from ``HV_GND`` contains all the contributions of all the electron-gun currents, and, in particular, the total number of electrons that leave the electron gun (i.e., those emitted by the ``FILAMENT``, minus those collected back at ``ANODE``, ``WEHNELT``, and lenses).

In :numref:`fig_ngleed_i0_circuit_beam_module_original`, the operational amplifier U4 holds ``HV_GND`` at the same potential as ``GND`` via the feedback, and acts as an inverting transimpedance amplifier for the :math:`I_0` current (with a gain of −100 kΩ = −1 × 10⁵ V⁠/⁠A). Together with the inverting buffer around U3, this corresponds to a gain of 1 × 10⁵ V⁠/⁠A, or 0.1 V⁠/⁠µA. The output of U3 is then amplified by another factor of 10 (on the ``CONTROL`` board, see :numref:`fig_ngleed_i0_control_and_buffer`), giving the overall 1 V/µA gain mentioned before. (Two more unity-gain stages follow --- see :secref:`ngleed_i0_modify_e0_buffer`.)

.. _warn_swapped_components:
.. warning::
    In our instrument, the resistors R15 and R17 are **swapped** relative to those shown in the official circuit diagram. This means, that the feedback resistor of stage U4 is actually R15, while R17 is the feedback resistor of stage U3. The resistor values are as shown in :numref:`fig_ngleed_i0_circuit_beam_module_original`, so stage U4 has transimpedance gain of −10 kΩ (rather than the −100 kΩ in :numref:`fig_ngleed_i0_circuit_beam_module_original`) while the second stage has a gain of −10. While the overall gain remains equal to 0.1 V/µA, this has important consequences for the accurate identification of how the modifications described below should be done. Diodes D10 and D11 were also **swapped** on our board.

The circuit in :numref:`fig_ngleed_i0_circuit_beam_module_original` however has an important issue that we have overlooked in the simple analysis above: there is a huge purely-capacitive (1.1 µF) input impedance on the inverting input of U4. This, combined with the fact that the operational amplifier is not just an infinite-gain, zero-output-impedance ideal one, gives strong instabilities, as it reduces the phase margin. :numref:`fig_ngleed_bode_instability` shows the Bode diagram of the circuit around U4, split into the forward gain (i.e., the gain of the `LF411 <https://www.jameco.com/jameco/products/prodds/835500.pdf>`_ op-amp U4) and the inverse of the feedback gain. The point where the two amplitude curves intersect corresponds to :math:`|G_\mathrm{open\,loop}| = 1`. If the phase of :math:`G_\mathrm{open\,loop}` at this point is close to 180° the feedback is positive, and the circuit is unstable (see `phase margin <https://en.wikipedia.org/wiki/Phase_margin>`_).

.. _fig_ngleed_bode_instability:
.. figure:: /_static/hardware/SpectaLEED_modification/i0_bode_instability.svg
    :align: center

    Bode diagram of the forward (black) and inverse-feedback (orange) gain of the configuration of U3. When the 40 Ω output impedance of the LF411 is neglected (dashed orange), the circuit appears stable: the phase margin (at 500 kHz) is approximately 90°. However, the output impedance, combined with the large input capacitance, gives an additional pole at ~22 kHz in the feedback gain. This degrades the phase margin to ~10° (at 104 kHz), making the circuit unstable. Adding a 1 kΩ resistor on the inverting input fixes the problem (blue) by introducing an extra feedback zero at ~145 Hz, and by moving the problematic pole down to ~835 Hz.

Neglecting the 40 Ω output impedance of the LF411, the inverse of the feedback gain follows the dashed orange line. This would suggest that the circuit is stable: the two curves intersect at 500 kHz with a 20 dB/dec slope difference and a 90° phase margin. However, the output impedance, combined with the large input capacitance, adds an extra pole at ~22 kHz. This moves the :math:`|G_\mathrm{open\, loop}| = 1` frequency to ~104 kHz: there, the two curves intersect with a 40 dB/dec slope difference and the phase margin is reduced to ~10°. This means that the original circuit design is **unstable** and will provide **meaningless** :math:`I_0` values.

.. _fig_ngleed_i0_circuit_beam_module_mod:
.. figure:: /_static/hardware/SpectaLEED_modification/i0_circuit_beam_module_mod.svg
    :align: center

    Schematic representation of the modification(s) to be performed on the portion of the ``BEAM`` HV module circuit responsible of producing the measurement of :math:`I_0`. Adding a 1 kΩ resistor makes the circuit stable (see text); replacing the LF411 operational amplifier with an OPA627 improves offsets and noise. The better operational amplifier does not need offset trimming; adding an extra feedback capacitor on U3 improves noise filtering and stability.

This very poor design decision can be fixed quite easily as shown in :numref:`fig_ngleed_i0_circuit_beam_module_mod`. An additional 1 kΩ resistor is placed on the inverting input, making the input impedance a low-pass with cutoff frequency of ~145 Hz. The new resistor also dominates the problematic pole: it is in series with the output impedance of the U4 operational amplifier in open-loop conditions. This means that the pole frequency (i) does not depend any longer on the specific value of the output impedance of U4, and (ii) is shifted down to ~835 Hz. The resulting inverse feedback gain is shown in :numref:`fig_ngleed_bode_instability` in blue. The :math:`|G_\mathrm{open\,loop}| = 1` frequency is moved to the unity-gain bandwidth of the operational amplifier (~3 MHz), and the phase margin is increased to ~45°, making the circuit stable.

A marginal side effect of adding the resistor is that ``HV_GND`` will be slightly different from ``GND``: the virtual connection at the input of U4 holds the "right side" (cf. :numref:`fig_ngleed_i0_circuit_beam_module_mod`) of the resistor to ``GND``; its "left side", i.e., ``HV_GND``, is at :math:`1\,\mathrm{kΩ} \cdot I_0`. Considering that :math:`I_0` is mostly in the 1 µA range, this means that ``HV_GND`` differs from ``GND`` by a negligible 1 mV.

Required Components
-------------------

Before you start with the modifications, make sure you have all required components and tools:
    - one 1 kΩ resistor,
    - one `OPA627BP <https://www.ti.com/lit/ds/symlink/opa627.pdf>`_ operational amplifier (e.g., RS code 660-4355),
    - one 1 nF and two 10 nF non-polarized radial capacitors (e.g., ceramic) --- or three 10 nF non-polarized radial capacitors if your ``BEAM`` HV module has a 10 kΩ feedback resistor on U3 (see :numref:`fig_ngleed_i0_circuit_beam_module_mod` and the :ref:`warning <warn_swapped_components>` in the previous section),
    - a suitable anti-static electronics workbench,
    - a suitable set of Phillips, Pozi, and flat-head screwdrivers,
    - a soldering iron and solder,
    - means to remove solder from PCB vias --- e.g., a manual or (better) pneumatic solder pump,
    - tweezers (not necessary, but useful),
    - a sharp blade (e.g., a paper cutter),
    - a multimeter,
    - small cable ties.

.. _ngleed_i0_beam_module:

Modifying the ``BEAM`` HV module
--------------------------------

Open up the electronics as described in :secref:`ngleed_opening`, then extract the ``BEAM`` HV board following the instructions in :secref:`ngleed_extract_hv_board`. The top panel in :numref:`fig_ngleed_i0_beam_module_mod` shows the region of the ``BEAM`` high-voltage module where the circuit in :numref:`fig_ngleed_i0_circuit_beam_module_mod` can be found.

.. _fig_ngleed_i0_beam_module_mod:
.. figure:: /_static/hardware/SpectaLEED_modification/i0_beam_board_mod.svg
    :align: center

    Overview of the ``BEAM`` high-voltage module (top) indicating the area where the circuit for :math:`I_0` measurement is located. A close-up view of the same area is shown in the bottom panels.

In order to insert the 1 kΩ resistor needed for stabilization of the amplifier (see :secref:`i0_instability` and :numref:`fig_ngleed_i0_circuit_beam_module_mod` for details) you will need to **cut** the thick track on the back of the board that connects ``HV_GND`` to the inverting input of U4. Before you proceed, make sure you understand the circuit, and take especially note of the :ref:`warning <warn_swapped_components>` in the previous section: cross check which components are actually connected to ``HV_GND`` and which ones to the inverting input of U4 (pin 2).

To interrupt the track, cut two segments across it using the blade knife. Cut as deep as needed to pass through the copper layer of the track: you will normally need a few passes. Additionally, cut a bit the board next to the track and parallel to it in order to separate the insulation above the track from the surroundings. Then, with the help of the tip of the knife (or some other sharp tool), remove the section of the track  between the two cuts. You should be able to lift away both the copper and the insulation above it. Tweezers or a small nose plier can help stripping the track off. Cross-check with a multimeter that the whole track has been removed by measuring the resistance (and not just using the continuity beeper!). You can see the cut track in the bottom-left panel of :numref:`fig_ngleed_i0_beam_module_mod`.

.. _fig_bent_resistors:
.. figure:: /_static/hardware/SpectaLEED_modification/bending_resistors.svg
    :align: center

    Resistors prepared for soldering in place of the cut track between ``HV_GND`` and the inverting input of U4.

:numref:`fig_bent_resistors` shows how to prepare the 1 kΩ resistor (or, for that matter, any other component) for soldering in place. As pointed out in :numref:`fig_ngleed_i0_beam_module_mod`, you can decide to place the resistor on either the back or front sides of the board: in the former case, you can bend the resistor leads into two loops, and solder them around the two pads, as shown in the bottom-left panel of :numref:`fig_ngleed_i0_beam_module_mod`. When mounting it on the front, you can bend the leads into hooks, and solder them around the leads of components (see the schematic indication in the bottom-right panel of :numref:`fig_ngleed_i0_beam_module_mod`). When choosing where to place the resistor, consider also the additional modifications mentioned below.

Since you already have the ``BEAM`` HV module out, it is worth performing a few more modifications that significantly improve the quality of the measurement of :math:`I_0`. Refer to the schematics in :numref:`fig_ngleed_i0_circuit_beam_module_mod`.

Replace the somewhat basic LF411 op-amp with the much better performing OPA627(BP) op-amp. The `OPA627BP <https://www.ti.com/lit/ds/symlink/opa627.pdf>`_ has: larger DC gain, larger unity-gain bandwidth, a factor of 10 less noise, almost a factor of 10 better input offset voltage, and a factor of 50 better input bias and offset currents. The OPA627 op-amp is a drop-in replacement for the LF411 with the exact same pinout. Before extracting the LF411, make sure to clearly note down the direction of the chip, i.e., which pin is the first one. Replacing the LF411 with the OPA627 op-amp also means you can (and **should**!) get rid of the R43 trimmer. The trimmer should probably not have been there in the first place, according to the datasheets of both `LF411 <https://www.jameco.com/jameco/products/prodds/835500.pdf>`_ and `OPA627BP <https://www.ti.com/lit/ds/symlink/opa627.pdf>`_. Trimming of the offset of :math:`I_0` is performed with a dedicated circuit on the ``E0 BUFFER`` board. See :secref:`ngleed_i0_modify_e0_buffer` for further details.

As an additional precaution, as well as an improvement of the noise level on the :math:`I_0` output, install an extra capacitor in parallel to the feedback resistor of the second amplification stage (U3; see :numref:`fig_ngleed_i0_circuit_beam_module_mod`). This capacitor will improve the rejection of high-frequency interference signals, for example those that couple in from the switching power supply of the NG LEED. You can choose the value of the feedback capacitor for a ~1.6 kHz cutoff frequency. The value of the capacitor depends on the value of the feedback resistor of U3. As :ref:`mentioned earlier <warn_swapped_components>`, our instrument differed from the 'official' schematics: we have a 100 kΩ resistor rather than a 10 kΩ one (and the 10 kΩ is on the feedback of U4). 1.6 kHz cutoff corresponds to a 1 nF capacitor in parallel to 100 kΩ, or to 10 nF in parallel to 10 kΩ. You can solder the capacitor around the leads of the feedback resistor, after having shaped its leads as two hooks, similar to :numref:`fig_bent_resistors`. The bottom-right panel of :numref:`fig_ngleed_i0_beam_module_mod` shows the capacitor mounted in parallel to R15 which, :ref:`for our instrument <warn_swapped_components>`, acts as the feedback resistor of U3.

.. _fig_bode_i0_mod_overall:
.. figure:: /_static/hardware/SpectaLEED_modification/i0_bode_overall.svg
    :align: center

    Bode diagrams of the transimpedance gain of the circuit that measures :math:`I_0`, after the modifications described in this guide. The bode diagrams for the cumulative gain of several stages of amplification are shown. The first transimpedance stage around op-amp U4 (black) and the second voltage-gain stage around op-amp U3 (blue) are on the ``BEAM`` HV module. They are followed by a tenfold amplification on the ``CONTROL`` board (green), and by two more unity-gain stages (orange) on the ``E0 BUFFER`` board.

:numref:`fig_bode_i0_mod_overall` shows the overall transimpedance gain of the U4 amplification stage (black), as well as the one of the combination of the U3 and U4 stages (blue) after the modifications mentioned above. The plot considers the case in which stages U4 and U3 have −10 kΩ and −10 gains respectively. The extra capacitor added on the feedback of U3 maintains the unity-gain bandwidth of the combination of the U3 and U4 stages in the same range as the one of U4 (i.e., ~10 kHz).

.. note::

    While you have the ``BEAM`` HV module unmounted, consider looking also into the modifications described in :secref:`ngleed_faster_beam`. They are especially suggested for users with a microchannel-plate LEED.

.. _ngleed_i0_modify_e0_buffer:

Modifying the ``E0 BUFFER`` board
---------------------------------

The :math:`I_0` output of the ``BEAM`` HV module is further processed in two additional boards within the NG LEED unit. :numref:`fig_ngleed_i0_control_and_buffer` summarizes the location and functionality of the two extra processing stages.

.. _fig_ngleed_i0_control_and_buffer:
.. figure:: /_static/hardware/SpectaLEED_modification/i0_control_and_buffer_boards.svg
    :align: center

    Location, functionality and modifications of the two stages that process the :math:`I_0` voltage output produced on the ``BEAM`` HV module. Location of the (a) ``E0 BUFFER`` and (b) ``CONTROL`` boards within the chassis. (c) ``E0 BUFFER`` board removed for the modifications. (d) Circuit diagram of the relevant part of the ``E0 BUFFER`` board, including the suggested modifications.

The transimpedance-amplified :math:`I_0` signal generated on the ``BEAM`` HV module reaches the ``CONTROL`` board (via the motherboard) through the wide ribbon cable on the right side of :numref:`fig_ngleed_i0_control_and_buffer`\ (b). There, it is amplified by an additional factor of 10 with a non-inverting configuration around one of the op-amps in U24 \[see inset of :numref:`fig_ngleed_i0_control_and_buffer`\ (b)\]. As mentioned at the beginning of section :secref:`i0_instability`, this extra amplification stage is necessary because the gain on the ``BEAM`` HV module is 0.1 V/µA, while the unit is designed for an overall gain of 1 V/µA.

As can be seen in :numref:`fig_ngleed_i0_control_and_buffer`\ (b), the ``CONTROL`` board is found right behind the front panel of the NG LEED unit. We decided to not modify this amplification stage, but you can choose to add a capacitor in parallel to the 18 kΩ feedback resistor R131 in case you experience excessive noise. You should be able to solder it after removing the few cables around (i.e., the ribbon cable and the connector for the ``SUPPRESSOR`` setpoint): there should be no need to remove the whole board. In the Bode diagram of :numref:`fig_bode_i0_mod_overall`, the green curve includes the contribution of this amplification stage in the unmodified state. It is a pure gain stage with the expected pole at 100 kHz, given the ~1 MHz unity-gain bandwidth of the `AD704 <https://www.analog.com/media/en/technical-documentation/data-sheets/AD704.pdf>`_ op-amp.

The last processing stage occurs on the ``E0 BUFFER`` board \[visible in :numref:`fig_ngleed_i0_control_and_buffer`\ (a)\], mounted on the side panel of the chassis, and, unfortunately, very close to the switching power supply \[top left in :numref:`fig_ngleed_i0_control_and_buffer`\ (a)\] as well as the high-voltage supplies --- in the region of the large heat sink on the motherboard. The :math:`I_0` signal reaches the ``E0 BUFFER`` board through the thinner ribbon cable of the ``CONTROL`` board \[left in :numref:`fig_ngleed_i0_control_and_buffer`\ (b)\]. After the processing on the ``E0 BUFFER`` board, the :math:`I_0` signal travels along the long ribbon cable toward the output BNC on the rear panel. As the ribbon cables are unshielded, they can easily pick up high-frequency noise: they are not the ideal choice for cables passing next to the switching power supplies. You can improve this by twisting them around, as visible in :numref:`fig_ngleed_i0_control_and_buffer`\ (a), and by adding an iron core around (at least) the short one --- which passes right above the high-voltage transformer. You will need to untie the long ribbon cable from the support bars in order to twist it all the way toward the rear panel. Use cable ties to keep it together.

The primary role of the processing of :math:`I_0` in the ``E0 BUFFER`` board is offset compensation. Offsets arise because of non-idealities of the op-amps, but should be almost negligible after the modifications described in this guide, especially the replacement of the LF411 op-amp suggested in :secref:`ngleed_i0_beam_module`. The offset correction happens in the first (inverting) unity-gain stage of the ``E0 BUFFER`` board, as visible in :numref:`fig_ngleed_i0_control_and_buffer`\ (d). U1.B adds a correction current :math:`I_\mathrm{correction}` to the the :math:`I_0` signal coming from the ``CONTROL`` board. The circuit in the bottom part of :numref:`fig_ngleed_i0_control_and_buffer`\ (d) generates the correction from a 10 V reference:

.. math::
    I_\mathrm{correction} = \frac{10\,\mathrm{V}}{R_1} \left(1 - \frac{R_1}{R_2} x\right),

where :math:`x` is the fractional position of the trimmer R3. Using :math:`R_1 \approx 2 R_2` gives an (approximately) symmetric offset adjustment range of :math:`\pm 10\,\mathrm{V}/R_1`. The output of U1.B is inverted once more with the unity-gain stage around U1.A. The final output is

.. math:: V_{I_0} + \frac{R_8}{R_1}\left(1 - \frac{R_1}{R_2} x\right)\,10\,\mathrm{V} .

In terms of the original current :math:`I_0`, the signal is then

.. math:: I_0 + \frac{R_8}{R_1}\left(1 - \frac{R_1}{R_2} x\right)\,10\,\mathrm{µA} ,

where we have used the 1 V/µA gain between :math:`I_0` and :math:`V_{I_0}` resulting from the combined transimpedance amplification of the ``BEAM`` and ``CONTROL`` boards. By picking R1 relative to R8, one can then tune the range of variability of the offset correction. With the original values in :numref:`fig_ngleed_i0_control_and_buffer`\ (d), the range of offset correction is (−2.42, +2.13) µA, which is probably larger than any reasonable offset. It is a good idea to improve the range of offset correction by increasing the values of R1 and R2. Using :math:`R_1=680\,\mathrm{kΩ}` and :math:`R_2=330\,\mathrm{kΩ}` gives a more reasonable range of (−156, +147) nA. Users with a microchannel-plate LEED may want to use even larger resistors.

.. note::
    Before choosing resistors R1 and R2 that are appropriate for your unit, we suggest that you measure your offset to evaluate which range makes most sense. You should leave this as the **last step** among all the edits suggested in this guide. Leave the chassis open and connect the mains. Short the interlock pins on the outside of the HV plug (see :numref:`fig_ngleed_hv_plug` and :numref:`fig_ngleed_suppressor_check_display`) with a short wire, and turn on the unit. Wait for at least 30 min to allow for warm-up. **Leave the beam energy control knob at zero on the front panel**. Connect a multimeter to the ``I0`` BNC of the rear panel. Expect voltages in the millivolt range. Using a small screwdriver, turn the trimmer R3 of the ``E0 BUFFER`` board. **Be careful** as the mains supply cables run somewhat close by. If you feel like you would need more resolution to be able to trim the value to zero, you need larger resistor values for R1 and R2.

Aside from modifying the range of offset adjustment, you should also use the two unity-gain stages on the ``E0 BUFFER`` board to include some more filtering of the high-frequency noise (which may have been picked up by the ribbon cable, as mentioned above). To this end, solder 10 nF capacitors in parallel to the feedback resistors of both gain stages, as shown in :numref:`fig_ngleed_i0_control_and_buffer`\ (d). Use the hints in :numref:`fig_bent_resistors` to prepare the capacitors. Adding the two capacitors gives an extra second-order low-pass filtering --- with a cutoff frequency of ~1.6 kHz --- to the overall transimpedance gain. The Bode diagram of the overall gain resulting from this modification is drawn in :numref:`fig_bode_i0_mod_overall` as an orange trace.

.. _ngleed_i0_rear_panel:

Checking the BNC output
-----------------------

Another major design fault exists in the NG LEED unit. According to the official schematics, a 1 µF capacitor should be present at the output BNC between the signal (center) and ground (shell). This is visible in :numref:`fig_ngleed_i0_rear`. This is problematic, as the very large capacitor is essentially on the output of stage U1.A of the ``E0 BUFFER`` board. See schematics in :numref:`fig_ngleed_i0_control_and_buffer`\ (d). The problem is very similar to the one that causes instability of the first transimpedance stage on the ``BEAM`` HV module (solved in :secref:`ngleed_i0_beam_module`): the 200 Ω output impedance of the AD704 op-amp, in series with the capacitor, introduces a pole in the feedback factor that reduces the phase margin and can make U1.A unstable. The `datasheet of AD704 <https://www.analog.com/media/en/technical-documentation/data-sheets/AD704.pdf>`_ indeed indicates that the op-amp can drive at most a 10 nF capacitive load. In principle, the addition of the 10 nF capacitor in parallel to the feedback resistor of U1.A --- described in :secref:`ngleed_i0_modify_e0_buffer` --- should maintain the phase margin large enough for stability. Nevertheless, we advise to **remove** the 1 µF capacitor.

.. _fig_ngleed_i0_rear:
.. figure:: /_static/hardware/SpectaLEED_modification/i0_rear_panel.svg
    :align: center

    Location of the :math:`I_0` BNC output on the rear panel of the NG LEED unit.

For this purpose, open up the rear panel of the unit, as described in :secref:`ngleed_remove_rear`. The incriminated capacitor can be seen in the right panel of :numref:`fig_ngleed_i0_rear`. It is not necessary to fully remove the capacitor: it's enough to unsolder its lead that is connected to the center conductor of the BNC.

.. note ::
    On our unit, while the capacitor was present, there has clearly been a (lucky) manufacturing error: Both leads of the capacitor were soldered to the stainless steel wire connected to the shell of the BNC plug --- i.e., the capacitor was connected between ground and... ground. This meant that, for our unit, there was no need to unsolder the capacitor lead. Cross check that your unit indeed has the problem before unsoldering.

As you have the rear panel open, consider also the modifications described in :secref:`ngleed_rewire_hv_plug` and :secref:`ngleed_shield_mains`, especially suggested for users with a microchannel-plate LEED.

.. _ngleed_improve_noise:

Reducing noise on :math:`I_0`
=============================

.. admonition:: TL;DR

    These modification are recommended for users with a microchannel-plate LEED to reduce the noise on :math:`I_0`.

      * On the ``ANODE``, ``FILAMENT``, ``L1/3``, and ``L2`` boards rewire the guard ring from ``HV_GND`` to ``GND``, as shown in :numref:`fig_swap_hvguard_on_boards`.
      * Rewire the HV connector as shown in :numref:`fig_ngleed_hv_plug`.
      * Add a shield to the mains supply cable as shown in :numref:`fig_ngleed_shield_mains`.

This modification of the NG LEED unit is strongly suggested for users with a microchannel-plate LEED, where significantly lower electron currents are used (:math:`I_0 \approx 1-30\,\mathrm{nA}`). Users with a standard LEED will normally have beam currents in the microampere range and should most likely not need to modify their unit. The modifications described in this section should be considered a second-order improvement of those in :secref:`ngleed_i0_mod_mandatory`.

Once the modifications in :secref:`ngleed_i0_mod_mandatory` have been carried out (with the exception of the modification of the range of :math:`I_{0,\mathrm{offset}}` adjustment), the next-worst source of noise on :math:`I_0` has to do with the generation of the high voltages. A more detailed description of how high voltages are generated in the HV modules can be found in :secref:`ngleed_faster_beam`. In short, each high-voltage module generates its voltage with a `Voltage multiplier <https://en.wikipedia.org/wiki/Voltage_multiplier>`_ fed by a transformer. The transformer separates the 'high-voltage' from the 'low-voltage' areas of each module. :numref:`fig_ngleed_beam_high_and_low_voltage` shows, for example the ``BEAM`` board --- which, as mentioned below, is one of the few where no modification is needed.

.. _fig_ngleed_beam_high_and_low_voltage:
.. figure:: /_static/hardware/SpectaLEED_modification/beam_board_high_and_low_voltage.svg
    :align: center

    View of the high- and low-voltage portions of the ``BEAM`` HV module. All other HV modules are structured similarly.

The low-voltage area is also surrounded by a guard ring, i.e., a relatively thick track held at ground that shields the low-voltage from the high-voltage side. Several boards (``ANODE``, ``FILAMENT``, ``L1/3``, ``L2``, and ``WEHNELT``) have also an optoinsulator module allowing communication between the two sides of the board. On all boards except for ``ANODE``, the optoinsulator is located in a smaller PCB mounted perpendicular to the module (not shown). The guard ring also shields the low-voltage side of the optoinsulator parts.

Importantly, the guard ring is also connected to the low-voltage side of the electrostatic shield of the transformer. (The high-voltage side is also shielded separately.) The main purpose of shielding is to attenuate as much as possible the common-mode noise between the two sides of the transformer: the shield on each side picks up the noise from the respective winding. It is then important to make sure that each shield is connected to the **correct ground**. Unfortunately this is done **incorrectly** in most of the boards of the NG LEED unit. The low-voltage shield (and guard ring) on ``ANODE``, ``FILAMENT``, ``L1/3`` and ``L2`` HV modules is connected to ``HV_GND`` rather than to power ``GND``. This means that ``HV_GND`` picks up noise from the low-voltage windings of each of these transformers. As discussed in :secref:`i0_instability`, :math:`I_0` is measured as the total current flowing from ``HV_GND`` to ``GND``. The noise picked up by ``HV_GND`` is then present also on :math:`I_0`.

The solution is to **reconnect** the guard ring and transformer shields of all the improperly connected boards. Follow the instructions in :secref:`rewire_guard_rings`.

Other poor design decisions can be fixed by :secref:`ngleed_rewire_hv_plug`, where components that are supposed to be on ``HV_GND`` are on ``GND`` instead.

:secref:`ngleed_shield_mains` also helps reducing the noise on :math:`I_0` by adding a metallic plate between the mains input and the ribbon cable that carries :math:`I_0` to the rear panel.

.. todo::
    @Michael: I'm not quite sure which one of these is also supposed to help with the dielectric relaxation that I haven't mentioned yet.

Required Components
-------------------

Before you start with the modifications, make sure you have all required components and tools:

    - a suitable anti-static electronics workbench,
    - a suitable set of Phillips, Pozi, and flat-head screwdrivers,
    - a soldering iron and solder,
    - a sharp blade (e.g., a paper cutter),
    - a multimeter.

Additionally, for the modifications in :secref:`ngleed_rewire_hv_plug`:

    - two UF4002 diodes,
    - insulated wire (e.g., 0.25 mm²) and means to strip its insulation at the ends,
    - thin insulating sheet (e.g., an old transparency foil),
    - double-sided tape or other means to glue the sheet,
    - wire cutter.

Finally, for the modification in :secref:`ngleed_shield_mains`:

    - aluminium (or other high-conductivity material) sheet metal,
    - metal-working tools (e.g., file, saw, drill).

.. _rewire_guard_rings:

Rewiring the low-voltage guard
------------------------------

The following boards need rewiring and should be removed as described in :secref:`ngleed_extract_hv_board`: ``ANODE``, ``FILAMENT``, ``L1/3``, ``L2``. The ``WEHNELT`` module is also incorrectly wired, but there is no transformer on it, so it should not contribute to the noise on :math:`I_0`. As most boards look very similar to one another, we suggest to proceed with one board at a time in oder not to confuse which board is which.

:numref:`fig_swap_hvguard_on_boards` shows, for each board, suggestions of where the connections should be interrupted by **cutting** the relevant tracks and where they can be reconnected to the correct ground lines. For tips on how to cut tracks, see :secref:`ngleed_i0_beam_module`.

.. _fig_swap_hvguard_on_boards:
.. figure:: /_static/hardware/SpectaLEED_modification/swap_hvguard_boards.svg
    :align: center

    Suggestions on how to correctly rewire the guard ring from ``HV_GND`` to ``GND`` on the faulty boards.

On all boards, the track to be cut (i.e., ``HV_GND``) is the one connected to the **first** pair of pins of the low-voltage connector. The next pair of pins of the same connector is on ``GND``, i.e., where the guard ring and transformer shields should be connected to.

The easiest way to reconnect the guard ring is by stripping away a narrow bit of insulation from it (and, in some cases, also from the track to which to connect). It should be enough to scratch away the insulation by pressing strongly against the track with a small flat screwdriver, then sliding it along the track while applying pressure. Usually a single pass is enough. **Be careful not to damage the copper connection underneath the insulation**. It is now very easy to join the exposed copper of the stripped track to ``GND``. You can use the remainder of the leads of one of the discrete components that you have used in :secref:`ngleed_i0_mod_mandatory` and solder it in place. There is no need to add insulation, as the track is at ground. Use :numref:`fig_swap_hvguard_on_boards` and a multimeter to identify the correct spot.

.. _ngleed_rewire_hv_plug:

Rewiring the HV plug
--------------------

The capacitors visible on the inside of the high-voltage connector are connected between the high-voltage outputs and ``GND``. This is **incorrect** for all the voltages related to the electron gun (i.e., ``FILAMENT``, ``ANODE``, ``WEHNELT``, ``L1/3``, and ``L2``) which should be referred to ``HV_GND`` instead. This affects :math:`I_0`, which, as described :ref:`earlier <ngleed_i0_mod_mandatory>`, is the total current draw between ``HV_GND`` and ``GND``.

The high-voltage connector can be accessed by :ref:`releasing the rear panel <ngleed_remove_rear>`. :numref:`fig_ngleed_hv_plug` shows the HV connector after the modifications described below, as well as a schematic of the wiring diagram.

.. _fig_ngleed_hv_plug:
.. figure:: /_static/hardware/SpectaLEED_modification/hv_plug.svg
    :align: center

    Modification of the wiring of the high-voltage connector on the rear panel of the NG LEED. The stainless steel wire serving as the ground for all eight 1 nF capacitors should be split in two. The bottom part, exclusively connected to ``RETARD``, should maintain the current connection to protective earth (i.e., ``GND``). The top part where the other electron-gun-related voltages are connected, should be moved to ``HV_GND``. Two antiparallel diodes can be used as a transient-voltage suppressor between the two grounds.

The 1 nF capacitors are ``GND``\ ed together via the stainless-steel wire surrounding the high-voltage connector. You will **not unsolder** the capacitors, but will need to **cut open** the wire to move the top six capacitors to ``HV_GND``. The wire should be cut midway between the bottom two capacitors on both sides of the high-voltage connector, as indicated in the schematics in :numref:`fig_ngleed_hv_plug`.

The newly created four ends of the wire can be bent outwards in order to install two antiparallel diodes between the two grounds. The diodes will act as a transient-voltage suppressor. They are not strictly needed, but useful as they also make the assembly sturdier. They can be prepared as the 'flat-mounted' resistors in :numref:`fig_bent_resistors`: the cut-and-bent ends of the stainless-steel wire can be inserted in the hoops and securely soldered in place. As you solder the cathode of diode D2, also add an insulated wire in the top hoop. You can then connect the other end of this wire to ``HV_GND``: use the rightmost small pin at the very bottom of the high-voltage connector, labelled '1F' in the schematics of :numref:`fig_ngleed_hv_plug`.

As a result of cutting the wire, the top part of the construction is only supported by the six capacitors. The bottom part, instead, is soldered to (and grounded by) the slug protruding from the shell at the bottom of the connector. To prevent inadvertent contact between the wire and the chassis --- which would render :math:`I_0` measurements impossible --- it is a good idea to (i) add an extra **insulated** support, and (ii) glue (e.g., with double-sided tape) a thin insulating sheet underneath the ``HV_GND``-connected portion of the wire. For the support, you can solder a short piece of insulated wire to one of the the solder slugs at the top of the connector (see the top-left panel of :numref:`fig_ngleed_hv_plug`). Before you glue the insulating foil, consider the modification in :secref:`ngleed_shield_mains`, as the mains plug is quite close.

.. _ngleed_shield_mains:

Shielding the mains plug
------------------------

.. todo:
    Check if this is also needed for the older units. Looks like the HV plug is far from the mains there. It may still be necessary to shield the I0.

The positioning of the mains (i.e., power) fuse box and cables in the NG LEED is a bit unfortunate: it is right next to the high-voltage output as well as to the the ``HV_GND`` connection cable (black in our unit, in the foreground of the center-top panel of :numref:`fig_ngleed_hv_plug`). This means that there can be significant capacitive coupling between the mains and the high-voltage cables, in turn showing up as noise at the line frequency. While our ViPErLEED hardware has a very effective suppression of the line frequency **TODO: ref the section about the filter of the AD7705 in our box**, it is a good idea to minimize the noise in the first place.

.. _fig_ngleed_shield_mains:
.. figure:: /_static/hardware/SpectaLEED_modification/mains_shielding.svg
    :align: center

    Addition of a metallic shielding surrounding the mains fuse box and cables. Dimensions in millimetres.

You can improve the shielding by adding a simple grounded metal plate between the mains fuse box and the high-voltage connector, as shown in :numref:`fig_ngleed_shield_mains`. An aluminium plate bent into an 'L' shape with the rough dimensions in the top panel of :numref:`fig_ngleed_shield_mains` should fit in between the fuse box and the high-voltage connector. It can be held in place using the nuts of two existing screws --- the rightmost fuse-box mounting screw, and an (unused) hex stud next to the high-voltage connector ---, as indicated in the bottom-right panel of :numref:`fig_ngleed_shield_mains`.

It is a good idea to drill/file the mounting holes with a bit of clearance. This is because:

    - The space between the fuse box and the HV connector is quite tight;
    - There are two protective earth (yellow--green) wires and one ``HV_GND`` cable that should reach the high-voltage connector passing on top of the motherboard edge. You will need to leave ~2 mm clearance between the shield plate and the motherboard to avoid damaging the wires.

After you have mounted the plate, make sure you take precautions to prevent contact of the stainless steel wire --- modified in :secref:`ngleed_rewire_hv_plug` --- with the plate. For example, add a small pad of insulating sheet under diode D2 in :numref:`fig_ngleed_hv_plug`.

.. _ngleed_faster_beam:

Improving the beam-energy time response
=======================================

.. admonition:: TL;DR

    These modifications are **not** strictly necessary, but they can improve the time response of the beam energy and :math:`I_0`. They are **strongly suggested** for MCP-LEED setups. This enables faster LEED-:math:`I(V)` measurements:

      * On the ``BEAM HV`` board replace capacitors C32 and C40 with two 470 nF capacitors and C23 with a 470 µF electrolytic capacitor.

The speed at which the beam energy is changed from one value to the next, and, particularly, the time it takes to stabilize a new value of the energy determines how quick a LEED-:math:`I(V)` measurement can be. In fact, a LEED optics is primarily a capacitive load for the controlling electronics: a step in energy requires adjusting the voltages of (at least) filament and lenses accordingly by (dis)charging them. Also, the shielded cables carrying the voltages to the optics are a primarily capacitive load. In turn, this means that the effect of stepping the energy has even more important consequences on the stabilization of the :math:`I_0` current, to which all (dis)charging currents contribute. (The current through a capacitive load is proportional to the time derivative of the voltage across it.)

.. _fig_ngleed_hv_i0_relaxation:
.. figure:: /_static/hardware/SpectaLEED_modification/hv_i0_time_response.svg
    :align: center

    Evolution of (top) beam energy and (bottom) :math:`I_0` current after a large (100 eV) step of the beam energy.

:numref:`fig_ngleed_hv_i0_relaxation` shows the evolution of the beam energy and of the :math:`I_0` current after a 100 eV step of the beam energy for an NG LEED electronics before the modifications described in the present section. The real beam energy overshoots significantly, then undershoots, and wiggles its way toward the final stable value. Very small oscillations in the beam energy are still visible after 500 ms. The impact on :math:`I_0` is dramatic, with strong oscillations in correspondence of the wiggles visible on the beam energy trace. It is important to note that the height of the step was exaggerated on purpose. Typical energy steps used for LEED-:math:`I(V)` measurements are 0.5 eV. However, measurements with smaller steps (not shown) reveal that the amount of swing is roughly proportional to the step height. With smaller energy steps, :math:`I_0` does not saturate in the first 250 ms but shows a pronounced undershoot.

.. note::
    Users with a standard LEED setup likely need not worry about the response of :math:`I_0`. With 0.5 eV steps, the maximum peak-to-peak swing is of the order of 50 nA, much smaller than the typical microampere-range currents used. Instead, the amount of swing is **problematic** for microchannel-plate LEED setups, where typical :math:`I_0` values are in the 1--20 nA range.

The modifications in this section aim at mitigating the effects visible in :numref:`fig_ngleed_hv_i0_relaxation` by modifying the ``BEAM`` HV module. It should be mentioned that, at the time of writing, **we have not yet managed to completely remove the transients** of :numref:`fig_ngleed_hv_i0_relaxation`. This currently **limits the speed** at which LEED-:math:`I(V)` measurements can be acquired on an Omicron SpectaLEED instrument (with MCP): at each energy step ~600 ms are necessary for an acceptable stabilization of :math:`I_0`.

.. _fig_ngleed_beam_board_flow:
.. figure:: /_static/hardware/SpectaLEED_modification/beam_board_circuit_flow.svg
    :align: center

    Schematic representation of the location (top) and interconnection (bottom) of the various functional blocks present on the ``BEAM`` HV module responsible for the production of the beam energy of the electrons. The modifications described in the current section only target the blocks marked with a star, i.e., the PI regulator and the supply for the high-voltage transformer.

:numref:`fig_ngleed_beam_board_flow` shows a schematic view of the ``BEAM`` HV module, indicating the various functional blocks that participate in producing the desired electron beam energy. The high-voltage output of the board is generated using a transformer and a three-stage `high-voltage cascade-multiplier <https://en.wikipedia.org/wiki/Voltage_multiplier>`_, parts of which are visible on the left side of :numref:`fig_ngleed_beam_board_flow`. The high-voltage signal is divided by a factor 100, lowpass-filtered (cutoff at ~340 Hz) and fed back to the circuitry controlling the transformer. The control of the high-voltage output is obtained via a proportional--integral (PI) regulator that, given the desired energy value and the (downscaled and filtered) output, generates a control signal for the supply of the transformer.

.. note::
    All the high-voltage and high-current modules in the NG LEED control unit use the same principle as in the bottom part of :numref:`fig_ngleed_beam_board_flow`: the output of a PI regulator is used as the control signal for the driver of a transformer that generates the high voltage/current. The output value is fed back to the PI regulator for control. The various modules conceptually only differ in the circuitry that generates the specific output from the transformer.

In order to improve the time response of the beam energy, we will modify the PI regulator as well as the circuit for the supply of the transformer.

Before you start with the modifications, make sure you have all required components and tools:
    - the ``BEAM`` HV module, removed following the instructions in :secref:`ngleed_extract_hv_board`,
    - two 470 nF capacitors (e.g., film capacitors),
    - one 470 µF electrolytic capacitor rated for at least 25 V and with **large maximum ripple current** (e.g., RS 725-8168).
    - a suitable anti-static electronics workbench,
    - a suitable set of Phillips, Pozi, and flat-head screwdrivers,
    - a soldering iron and solder,
    - means to remove solder from PCB vias --- e.g., a manual or (better) pneumatic solder pump,
    - tweezers (not necessary, but useful).

.. _fig_ngleed_beam_board_response_time:
.. figure:: /_static/hardware/SpectaLEED_modification/beam_board_response_time.svg
    :align: center

    Top panels: position of the components that should be replaced in the PI regulator circuit (blue) and in the supply for the HV transformer (white). The insets show detailed view of the same components. Bottom: schematic circuit diagrams of the two functional blocks.

:numref:`fig_ngleed_beam_board_response_time` can be used as a guide to identify the location of the components that will be replaced. As can be seen in the bottom-left part of :numref:`fig_ngleed_beam_board_response_time`, the PI regulator takes the desired value of the beam energy (BEAM_SET) and subtracts the down-scaled version of the current energy (BEAM_MON) to generate the control voltage

.. math:: V_\mathrm{PI} = - K_\mathrm{p} e(t) - \frac{1}{T_\mathrm{int}}\int_0^t e(\tau) \mathrm{d}\tau.

The error signal is given by

.. math:: e = \mathrm{BEAM\_SET} - \mathrm{BEAM\_MON}.

It contributes to the PI output voltage via the proportional term :math:`K_\mathrm{p} = R_{35} / R_{39}`, and via its time integral, scaled by the integral time constant :math:`T_\mathrm{int} = C_{32} R_{39}`. (In these relations, the nominal ~0.2% difference between :math:`R_{39}` and :math:`R_{37} + R_{38}` is neglected.) The correct choice of proportional gain and integral time constant is critical for the stability of the regulated system. The two values should be chosen to match the system to be controlled. We have found by experimenting that the factory default for the integral time constant of the PI regulator (:math:`T_\mathrm{int} = 10\,\mathrm{ms})` is not ideal. Reducing :math:`C_{32}` --- and the integral time --- by roughly a factor of two (:math:`C_{32}=470\,\mathrm{nF}`) gives a stability improvement. Effectively this makes the feedback faster, and increases the weight of the integral component in :math:`V_\mathrm{PI}`.

The output of the PI regulator controls the DC--DC step-down converter that generates the supply voltage for the HV transformer, as shown in the bottom-right part of :numref:`fig_ngleed_beam_board_response_time`. The step-down supply is configured as a `buck converter <https://en.wikipedia.org/wiki/Buck_converter>`_. The output of the `LT1074 switching regulator <https://www.analog.com/media/en/technical-documentation/data-sheets/1074fds.pdf>`_ (VSW) is low-pass LC-filtered through L1 and (C21 + C23), generating the supply voltage for the transformer. The same DC voltage is fed back to the FB pin (via R24 and C40), where it adds to the output of the PI regulator. This way, an error in the HV output of the board, which causes a :math:`V_\mathrm{PI}` voltage, unbalances the feedback that generates the DC voltage. The LT1074 will then changed its DC output to restore the balance in the feedback, causing a regulation of the HV output of the board.

The reactivity of the feedback around the LT1074 step-down converter depends on the time constants of the (L1, C21 + C23) and (R23 || R24, C40) networks. By experimenting, we have found that reducing these time constants -- by reducing the C23 and C40 capacitors --- improves the time response of the beam energy. Good values for the capacitances are C40 = 470 nF and C23 = 470 μF. It is important that the electrolytic capacitor that will replace C23 can sustain large ripple currents at 100 kHz (the frequency of the LT1074), as it is responsible for filtering out the ripple. Before removing C23, make sure you clearly identify the polarization direction! It should be possible to replace all components (see :numref:`fig_ngleed_beam_board_response_time`) without removing the large heat sink where LT1074, D12, and the MOSFET transistors driving the transformer are mounted.

.. todo::
    - Try changing also R35 to modify the P gain. The best way would be: (i) short temporarily C32 and increase R35 till oscillations appear on I0. Pick R35 as half the value. Then try swapping around C32.
    - Figure out whether the problem is when the thing is **loaded** by disconnecting all boards (except ANODE that makes the 200V), then reconnect one at a time.
    - Measure response curve after all modifications, and add a figure. This way people can judge whether they feel the modifications are worth the effort.

.. _ngleed_suppressor_range:

Modifying the ``SUPPRESSOR`` range
==================================

.. admonition:: TL;DR

    The factory rage of the ``SUPPRESSOR`` dial is not ideal. This can be improved as follows:

      * Adjust potentiometer R24 as described in :secref:`ngleed_suppressor_feedback` and shown in :numref:`fig_ngleed_suppressor_feedback`.
      * Add a 100 kΩ resistor in parallel to R125, as shown in :numref:`fig_ngleed_suppressor`.


Every LEED setup uses (at least one) grid to repel inelastically scattered electrons. The grid is biased at a voltage that is normally slightly smaller (in absolute value) than the beam energy. For LEED-:math:`I(V)` measurements we have found that the optimal retarding voltage is between 80% and 95% of the beam energy. NG LEED units have a design fault for what concerns the range of retarding voltages that can be accessed. On our unit, the factory range was between 81% and 111% of the beam energy. These values are close to the ones specified by design (80--110%). Setting the suppressor voltage to energies larger than the beam energy makes no sense: there are no electrons going through anyway. This means that:

#. Approximately one third of the range of the one-turn potentiometer is useless;
#. At least part of the good range of suppressor voltages is outside the accessible range. On our unit, the 'best-looking' settings used to be very close to where the potentiometer switches on. (In the off position, the retarding voltage is set to 0% of the beam energy);
#. It is hard to control the actual voltage because of the large sensitivity of the setting to the potentiometer position.

We can fix all of this very simply, as described in :secref:`ngleed_suppressor_setpoint`. A few preliminary checks are described in :secref:`ngleed_suppressor_front_panel` and :secref:`ngleed_suppressor_feedback`.

The ``RETARD`` HV module works conceptually identically to all other HV modules in the NG LEED unit. As described in :secref:`ngleed_faster_beam` (especially :numref:`fig_ngleed_beam_board_flow`), the high-voltage is generated in a closed-loop configuration: a scaled version of the high-voltage output is subtracted from the 0--10 V setpoint and fed to a PI regulator, which generates the control signal for the portion of the circuit that drives the HV transformer. The primary difference between the ``BEAM`` and ``RETARD`` HV modules is in the range of accessible voltages. The ``BEAM`` HV module generates up to 1000 V in LEED mode and up to 3500 V in Auger mode. The ``RETARD`` HV module has a nominal maximum output of 2000 V, irrespective of the mode. The different maximum outputs are realized in a very simple manner. With reference to :numref:`fig_ngleed_beam_board_flow`, it is enough to change the gain of the feedback in order to modify the maximum output voltage. As already mentioned in :secref:`ngleed_faster_beam`, in LEED mode the DC feedback gain of the ``BEAM`` HV module is :math:`G_\mathrm{feedback}^\mathrm{BEAM, LEED}=1/100`, giving an output voltage

.. math:: V_\mathrm{HV,max}^\mathrm{BEAM, LEED} = -\frac{V_\mathrm{setpoint, max}}{G_\mathrm{feedback}^\mathrm{BEAM, LEED}} = -\frac{10\,\mathrm{V}}{0.01} = -1000\,\mathrm{V}.

In Auger mode, :math:`G_\mathrm{feedback}^\mathrm{BEAM, Auger}=1/350`, so that :math:`V_\mathrm{HV,max}^\mathrm{BEAM, Auger}=-3500\,\mathrm{V}`. In the ``RETARD`` HV module, the DC feedback gain is fixed to :math:`G_\mathrm{feedback}^\mathrm{RETARD}\approx 1/200`, giving :math:`V_\mathrm{HV,max}^\mathrm{RETARD}\approx-2000\,\mathrm{V}` at :math:`V_\mathrm{setpoint}=10\,\mathrm{V}`. See :secref:`ngleed_suppressor_feedback` for more details on this feedback gain.

.. _ngleed_suppressor_front_panel:

Checking the front-panel reading
--------------------------------

Before proceeding, make sure you have all the necessary tools:
    - two accurate digital multimeters,
    - two sets of pin grabbers,
    - a suitable anti-static electronics workbench,
    - a suitable set of Phillips, Pozi, and flat-head screwdrivers.

Before modifying the circuit that produces the suppressor voltage, it is worth cross-checking that the reading on the front panel is indeed correct. In fact, you will normally use the value displayed to set the suppressor voltage.

.. warning::
    Checking the correctness of the value displayed on the front panel requires to temporarily override safety features of the NG LEED unit and to measure high voltages. The ViPErLEED team takes no responsibility for the damage to people and/or equipment that may result from the application of inappropriate safety measures when dealing with high voltages.

It is important to note that, in LEED mode, the front-panel display shows **the difference between the retarding voltage and the beam energy** \[see the schematic circuit diagram in :numref:`fig_ngleed_suppressor_check_display`\ (c)\]. This means that, for example, a displayed value of −50 V at 500 V beam energy corresponds to 90% retarding.

.. _fig_ngleed_suppressor_check_display:
.. figure:: /_static/hardware/SpectaLEED_modification/suppressor_check_display.svg
    :align: center

    Identification of high- and low-voltage versions of the suppressor retarding voltage. (a) HV connector on the rear panel where the actual suppressor voltage should be measured at high voltages. (b) portion of the ``CONTROL`` board that produces the low-voltage version of the suppressor voltage as well as the value displayed on the front panel. The "Suppressor low-V" can be measured between the right lead of R38 and the bottom lead of R26. (c) Circuit diagram of the relevant part of the ``CONTROL`` board.

To ascertain the correctness of the value displayed, it is easiest to compare the high-voltage output of the ``RETARD`` HV module with the low-voltage version used to generate the input of the display. This is named "Suppressor low-V" in :numref:`fig_ngleed_suppressor_check_display`, and it is generated on the ``CONTROL`` board. See :numref:`fig_ngleed_i0_control_and_buffer`\ (b) for the location of the board.

The "Suppressor low-V" can be accessed by measuring the voltage at R38, relative to ground. :numref:`fig_ngleed_suppressor_check_display`\ (b) shows the location of the correct lead of R38 (right lead, see white circle) and an example of where the ground potential can be accessed (lower lead of R26). Expect voltages in the range of few hundred millivolts.

The high-voltage output of the suppressor board can be accessed at the high-voltage connector on the rear panel, visible in :numref:`fig_ngleed_suppressor_check_display`\ (a). Both pins on the bottommost row of the high-voltage outputs carry the suppressor voltage. It is referred to chassis ground, accessible at either the PE screw or as the second (counting from the left) of the low-voltage pins on the HV connector. The NG LEED unit will not output any high voltage unless the interlock pins on the HV connector are shorted.

We strongly advise using two multimeters, prepared to measure the low- and high-voltage versions of the suppressor voltage **before** plugging in the mains cord and turning on the unit. Then, only operate the unit from the front panel. Before turning the unit on, make sure that the beam-energy knob is turned to the minimum. You can use the beam-energy setting to control the exact suppressor voltage, as the two are proportional.

.. warning::
    Most multimeters are rated at max 1000 V with proper cabling. It is a good idea **not to exceed 500 V** to limit the risk of damaging the multimeter (and yourself).

As visible in :numref:`fig_ngleed_suppressor_check_display`\ (a), the ``RET_MON`` signal coming from the ``RETARD`` HV module is downscaled by roughly a factor of 5 using R25, R26, and R27. As ``RET_MON`` is roughly 1/200 of the high-voltage output (for more details see :secref:`ngleed_suppressor_feedback`), "Suppressor low-V" should be 500 mV when the high-voltage output is −500 V. The R27 potentiometer can be adjusted so that this is the case.

.. note::
    On our unit there was no need to modify the R27 potentiometer, so there is a good chance it has been correctly factory-adjusted on other units too.

.. _ngleed_suppressor_feedback:

Checking the feedback gain
--------------------------

.. warning::
    Checking the feedback gain of the ``RETARD`` HV module requires to temporarily override safety features of the NG LEED unit and to measure high voltages. The ViPErLEED team takes no responsibility for the damage to people and/or equipment that may result from the application of inappropriate safety measures when dealing with high voltages.

Before proceeding, make sure you have all the necessary tools:
    - two accurate digital multimeters,
    - two sets of pin grabbers,
    - a suitable anti-static electronics workbench,
    - a suitable set of Phillips, Pozi, and flat-head screwdrivers.

This check may be more interesting for those that use their SpectaLEED setup also as a retarding-field analyser for Auger spectroscopy. There, the retarding voltage determines the accurateness of the energy scale. In practice, as explained at the end of the introduction to :secref:`ngleed_suppressor_range`, checking the correctness of the feedback gain equals to ensuring that the high-voltage output range is the nominal 0 to −2000 V.

.. _fig_ngleed_suppressor_feedback:
.. figure:: /_static/hardware/SpectaLEED_modification/suppressor_feedback.svg
    :align: center

    \(a\) Suggested contacts on the ``CONTROL`` board on which to measure the feedback signal for the ``RETARD`` HV module: ``RETARD_MON`` on the lower lead of R25, ground on the lower lead of R26. (b) Potentiometer that is used to set the feedback gain to 1/200. (c) Circuit diagram of the portion of the ``RETARD`` HV module that is responsible for feeding back the high-voltage output to the input of the PI regulator. An inverted version of the same signal goes to the ``CONTROL`` board as ``RETARD_MON``.

:numref:`fig_ngleed_suppressor_feedback` can be used as a guide to check the correctness of the feedback gain. The diagram in :numref:`fig_ngleed_suppressor_feedback`\ (c) shows the feedback portion of the circuit of the ``RETARD`` HV module: The ``RETARD_HV`` high-voltage output is low-pass filtered (cutoff ~340 Hz) and scaled by a factor :math:`\approx-4.3\times10^{-3}` with the first inverting feedback stage around op-amp U9. The second inverting stage around U6 is used to adjust the overall gain to :math:`5\times10^{-3} = 1/200`. This is accomplished by adjusting potentiometer R24. The location of this potentiometer on the ``RETARD`` HV module is shown in :numref:`fig_ngleed_suppressor_feedback`\ (b). Notice that, should adjustments be needed, they can be performed without removing the ``RETARD`` HV module.

The output of stage U6 is added (with an essentially one-to-one ratio) to the ``RETARD_SET`` setpoint voltage, and fed to the PI regulator (:math:`K_\mathrm{P} = 1`, :math:`T_\mathrm{int} = 22\,\mathrm{ms}`) that generates the control signal for the transformer and HV cascade, which in turn produce the high-voltage output (see also :numref:`fig_ngleed_beam_board_flow` for the block diagram of HV modules).

The output of U6 is also accessible on the ``CONTROL`` board, as stage U7 produces an inverted version that is used for the front-panel reading, labelled ``RETARD_MON`` in :numref:`fig_ngleed_suppressor_feedback`\ (c). As such, it is easy to verify the correctness of the feedback gain by measuring at the same time the ``RETARD_MON`` signal on the ``CONTROL`` board and the high-voltage output at the HV connector. For the high-voltage part, follow the indications in :secref:`ngleed_suppressor_front_panel`. ``RETARD_MON`` can be measured as indicated in :numref:`fig_ngleed_suppressor_feedback`\ (a): it is the voltage on the bottom lead of R25 relative to ground (bottom lead of R26). See :numref:`fig_ngleed_suppressor_check_display`\ (c) for the circuit diagram of the relevant ``CONTROL`` board section.

Since the correct gain is 1/200, expect ~2.5 V with −500 V on the output. There probably is not much point in changing the position of the R24 potentiometer if your gain is somewhat close to 1/200. On our unit, we deemed acceptable a feedback voltage of 2.482 V with −500 V on the output. This corresponds to a gain error of 0.7%, or to a maximum output voltage of −2014 V.

.. note::
    If you modify the gain of the feedback, keep in mind that, as described in :secref:`ngleed_suppressor_front_panel`, the same signal also produces the values shown on the front-panel display. It is then important to :ref:`re-adjust also the front-panel reading<ngleed_suppressor_front_panel>`.

.. _ngleed_suppressor_setpoint:

Modifying the suppressor setpoint
---------------------------------

.. note::
    The modifications described in this section will **not** affect the retarding voltages used in "Auger" mode. The setup can still be safely used as a retarding-field analyzer. No recalibration is necessary.

Before proceeding, make sure you have all the necessary components and tools:
    - one 100 kΩ resistor,
    - a suitable anti-static electronics workbench,
    - a suitable set of Phillips, Pozi, and flat-head screwdrivers,
    - a soldering iron and solder,
    - tweezers (not necessary, but useful).

.. _fig_ngleed_suppressor:
.. figure:: /_static/hardware/SpectaLEED_modification/suppressor.svg
    :align: center

    \(a\) Diagram of the portion of the ``CONTROL``-board circuit that generates the setpoint for the retarding voltage, and (b) location of R125 that should be modified with an additional 100 kΩ resistor in parallel.

:numref:`fig_ngleed_suppressor`\ (a) shows the portion of the circuit of the NG LEED unit that is responsible for producing the setpoint of the retarding voltage (found on the ``CONTROL`` board): the setpoint for the energy is divided using two fixed resistors (R125 and R126) and the potentiometer on the front panel (R124). The ``RETARD_SET`` voltage is

.. math::
    V_\mathrm{RETARD\_SET} = \frac{1 + \frac{R_{124}}{R_{125}}(1-x)}{1 + \frac{R_{124}+R_{126}}{R_{125}}} V_\mathrm{BEAM\_SET} = \frac{1 + \frac{10}{R_{125}}(1-x)}{1 + \frac{40}{R_{125}}} V_\mathrm{BEAM\_SET},

where :math:`x` is the fractional position of the front-panel potentiometer, and :math:`R_{125}` is in kilo-ohms in the last relation.

As mentioned in the introductory part of :secref:`ngleed_suppressor_range` and described in detail in :secref:`ngleed_suppressor_feedback`, the ``RETARD`` HV module has twice the gain of the ``BEAM`` module (in LEED mode). This means that :math:`V_\mathrm{RETARD\_SET} = 5\,\mathrm{V}` gives −1 kV suppressor, or 100% retarding voltage. With the factory value of R125 (27 kΩ), the ``SUPPRESSOR`` knob gives :math:`V_\mathrm{RETARD\_SET}` between 4.03 V and 5.52 V, corresponding to the 80.6--110.4% range of retarding voltages mentioned earlier.

Modifying the range, is then merely a matter of modifying the R125 resistor: 20 kΩ should give a maximum suppressor retardation of exactly 100% if the ratio of the gains of the ``RETARD`` and ``BEAM`` HV modules is exactly 2. However, replacing R125 is tricky as it requires removing the whole ``CONTROL`` board --- including the potentiometers on the front panel. It is significantly easier to add a 100 kΩ resistor in parallel to R125. This gives an equivalent resistance of :math:`27\mathrm{k}\parallel100\mathrm{k}=21.3\mathrm{k}`, and a theoretical 69.4--102% suppressor range. Using a 100 kΩ resistor has also the advantage of leaving a little room in case :ref:`the feedback gain is a little off<ngleed_suppressor_feedback>`. See :numref:`fig_bent_resistors` for suggestions on how to prepare the new resistor.

.. _ngleed_beamhv_port:

Adding a ``BEAM HV`` port
=========================

This section is only specific to the **older** versions of the NG LEED control units that do not already have an SHV connector labelled BEAM HV on the rear panel.

.. todo::
    For older units only, like the one in Omega

.. _ngleed_control_cable:

Preparing a cable for controlling the LEED energy
=================================================

Before proceeding, make sure you have all the necessary components and tools:

    - **Depending on your unit**: one 5-pole male DIN-connector plug with 45° contacts (e.g., RS 491-011) **or** one 3-pole male DIN-connector plug with 90° contacts (e.g., RS 786-3439). Take a look at the ``BEAM ENERGY PROGRAMMING INPUT`` socket, right under the ``INTERNAL/EXTERNAL`` switch;
    - one cable-mount male BNC connector (e.g., RS 112-1669 \[50 Ω\] or 112-1675 \[75 Ω\]); **TODO @Michael: which connector do we mount on the ViPErLEED box?**
    - one coaxial cable suitable for the BNC connector above. Typically RG58 (e.g., RS 176-2081 or 240-8087) for 50 Ω connectors, RG59 (e.g., RS 393-024) for 75 Ω connectors. Length: 1--2 m. The ViPErLEED interface box should be close the the NG LEED to reduce as much as possible the capacitive load on the unit;
    - soldering iron and solder;
    - insulation-stripping tools for coaxial cable.

.. _fig_ngleed_program_cable:
.. figure:: /_static/hardware/SpectaLEED_modification/beam_egy_cable.svg
    :align: center
    
    Schematic connections needed for the cable between the ViPErLEED interface box and the NG LEED unit in order to control the beam energy for LEED-:math:`I(V)`.

In order to control the LEED unit to perform LEED-:math:`I(V)` measurements, the ViPErLEED interface box generates a 0--10 V signal at the ``OUT 0--10 V`` BNC connector. **TODO: ref to section where we explain the box** This signal needs to be brought to the NG LEED unit that will nominally generate a 0--1000 V electron beam energy. You can use a simple coaxial cable assembled as suggested in :numref:`fig_ngleed_program_cable`. Notice that some NG LEED units may have a 5-pole DIN connector instead of the 3-pole one mentioned in the manual.
