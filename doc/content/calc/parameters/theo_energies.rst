.. include:: /substitutions.rst

.. _theo_energies:

THEO_ENERGIES
=============

THEO_ENERGIES can be defined to set the electron energy range on which the
reference calculation, as well as the search algorithm will be performed.

**Default**: if no experimental data file is available,
``THEO_ENERGIES = 20 800 3``, otherwise ``THEO_ENERGIES = Efrom Eto 3``,
with ``Efrom = min(exp. energy) - 3`` and ``Eto = max(exp. energy) + 3``

**Syntax:**

::

   THEO_ENERGIES = 20 600 3
   THEO_ENERGIES = 50
   THEO_ENERGIES = _ _ 2   ! get Efrom and Eto from experimental data, but use step 2 eV

**Acceptable values:** Three positive floats for explicit definition, OR one
float to calculate only one energy, OR use underscores to leave some values
but not others at default (e.g. to be defined by experimental beams).

Notice that, generally, it is a good idea to choose the start and end energies
such that the theoretical calculation interval is a few eV larger than the
experimental energy range. During the determination of the |R factor|\ s, as
well as in the search of the best structure, the theoretical spectra will be
shifted a bit along the energy axis to find the best match with experiment
(see :ref:`IVSHIFTRANGE` to control this shift).

If ``Eto - Efrom`` is not an integer multiple of ``Estep``, ``Efrom`` will
be modified to a slightly lower energy such that this is the case.

.. todo::
    @michele-riva said:
    Perhaps we should have a different default in case EXPBEAMS is there:
    Look at the IV_SHIFT_RANGE to decide how much to 'expand' the THEO_ENERGIES
    rather than having a fixed 3 eV expansion on both sides?
    @scmhid-iap said:
    Should probably be slightly more than IV_SHIFT_RANGE since interpolation
    is rather bad at the ends.

    @michele-riva said:
    The above is still not done as of June 2024. The complication with using
    IV_SHIFT_RANGE too is that the latter may have its bounds modified in
    iorfactor.py, as TensErLEED requires the bounds to be integer multiples
    of the step, which may be undefined till we reach an R-factor calculation.
    We could be a bit generous, as suggested by @scmhid-iap, and expand the
    energy range by the IV_SHIFT_RANGE limits + XXX eV on both sides. We have
    to choose XXX. Probably something >= 3eV. Perhaps XXX could be taken from
    the step of IV_SHIFT_RANGE, if available.
