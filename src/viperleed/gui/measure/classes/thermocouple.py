"""Module thermocouple of viperleed.gui.measure.classes.

Defines the base _Thermocouple class as well as a bunch of thermocouple
classes for the coefficients found in thermocouple_coefficients.txt.
The file was downloaded on 2022-07-26 from the NIST web page at
https://srdata.nist.gov/its90/download/download.html. Should other
thermocouple types be needed, their coefficients can be added to the
file, and will be generated automatically. If this is the case, the
format should be kept strictly as in the NIST version. In case new
coefficients are added, make sure all symbols used are UTF-8. For
example, saving the NIST file for all coefficients gives non-UTF-8
symbols for the "degrees".
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2022-07-26'
__license__ = 'GPLv3+'

from math import exp
from pathlib import Path
import sys


__all__ = ['read_coefficients_file', 'THERMOCOUPLES']


_COEFFICIENTS = {}
THERMOCOUPLES = ()

# TODO: will not work with pyinstaller
_coeff_fname = (Path(__file__).parent.resolve()
                / "thermocouple_coefficients.txt")


def read_coefficients_file():
    """Return thermocouple coefficients from thermocouple_coefficients.txt.

    Calling this function also updates the list of known thermocouples

    Returns
    -------
    None.
    """
    # pylint: disable=global-statement
    # We need the global because this method is intended to actually
    # update attributes of the module itself based on the coeff. info
    global _COEFFICIENTS, THERMOCOUPLES
    module = sys.modules[__name__]

    for type_ in _COEFFICIENTS:
        delattr(module, f'{type_}_Thermocouple')

    _COEFFICIENTS = {}
    with open(_coeff_fname, 'r', encoding='utf-8') as fproxy:
        # Look for a line that reads "type: ..."
        for line in fproxy:
            if line.startswith("type:"):
                type_ = line.split("type:")[1].strip()
                _COEFFICIENTS[type_] = _read_coefficients_block(fproxy, type_)

    tcs = []
    for type_ in _COEFFICIENTS:
        tc_name = f'{type_}_Thermocouple'
        tcouple = Thermocouple(type_)
        tcs.append(tcouple)
        setattr(module, tc_name, tcouple)
        __all__.append(tc_name)

    THERMOCOUPLES = tuple(tcs)


def _read_coefficients_block(fproxy, tc_type):
    """Read the block for thermocouple tc_type from fproxy."""
    # Check temperature and voltage units (next two lines)
    if next(fproxy).split('units:')[1].strip().lower() != "°c":
        raise RuntimeError("Unexpected temperature units for "
                           f"thermocouple {tc_type}. Expected '°C'")
    if next(fproxy).split('units:')[1].strip().lower() != "mv":
        raise RuntimeError("Unexpected voltage units for "
                           f"thermocouple {tc_type}. Expected 'mV'")

    # Now a series of blocks for the direct coefficients (°C->mV).
    direct, special = _read_direct_coefficients_block(fproxy)

    # Finally, read the inverse ones (mV->°C)
    inverse = _read_inverse_coefficients_block(fproxy, tc_type)

    if not direct:
        raise RuntimeError(f"No direct coefficients for {tc_type}-TC")
    if not inverse:
        raise RuntimeError(f"No inverse coefficients for {tc_type}-TC")

    if tc_type == 'K':
        if not special:
            raise RuntimeError("No exponential coefficients for K-TC")
        return _KThermocoupleCoefficients(direct, inverse, special)
    return _ThermocoupleCoefficients(direct, inverse)


def _read_direct_coefficients_block(fproxy):
    """Return the direct coefficients read from fproxy."""
    # Each block begins with "range: min, max, no. coefficients - 1"
    direct_coefficients = []
    exp_coefficients = []
    for line in fproxy:
        if not line.strip():
            break
        if line.startswith('range:'):
            min_, max_, n_coeff = line[6:].split(', ')
            n_coeff = int(n_coeff) + 1
            coefficients = [float(next(fproxy)) for _ in range(n_coeff)]
            direct_coefficients.append(
                (float(min_), float(max_), tuple(coefficients))
                )
            continue
        if line.startswith('exponential:'):
            # Special case, currently known only for K thermocouple
            for _line in fproxy:
                if not _line.split():
                    break
                exp_coefficients.append(float(_line.split('=')[1]))
    return direct_coefficients, exp_coefficients


def _read_inverse_coefficients_block(fproxy, tc_type):
    """Return the inverse coefficients for tc_type read from fproxy."""
    # The block of inverse coefficients starts with a line reading:
    begin_inverse = f"Inverse coefficients for type {tc_type}:"
    error = f"Could not find inverse coefficients for thermocouple {tc_type}"
    for line in fproxy:
        if line.startswith(begin_inverse):
            break
        if line.startswith("type:"):
            # This is the next thermocouple
            raise RuntimeError(error)
    else:
        raise RuntimeError(error)

    # Now go on till we reach a line that says
    # "  Voltage   <N numbers>\n", followed by
    # "  Range:    <N numbers>\n".
    # These are the min and max of the ranges
    for line in fproxy:
        if '*' in line:
            # next block!
            raise RuntimeError(error)
        if "Voltage" in line:
            break
    else:
        raise RuntimeError(error)
    minima = [float(m) for m in line.strip().split("Voltage")[1].split()]
    maxima = next(fproxy)
    if "Range:" not in maxima:
        raise RuntimeError(error)
    maxima = [float(m) for m in maxima.strip().split("Range:")[1].split()]
    assert len(minima) == len(maxima)

    # Finally, the coefficients
    coefficients = []
    for line in fproxy:
        if not line.split():
            continue
        if "Error" in line or "*" in line:
            # Read all
            break
        coefficients.append([float(c) for c in line.strip().split()])

    # Make sure all coefficients are the same length
    assert all(len(c) == len(minima) for c in coefficients)

    # Notice the zip* to swap columns and rows
    return [(min_, max_, tuple(c))
            for min_, max_, c in zip(minima, maxima, zip(*coefficients))]


class _ThermocoupleCoefficients:

    def __init__(self, temp_to_volt, volt_to_temp, *_):
        """Initialize instance."""
        self.__temp_to_volt = temp_to_volt
        self.__volt_to_temp = volt_to_temp

    @property
    def temperature_range(self):
        """Return the min and max temperature (°C) representable."""
        minima, maxima = zip(*((mi, ma) for mi, ma, _ in self.__temp_to_volt))
        return min(minima), max(maxima)

    @property
    def voltage_range(self):
        """Return the min and max voltages (mV)."""
        minima, maxima = zip(*((mi, ma) for mi, ma, _ in self.__volt_to_temp))
        return min(minima), max(maxima)

    def temp_to_volt(self, temperature):
        """Return the correct coefficients for this temperature."""
        for t_min, t_max, coefficients in self.__temp_to_volt:
            if t_min <= temperature <= t_max:
                return coefficients
        raise ValueError("No coefficients for this temperature.")

    def volt_to_temp(self, voltage):
        """Return the correct coefficients for this temperature."""
        for v_min, v_max, coefficients in self.__volt_to_temp:
            if v_min <= voltage <= v_max:
                return coefficients
        raise ValueError("No coefficients for this voltage.")


class _KThermocoupleCoefficients(_ThermocoupleCoefficients):

    def __init__(self, temp_to_volt, volt_to_temp, exp_coefficients, *_):
        """Initialize instance."""
        super().__init__(temp_to_volt, volt_to_temp)
        self.__exp_coefficients = exp_coefficients

    @property
    def exp_coefficients(self):
        """Return coefficients for the exponential part of °C->mV."""
        return self.__exp_coefficients


class Thermocouple:
    """Base class representing a thermocouple."""

    def __init__(self, type_):
        """Initialize instance."""
        if type_ not in _COEFFICIENTS:
            raise ValueError(
                f"Unkown thermocouple type {type_}. Make sure its "
                "coefficients are present in thermocouple_coefficients.txt"
                )
        self.__type = type_
        self.__coeff = _COEFFICIENTS[type_]
        self.__temperature_range = self.__coeff.temperature_range
        self.__voltage_range = self.__coeff.voltage_range

    def __str__(self):
        """Return a string representation of self."""
        return f"{self.type_}-type thermocouple"

    def __repr__(self):
        """Return a string representation of self."""
        return f"{self.__class__.__name__}({self.type_})"

    @property
    def coefficients(self):
        """Return the _ThermocoupleCoefficients of this thermocouple."""
        return self.__coeff

    @property
    def temperature_range(self):
        """Return the minimum and maximum temperatures (°C)."""
        return self.__temperature_range

    @property
    def type_(self):
        """Return the thermocouple type as a string."""
        return self.__type

    def mvolt_at(self, temperature):
        """Return the thermocouple voltage (mV) at a given temperature (°C)."""
        min_t, max_t = self.temperature_range
        if temperature < min_t or temperature > max_t:
            raise ValueError(f"Temperature {temperature} outside "
                             f"acceptable range ({min_t}, {max_t})")
        coefficients = self.coefficients.temp_to_volt(temperature)
        mvolt = sum(ck * temperature**k for k, ck in enumerate(coefficients))

        if self.type_ == 'K' and temperature > 0:
            # K thermocouple is the only one with a different equation
            coeff = self.coefficients.exp_coefficients
            mvolt += coeff[0] * exp(coeff[1] * (temperature - coeff[2])**2)
        return mvolt

    def temperature(self, voltage, cj_temperature=None):
        """Return the temperature for a given voltage.

        Parameters
        ----------
        voltage : float
            Thermocouple voltage in millivolts
        cj_temperature : float, optional
            Cold-junction compensation temperature in degrees
            centigrade. If not given the value returned is in
            fact a temperature DIFFERENCE.

        Returns
        -------
        temperature : float
            Temperature in degrees centigrade

        Raises
        ------
        ValueError
            If voltage, after applying the cold-junction compensation
            if given, lies outside the range of acceptable voltages
            for the known coefficients.
        """
        min_v, max_v = self.__voltage_range
        if cj_temperature:
            voltage += self.mvolt_at(cj_temperature)

        if voltage < min_v or voltage > max_v:
            raise ValueError(
                f"(Corrected) thermocouple voltage {voltage} outside "
                f"range of known coefficients ({min_v}, {max_v})"
                )

        coefficients = self.coefficients.volt_to_temp(voltage)
        return sum(ck * voltage**k for k, ck in enumerate(coefficients))


if not _COEFFICIENTS:
    # Update list of known thermocouples. Only the first time
    # the module is imported. If necessary, another explicit
    # call will update the information.
    read_coefficients_file()
