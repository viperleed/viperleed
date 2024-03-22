"""Error codes for the rfactor and interpolation."""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2022-01-20'
__license__ = 'GPLv3+'

import logging
import warnings

error_codes = {
    # Negative values correspond to warnings and non-fatal errors
    # Warnings
    -1: "Undefined Warning",
    -703: "At least one beam has not enough datapoints for interpolation. Consider using a lower interpolation degree or dropping short beams.",
    -812: "At least one numerator or denominator in Pendry R-factor caclulation was NaN. This may be because at least one beams did not contain usable data",
    -903: "R factor grouping encountered denominator of 0 for at least one beam",
    -904: "R factor grouping for R2 encountered group with 0 total energy overlapp",
    # No error, normal exit
    0: "No error",
    # 1... General Errors
    1: "Undefined error",
    # Compiler test failure
    10: 'Compilation check failed - NaN values are not recognized. This may have unexpected consequences. Make sure compiler option "-fno-finite-math-only" is enabled.',
    # 2 ... prepare_beams
    # 21.. limit range
    # 211: "???", # removed
    212: "???",
    # 22 ... averaging
    220: "given n_beams_out > n_beams",  # deprecated, now unused
    221: "Common averaging interval for (at least) one beam is too short",
    223: "Averaging skipped but given n_beams_out != n_beams",
    # 23 ... smoothing
    30: "",
    40: "",
    # 6xx : Interpolation
    611: "Origin grid not sorted in ascending order",
    612: "Target grid not sorted in ascending order",
    621: "Error in pre-factorize LHS",
    631: "LAPACK error in solve_coefficients",
    # 7xx : other in prepare
    701: "Unknown R factor requested. Only Pendry(1) and R2(2) are allowed.",
    702: "Invalid data passed to V0r optimization. Check if correct R factor is selected.",
    703: "In prepare_beams: more derivatives requested than are available.",
    704: "2nd derivative is not implemented yet",  # 2nd derivative will likely be needed in the future
    # 8: R-factor
    811: "At least one R-factor returned as NaN",
    851: "V0r optimization range too small",
    860: "Lapack problem during parabola fit",
    852: "All V0r steps sampled during fast search - V0r range may be too small",
    853: "Not able to find a suitable point to evaluate R(V0r)",
    854: "Fit range moved outside of V0r range. Resorted to V0r brute force.",
    856: "Minimum found in V0r shift search not well behaved. Consider changeing tolerances.",
    # 9xx : other
    902: "R factor grouping was given invalid group indices",
    903: "",  # changed to warning -903
    904: "",  # changed to warning -904
    921: "Error in apply_beamset_shift, check value for fill_outside",
}


def check_ierr(ierr, logger=None):
    if ierr and ierr > 0:
        # Error
        err_msg = f"ViPErLEED Fortran error code {ierr}: {error_codes[ierr]}"

        if logger:
            logger.error(err_msg)
        else:
            raise RuntimeError(err_msg)

    elif ierr and ierr < 0:
        # Warning
        warn_msg = f"ViPErLEED Fortran warning code {ierr}: {error_codes[ierr]}"

        if logger:
            logger.warn(warn_msg)
        else:
            warnings.warn(warn_msg, category=RuntimeWarning)
