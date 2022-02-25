"""
Error codes for the rfactor and interpolation

"""

error_codes = {
    # 1... General Errors
    1: "Undefined error",
    10: "",
    # 2 ... prepare_beams
    # 21.. limit range
    # 211: "???", # removed
    212: "???",
    # 22 ... averaging
    220: "given n_beams_out > n_beams",
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
    903: "R factor grouping encountered denominator of 0 for at least one beam",
    904: "R factor grouping for R2 encountered group with 0 total energy overlapp",

    921: "Error in apply_beamset_shift, check value for fill_outside",
}
