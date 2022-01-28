"""
Error codes for the rfactor and interpolation

"""

error_codes = {
    # 1... General Errors
    1: "Undefined error",
    10: "",
    # 2 ... prepare_beams
    # 21.. limit range
    211: "???",
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
    631: "LAPACK error in solve_coefficients"
}
