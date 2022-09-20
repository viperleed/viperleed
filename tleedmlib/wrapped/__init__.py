import warnings

try:
    from .rfactor import r_factor_new as rf
    from .interpolation import interpolation as intpol
except RuntimeError:
    warnings.warn("F2PY compiled module rfactor and/or interpolation failed to import.")
    warnings.warn("Check that these modules were compiled with the correct version of numpy/F2PY.")
    raise