#!/bin/bash


f2py -m interpolation interpolation.f90 -h interpolation.pyf --overwrite-signature --debug-capi only: get_array_sizes de_Boor_size single_calc_spline pre_eval_input_grid calc_spline_with_pre_eval single_interpolate_coeffs_to_grid get_intervals calc_deBoor eval_bspline_fast
f2py -c --f90flags="-g -Og -llapack" interpolation.pyf interpolation.f90