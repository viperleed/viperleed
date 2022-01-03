program intpol_test
    
    use, intrinsic :: iso_fortran_env

    use interpolation, only : equidist_quint_spline_coeffs
    implicit none
    
    integer, parameter :: sp = REAL32
    integer, parameter :: dp = REAL64
    integer, parameter :: qp = REAL128

    real(dp) :: data(11), data_out(30)
    data = (/0,1,2,3,-3,-1,0,1,-2,3,0/)

    print*, data

end program intpol_test

