program intpol_test
    
    use, intrinsic :: iso_fortran_env
    integer, parameter :: sp = REAL32
    integer, parameter :: dp = REAL64
    integer, parameter :: qp = REAL128
    use interpolation

    implicit none

    real(dp) :: data(11), data_out(30)
    data = (0,1,2,3,-3,-1,0,1,-2,3,0)


end program intpol_test

