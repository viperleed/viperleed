! This module contains subroutines that are useful for the calculations of
! R-factors, and are based on the original TensErLEED code (v1.6). This module
! is compiled to a python module by using f2py with the following
! 
! WINDOWS:
!    f2py -c rfactor.pyf rfactor.f90  --fcompiler=gnu95 --compiler=mingw32      << REMEMBER TO CREATE .pyf WITH f2py -m rfactor rfactor.f90 -h rfactor.pyf BEFORE
!
! UNIX:
!
! MACOS:
!
! test with gfortran rfactor.f90 -c
!
! TODO: consider using optimization options for compilation with --opt=<string> flag
!       and libraries (lapack/blas) with -l<libname>
!
! One can alternatively specify other compiler options depending on the system,
! --fcompiler is the Fortran compiler, --compiler is the C compiler that is used
! by f2py to create a wrapper around the Fortran functions
!
! Author: Michele Riva
! Created: 2021-02-28

subroutine derivative(data_in, dx, deriv, n_data)
!   subroutine derivative computes the first derivative of the data array given
!   as input, using a step size dx of the independent variable. The routine uses
!   the following numerical approximants for the derivative, which depend on the
!   number of data points in the original array (as per TensErLEED code):
!   - if more than 22 points, use 3rd order numerical derivatives, i.e.:
!       3 points to left and right in the "inside" of the array, 7 pts in total
!       3 points to the right (left) at the left (right) edge, 4 points in total
!   - if between 2 and 22 points, use 2nd order numerical derivatives, i.e.:
!       1 point to left and right in the "inside" of the array, 3 pts in total
!	    2 points to the right (left) at the left (right) edge, 3 pts in total
!   - if only two points, use 1st-order forward/backward at edges
!
!   See also: https://en.wikipedia.org/wiki/Finite_difference_coefficient
!   for the coefficients

    implicit none
    integer, intent(in) :: n_data
    real(8), intent(in), dimension(0:n_data-1) :: data_in
    real(8), intent(in) :: dx
!f2py real(8), optional,intent(in) :: dx=1.
    real(8), intent(inout), dimension(0:n_data-1) :: deriv
!f2py real(8), optional, intent(in,out), dimension(0:n_data-1) :: deriv
    real(8) :: scaling
    integer :: n7, n6, n5, n4, n3, n2, n1
    
    n7 = n_data - 7
    n6 = n_data - 6
    n5 = n_data - 5
    n4 = n_data - 4
    n3 = n_data - 3
    n2 = n_data - 2
    n1 = n_data - 1
    
    deriv = 0.
    if (n_data == 1) then
        goto 10            ! TODO: probably replace with some complaint
    else if (n_data > 22) then
!       do first the interior of the array
        scaling = 1./(60*dx)
        deriv(3:n2) = (data_in(6:n1) - 9.*data_in(5:n2) + 45.*data_in(4:n3) &
                       - 45.*data_in(2:n5) + 9.*data_in(1:n6) - data_in(0:n7))*scaling
!       then the edges
        scaling = 1./(6*dx)
        deriv(0:2) = (2.*data_in(3:5) - 9.*data_in(2:4) &
                      + 18.*data_in(1:3) - 11.*data_in(0:2))*scaling
        deriv(n3:n1) = (11.*data_in(n3:n1) - 18.*data_in(n4:n2) &
                        + 9.*data_in(n5:n3) - 2.*data_in(n6:n4))*scaling
    else if (n_data > 2) then
!       do first the interior of the array
        scaling = 0.5/dx
        deriv(1:n2) = (data_in(2:n1) - data_in(0:n3))*scaling
!       then the two elements at the edges
        deriv(0) = (-3.*data_in(0) + 4.*data_in(1) - data_in(2))*scaling
        deriv(n1) = (3.*data_in(n1) - 4.*data_in(n2) + data_in(n3))*scaling
    else
        deriv(0) = (data_in(1) - data_in(0))/dx
        deriv(1) = deriv(0)
    end if

!   to be seen if it makes more sense to have derivative_out passed from outside
!   or not. This may save a bit of computation time.

10  continue
end subroutine derivative


subroutine pendry_y(intensity, de, v0i, i_prime, y_func, n_data)
! Calculate the Pendry Y function starting from an array of intensities and an
! energy step

    implicit none
    integer, intent(in) :: n_data
    real(8), intent(in) :: de
!f2py real(8), optional,intent(in) :: de=1.
    real(8), intent(in), dimension(0:n_data-1) :: intensity
    real(8), intent(in) :: v0i                 ! imaginary part of inner potential
!f2py real(8), optional,intent(in) :: v0i=5.
    real(8), dimension(0:n_data-1) :: i_prime  ! derivative of intensity
!f2py intent(hide) :: i_prime                  ! invisible from python interface
    real(8), dimension(0:n_data-1) :: y_func   ! where the y function is stored
!f2py real(8), optional, intent(in, out), dimension(0:n_data-1) :: y_func      ! leave the caller the option to provide buffer
    external :: derivative

    call derivative(intensity, de, i_prime, n_data)
    y_func = intensity * i_prime / (intensity**2 + v0i**2 * i_prime**2)  ! TODO: Maybe can be better optimized?

end subroutine pendry_y