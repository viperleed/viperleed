! Created by Alexander M. Imre on 04.12.21.

module interpolation
    use, intrinsic :: iso_fortran_env
    implicit none
    
    integer, parameter :: sp = REAL32
    integer, parameter :: dp = REAL64
    integer, parameter :: qp = REAL128
    
    
    contains

    subroutine equidist_quint_spline_coeffs(in_data, n_in, coeffs)
        integer, intent (in) :: n_in
        real(dp), intent (in) :: in_data(n_in)

        real(dp), intent (out) :: coeffs(5*n_in)




    end subroutine equidist_quint_spline_coeffs

    subroutine interpolate(x, y, deg, bc_type, do_checks, ierr)
        
        integer, INTENT(IN) :: n ! length of x, y#
        integer, INTENT(IN) :: deg ! degree
        real(dp), INTENT(IN) :: x(n), y(n), knots
        integer, INTENT(IN) :: bc_type, do_checks
        integer, INTENT(OUT) :: ierr ! error code

        ! Checks: x sorted ascending, deg sensible
        ! Not checked: x same size as y
    end subroutine interpolate

    subroutine interpolte_knots(x, y, n, knots, bc_type, do_checks, ierr)

        ! called with known knots

        !TODO implement checks

        ierr = 10 ! undefined error

        !allocate(Ab) stuff
        
    end subroutine interpolte_knots

    ! bc_types: 0 = "not-a-knot", 1 = natural, 2 = clamped

    integer pure function get_n_knots(deg, n) result(n_knots)
        
        n_knots = n + deg + 1
        RETURN
    end function get_n_knots

pure function get_knots_bc_0(x, n, deg) result(knots)
        ! in
        integer, intent(in) :: n, deg
        real(dp), intent(in) :: x(n)
        ! out
        real(dp) :: knots(n + deg + 1)
        ! internal
        integer :: m

        m = (deg - 1) / 2 ! integer divison

        knots(m:m+n) = x
        knots(1:m-1) = x(1)
        knots(m+n+1:n+deg+1) = x(n)

    end function get_knots_bc_0





end module interpolation