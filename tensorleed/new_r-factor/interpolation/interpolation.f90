! Created by Alexander M. Imre on 04.12.21.

module interpolation
    use, intrinsic :: iso_fortran_env
    implicit none
    
    integer, parameter :: sp = REAL32
    integer, parameter :: dp = REAL64
    integer, parameter :: qp = REAL128
    
    
    contains

    subroutine cube_spline()
        ! input: x,y data, new grid
        ! output: y data on new grid
    end subroutine cube_spline

    subroutine equidist_spline()
        ! input: x,y data, new grid
        ! output: y data on new grid
    end subroutine equidist_spline

    subroutine equidist_quint_spline_coeffs()
        ! input: x,y data, new grid
        ! output: coefficiencts and knots for evaluation

        integer, intent (in) :: n_in
        real(dp), intent (in) :: in_data(n_in)

        real(dp), intent (out) :: coeffs(5*n_in)

    end subroutine equidist_quint_spline_coeffs



    subroutine interpolate(x, y, n, knots, n_knots, deg, bc_type, do_checks, ierr)
        
        integer, INTENT(IN) :: n, n_knots ! length of x, y, knots
        integer, INTENT(IN) :: deg ! degree
        real(dp), INTENT(IN) :: x(n), y(n), knots(n_knots)
        integer, INTENT(IN) :: bc_type, do_checks
        integer, INTENT(OUT) :: ierr ! error code

        ! Checks: x sorted ascending, deg sensible
        ! Not checked: x same size as y

        
    end subroutine interpolate

    subroutine interpolate_knots(x, y, n, knots, n_knots, deg, bc_type, do_checks, ierr)
        integer, INTENT(IN) :: n, n_knots ! length of x, y, knots
        integer, INTENT(IN) :: deg ! degree
        real(dp), INTENT(IN) :: x(n), y(n), knots(n_knots)
        integer, INTENT(IN) :: bc_type, do_checks
        integer, INTENT(OUT) :: ierr ! error code
        ! called with known knots

        !TODO implement checks

        ierr = 10 ! undefined error

        !allocate(Ab) stuff
        
    end subroutine interpolate_knots

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

    subroutine build_colloc_matrix(x, knots, deg, AB)
        ! -> bspl._colloc from scipy
    !kl = ku =k
        nt = size(knots) - deg - 1


        allocate(work(2*k+2))

        left = deg
        ! TODO fix indices
        do j=0, n:
            left = find_interval(t,k, x_val, left, extrapolate)

            deBoor_D(t, x_val, deg, left, 0, work)

            do a = 0, k+1:
                clmn = left - k + a
                AB(kl + ku + j + offset - clmn, clmn) = work(a)

    end subroutine build_colloc_matrix

    pure subroutine deBoor_D()
        !analog to _deBoor_D from fitpack.h in scipy

    end subroutine deBoor_D

    pure subroutine find_interval(knots, n_knots, deg, x_val, prev_l, extrapolate)
        integer, INTENT(IN) :: n_knots, deg, extrapolate
        real(dp), INTENT(IN) :: knots(n_knots)

        ! internal
        integer :: l, n
        real(dp) :: tb, te

        ! TODO Fortran indices fix!
        n = n_knots - deg - 1
        tb = knots(deg)
        te = knots(n)

        ! reference code checks for NaN in x_val -> skipped here

        if (((x_val < tb).or.(x_val > te)).and.(extrapolate == 0)):
            ! error?
        end if 

        if k < pref_l then
            l = prev_l
        else
            l = k
        end if

    end subroutine find_interval

