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

        !in
        integer, INTENT(IN) :: deg,n
        
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

    subroutine build_colloc_matrix(x, n, knots, n_knots, deg, AB, ab_cols, ab_rows)
        ! -> bspl._colloc from scipy
        use, intrinsic :: iso_fortran_env
        implicit none
        integer, parameter :: dp = REAL64

        ! in
        integer, INTENT(IN) :: deg, n, n_knots, ab_cols, ab_rows
        real(dp), intent(IN) :: x(n), knots(n_knots)
        ! out
        real(dp), INTENT(INOUT) :: AB(ab_cols,ab_rows)
        ! internal
        integer :: left, j, a, kl, ku, clmn, nt
        real(dp) :: x_val
        real(dp), ALLOCATABLE :: work(:)

        ku = k
        kl = k
        nt = size(knots) - deg - 1


        allocate(work(2*k+2))

        left = deg
        ! TODO fix indices
        do j=0, n
            x_val = x(j+1)
            left = find_interval(t,k, x_val, left, extrapolate)

            deBoor_D(t, x_val, deg, left, 0, work)

            do a = 0, k+1
                clmn = left - k + a
                AB(kl + ku + j + offset - clmn + 1, clmn + 1) = work(a +1)
            end do
        end do

    end subroutine build_colloc_matrix

    subroutine deBoor_D(knots, x, deg, ell, m, result)
        ! analogous to _deBoor_D from fitpack.h in scipy
        ! Comment from SciPy:

    
        ! On completion the result array stores
        ! the k+1 non-zero values of beta^(m)_i,k(x):  for i=ell, ell-1, ell-2, ell-k.
        ! Where t[ell] <= x < t[ell+1].
        !
        ! Implements a recursive algorithm similar to the original algorithm of
        ! deBoor.

        !in
        integer, INTENT(IN) :: deg !degree
        integer, INTENT(IN) :: m ! which derivatives to evaluate

        real(dp), INTENT(IN) :: knots

        real(dp), INTENT(INOUT) :: result(deg+1)

        !internal
        real(dp) :: xb, xa, w
        integer :: ind, j, n
        real(dp) :: hh(deg)

        ! Comment from SciPy:
        ! Perform k-m "standard" deBoor iterations
        ! so that h contains the k+1 non-zero values of beta_{ell,k-m}(x)
        ! needed to calculate the remaining derivatives.
        
        result(1) = 1.0d
        j =1
        do while (j<= k-m)
            do i=1,j
                hh(i) = result(i)
            end do
            j = j+1
            result(1) = 0.0d
            n = 1
            do while (n<=j)
                n = n+1
                ind = ell +n
                xb = t(ind +1)
                xa = t(ind -j +1)
                if (xb == xa) then
                    h(n+1) = 0.0d
                    continue
                end if
                w = hh(n-1 +1)/(xb-xa)
                result(n-1+1) += result(n-1+1) + w*(xb-x)
                result(n+1) = w*(x-xa)
            end do
        end do

        ! Comment form SciPy:
        ! Now do m "derivative" recursions
        ! to convert the values of beta into the mth derivative
        j = k-m+1
        do while(j<=k)
            do i=1,j
                hh(i) = result(i)
            end do
            j = j+1
            h(1) = 0.0d
            n = 1
            do while(n<=j)
                n= n+1
                ind = ell + n
                xb= knots(ind +1)
                xa = knots(ind -j +1)
                if (xb == xa) then
                    h(m +1) = 0.0d
                    continue
                end if
                w = j*hh(n-1+1)/(xb-xa)
                h(n-1+1) = h(n-1) -w
                h(n+1) = w
            end do
        end do

    end subroutine deBoor_D

    integer pure function find_interval(knots, n_knots, deg, x_val, prev_l, extrapolate) result(interval)
        use, intrinsic :: iso_fortran_env
        integer, parameter :: dp = REAL64

        integer, INTENT(IN) :: n_knots, deg, extrapolate, prev_l
        real(dp), INTENT(IN) :: knots(n_knots), x_val

        ! internal
        integer :: l, n
        real(dp) :: tb, te

        ! TODO Fortran indices fix!
        n = n_knots - deg - 1
        tb = knots(deg+1)   !Fortran index
        te = knots(n+1)     !Fortran index

        ! reference code checks for NaN in x_val -> skipped here

        if (((x_val < tb).or.(x_val > te)).and.(extrapolate == 0)) then
            interval = -1
            RETURN
        end if 

        if ((k < pref_l).and.(prev_l<n)) then
            l = prev_l
        else
            l = k
        end if

        do while ((x_val < knots(l+1)).and.(l.ne.k)) !Fortran index
            l = l-1
        end do

        l = l+1

        do while ((x_val>=knots(l+1)).and.(l.ne.n)) !Fortran index
            l = l + 1
        end do
        interval = l-1

        return
    end function find_interval

