! Created by Alexander M. Imre on 04.12.21.
!
! Fast Bspline Interpolation
! A part of ViPErLEED
!
! This module is designed to allow for a fast and repeated interpolation of values to a denser subgrid using natural boundary conditions.
! It allows to separate the setup of the system of equations form the actual interpolation.
! 
!

module interpolation
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic
    implicit none
    
    integer, parameter :: sp = REAL32
    integer, parameter :: dp = REAL64
    integer, parameter :: qp = REAL128
    
    contains

! inputs : x, y, n | deg, deriv | x_new, n_new
! outputs : y_new 
! setup LHS
    subroutine interp_equidist_single(n, x_start, x_step, y, n_new, x_new_start, x_new_step, do_checks, y_new, ierr)
        ! Used for ...
        ! Calculates x and x_new from given steps and calls interpolate_single
        real(dp), INTENT(IN) :: x_start, x_step
        real(dp), INTENT(IN) :: x_new_start, x_new_step

        real(dp) :: x(n), x_new(n_new)

        do i = 1, n
            x(i) = x_start + x_step*(i-1)
        end do
        do i = 1, n_new
            x_new(i) = x_new_start + x_new_step*(i-1)
        end do

        call interpolate_single(n, x, y, n, deg, n_new, x_new, do_checks, y_new, ierr)
        RETURN        
    end subroutine interp_equidist_single

    subroutine interpolate_single(n, x, y, n, deg, n_new, x_new, do_checks, y_new, ierr)
        integer, INTENT(IN) :: n ! length of x, y
        integer, INTENT(IN) :: deg ! degree
        real(dp), INTENT(IN) :: x(n), y(n)
        integer, INTENT(IN) :: do_checks
        
        integer, INTENT(IN) :: n_new ! number of new points
        real(dp), INTENT(IN) :: x_new(n_new) !x_new positions where to evaluate

        real(dp), INTENT(OUT) :: y_new(n_new) ! interpolated y values at x_new positions
        integer, INTENT(OUT) :: ierr ! error code

        ! internal variables and arrays
        integer :: n_knots
        real(dp), ALLOCATABLE :: knots(:)
        
        ierr = 0
        if (do_checks.ne.0) then
            call perform_checks()
        end if

        
        
    end subroutine interpolate_single

! *****************************************************************************************
! Pre-Evaluation
 
    subroutine equidist_pre_evaluate_grid()
    end subroutine equidist_pre_evaluate_grid


    subroutine pre_evaluate_grid(x, n, x_new, n_new, deg, nt, kl, ku, &
        knots, n_knots, LHS, LHS_rows, LHS_cols, RHS, RHS_cols, ipiv, intervals, deBoor_matrix)
        ! Sets up and pre-evaluate parts of the caluclation related to ...

        ! IN
        integer, INTENT(IN) :: deg              !! order of spline
        integer, INTENT(IN) :: n, n_new         !! # points in, out
        real(dp), INTENT(IN):: x(n)             !! origin grid values
        real(dp), INTENT(IN):: x_new(n_new)     !! result grid values


        ! OUT
        integer, INTENT(OUT)                :: n_knots              !! # knots
        real(dp), INTENT(OUT), ALLOCATABLE  :: knots(:)             !! knot points
        integer, INTENT(OUT)                :: nt, kl, ku           !! matrix shape specifiers
        integer, INTENT(OUT)                :: LHS_rows, LHS_cols   !! shape LHS
        integer, INTENT(OUT)                :: RHS_cols             !! shape RHS
        real(dp), INTENT(OUT), ALLOCATABLE  :: LHS(:,:)
        real(dp), INTENT(OUT), ALLOCATABLE  :: RHS(:)
        integer, INTENT(OUT), ALLOCATABLE   :: ipiv(:)
        integer, INTENT(OUT)                :: intervals(n_new)
        real(dp), INTENT(OUT), ALLOCATABLE  :: deBoor_matrix(:,:)


        ! Internal


        ! Compute knot vector and size
        call get_natural_knots(x, n, deg, knots, n_knots)

        ! Get the correct derivatives at the edge for the degree givennatural boundary conditions
        call prepare_derivs_nat_bc(deg, derivs_known_l, derivs_known_r, nleft, nright,&
        derivs_l_ord, derivs_r_ord, derivs_l_val, derivs_r_val)

        ! Using the derivatives and the x_values, we build up the LHS of the equation
        call build_LHS(x, n, knots, n_knots, deg, derivs_known_l, derivs_known_r, nleft, nright, &
        derivs_l_ord, derivs_r_ord, AB, ab_rows, ab_cols, nt, kl, ku, ierr)

        ! with this done, we can now prefactorize the LHS
        call pre_factorize_LHS(AB, ab_rows, ab_cols, nt, kl, ku, ipiv, info)

        if (info.ne.0) then
            ! Error in pre-factorization
            ierr = 25
        end if

        ! Now set up the RHS
        call prepare_RHS(nt, nleft, nright, derivs_l_val, derivs_r_val, rhs)

        ! Build up the deBoor Matrix and the intervals array
        call pre_eval_bspline(knots, n_knots, x_new, n_new, deg, 0, intervals, deBoor_matrix)

        RETURN
    end subroutine pre_evaluate_grid


    subroutine perform_checks(x, y, n, x_new, n_new, ierr)
        implicit none
        integer, intent(in) :: n, n_new
        real(dp), INTENT(IN) :: x(n), y(n), x_new(n)
        integer,intent(out) ::  ierr

        integer :: i ! loop var

        ierr = 0

        ! check if x values are sorted
        do i = 2, n
            if (x(i).le.x(i-1)) then
                ! origin grid not sorted
                ierr = 11
                RETURN
            end if
        end do

        do i = 2, n_new
            if (x_new(i).le.x_new(i-1)) then
                ! target grid not sorted
                ierr = 12
                RETURN
            end if
        end do

        RETURN
    end subroutine perform_checks

! *****************************************************************************************
! Fast Evaluation using prepared values


    subroutine interpolate_new(x, y, n, deg, do_checks, x_new, y_new, n_new, knots, n_knots, &
                                nleft, nright, LHS, RHS, LHS_rows, LHS_cols, RHS_cols, nt, kl, ku,&
                                ipiv, intervals, deBoor_matrix, ierr)
        integer, INTENT(IN) :: n, do_checks ! length of x, y, knots
        integer, INTENT(IN) :: deg ! degree
        real(dp), INTENT(IN) :: x(n), y(n)
        
        integer, INTENT(IN) :: n_new ! number of new points
        real(dp), INTENT(IN) :: x_new(n_new) !x_new positions where to evaluate

        integer, INTENT(OUT) :: n_knots
        real(dp), intent(OUT), ALLOCATABLE :: knots(:)
        real(dp), INTENT(OUT) :: y_new(n_new) ! interpolated y values at x_new positions
        integer, INTENT(OUT) :: ierr ! error code

        integer :: ab_cols, ab_rows
        integer, INTENT(out) :: kl, ku, nt
        real(dp), ALLOCATABLE :: AB(:,:)

        integer, INTENT(out) :: nleft, nright
        integer :: derivs_known_l, derivs_known_r
        integer, ALLOCATABLE :: derivs_l_ord(:), derivs_r_ord(:)
        real(dp), ALLOCATABLE :: derivs_l_val(:), derivs_r_val(:)

        integer :: LHS_rows, LHS_cols, RHS_cols
        real(dp), ALLOCATABLE :: rhs(:), coeffs(:), LHS(:,:)


        !for lapack
        integer :: info
        integer, INTENT(out), ALLOCATABLE :: ipiv(:)

        ! for eval bspline
        integer,  INTENT(out) :: intervals(n_new)
        real(dp), ALLOCATABLE, INTENT(out) :: deBoor_matrix(:,:)

        ierr = 0 ! undefined error

        !TODO implement checks
        if (do_checks.ne.0) then
            print*, "Checks not yet implemented"
            ierr = 99
        end if
        ! Compute knots, then hand off to interpolate_knots
        call get_natural_knots(x, n, deg, knots, n_knots)


        call prepare_derivs_nat_bc(deg, derivs_known_l, derivs_known_r, nleft, nright,&
        derivs_l_ord, derivs_r_ord, derivs_l_val, derivs_r_val)

        call build_LHS(x, n, knots, n_knots, deg, derivs_known_l, derivs_known_r, nleft, nright, &
                derivs_l_ord, derivs_r_ord, AB, ab_rows, ab_cols, nt, kl, ku, ierr)

        ! Set up right hand side
        ! extradim in scipy can be ignored since we are interpolating only one set of y values

        call prepare_RHS(nt, nleft, nright, derivs_l_val, derivs_r_val, rhs)
        call assign_y_to_RHS(y, rhs, nt, nleft, nright)

        !!!! Everything above does not need to be repeated if x and grid stay the same!!!


        call pre_factorize_LHS(AB, ab_rows, ab_cols, nt, kl, ku, ipiv, info)

        ! Now, here, we need to interject: RHS and LHS would be overwritten
        RHS_cols = nt
        LHS_rows = ab_rows
        LHS_cols = ab_cols

        ALLOCATE(LHS(LHS_rows, LHS_cols), coeffs(RHS_cols))

        LHS = AB
        coeffs = RHS


        call solve_coeffcients(nt, kl, ku, ab_rows, ab_cols, AB, coeffs, ipiv, info)

        if ((info >0).or.(info<0)) then
            print*, "LAPACK ERROR"
        end if


        ! Nice, if you made it here, you have all coefficients!
        
        ! Now use them to evaluate at x_new positions

        call pre_eval_bspline(knots, n_knots, x_new, n_new, deg, 0, intervals, deBoor_matrix)

        call evaluate_bspline(knots, n_knots, x_new, y_new, n_new, deg, 0, coeffs, nt)


        RETURN
    end subroutine interpolate_new

    subroutine interpolate_fast(deg, nt, kl, ku, knots, n_knots, LHS, RHS, LHS_rows, LHS_cols, RHS_cols, ipiv, x_new, n_new,&
                                intervals, deBoor_matrix, y_new)
        integer, INTENT(in) :: deg
        integer, INTENT(in) :: LHS_rows, LHS_cols, RHS_cols, n_new, n_knots
        integer, INTENT(INOUT) :: ipiv(nt)
        integer :: info

        real(dp), intent(in) :: LHS(LHS_rows, LHS_cols), RHS(RHS_cols), x_new(n_new), knots(n_knots)

        ! for eval bspline
        integer,  INTENT(in) :: intervals(n_new)
        real(dp), ALLOCATABLE, INTENT(in) :: deBoor_matrix(:,:)

        ! OUT
        real(dp), intent(out) :: y_new(n_new)
        ! INTERNAL
        real(dp) :: LHS_copy(LHS_rows, LHS_cols), coeffs(RHS_cols)
        integer  :: nt, kl, ku

        
        ! Copy LHS, RHS:
        LHS_copy = LHS
        coeffs = RHS

        call assign_y_to_RHS()

        call solve_coeffcients(nt, kl, ku, lhs_rows, lhs_cols, LHS_copy, coeffs, ipiv, info)
        if ((info >0).or.(info<0)) then
            print*, "LAPACK ERROR"
        end if
        
        ! Now use them to evaluate at x_new positions
        call eval_bspline_fast(deg, coeffs, intervals, deBoor_matrix, n_new, y_new)

        RETURN
    end subroutine interpolate_fast

! *****************************************************************************************
! Internal routines for setup

    subroutine prepare_derivs_nat_bc(deg, derivs_known_l, derivs_known_r, nleft, nright,&
        derivs_l_ord, derivs_r_ord, derivs_l_val, derivs_r_val)
        integer, INTENT(IN) :: deg
        

        integer, INTENT(OUT) :: derivs_known_l, derivs_known_r, nleft, nright
        integer, ALLOCATABLE, INTENT(OUT) :: derivs_l_ord(:), derivs_r_ord(:)
        real(dp), ALLOCATABLE, INTENT(OUT) :: derivs_l_val(:), derivs_r_val(:)

        if (deg==3) then
            derivs_known_l = 1
            derivs_known_r = 1
            ALLOCATE(derivs_l_ord(derivs_known_l), derivs_r_ord(derivs_known_r))
            ALLOCATE(derivs_l_val(derivs_known_l), derivs_r_val(derivs_known_r))

            derivs_l_ord = (/2/)
            derivs_l_val = (/0d0/)

            derivs_r_ord = (/2/)
            derivs_r_val = (/0d0/)

            
        else if (deg == 5) then
            derivs_known_l = 2
            derivs_known_r = 2
            ALLOCATE(derivs_l_ord(derivs_known_l), derivs_r_ord(derivs_known_r))
            ALLOCATE(derivs_l_val(derivs_known_l), derivs_r_val(derivs_known_r))
            derivs_l_ord = (/3, 4/)
            derivs_l_val = (/0d0, 0d0/)

            derivs_r_ord = (/3, 4/)
            derivs_r_val = (/0d0, 0d0/)
        else
            print*, "Only order 3 and 5 supported so far"
            stop
        end if

        nleft = size(derivs_l_ord)
        nright = size(derivs_r_ord)
        RETURN

    end subroutine prepare_derivs_nat_bc

    ! Called with known knots
    subroutine interpolate_knots(x, y, n, knots, n_knots, deg, x_new, y_new, n_new, ierr)
        integer, INTENT(IN) :: n, n_knots ! length of x, y, knots
        integer, INTENT(IN) :: deg ! degree
        real(dp), INTENT(IN) :: x(n), y(n), knots(n_knots)
        
        integer, INTENT(IN) :: n_new ! number of new points
        real(dp), INTENT(IN) :: x_new(n_new) !x_new positions where to evaluate

        real(dp), INTENT(OUT) :: y_new(n_new) ! interpolated y values at x_new positions
        integer, INTENT(OUT) :: ierr ! error code

        integer :: ab_cols, ab_rows, kl, ku, nt
        real(dp), ALLOCATABLE :: AB(:,:)

        integer :: derivs_known_l, derivs_known_r, nleft, nright
        integer, ALLOCATABLE :: derivs_l_ord(:), derivs_r_ord(:)
        real(dp), ALLOCATABLE :: derivs_l_val(:), derivs_r_val(:)

        real(dp), ALLOCATABLE :: rhs(:)

        !for lapack
        integer :: info
        integer, ALLOCATABLE :: ipiv(:)

        ierr = 0 ! undefined error

        call prepare_derivs_nat_bc(deg, derivs_known_l, derivs_known_r, nleft, nright,&
        derivs_l_ord, derivs_r_ord, derivs_l_val, derivs_r_val)

        ! Set up left hand side
        nt = n_knots - deg - 1
        kl = deg
        ku = deg
        ab_rows =2*kl + ku +1
        ab_cols = nt
        allocate(AB(ab_rows, ab_cols))
        AB = 0.0d0

        call build_colloc_matrix(x, n, knots, n_knots, deg, AB, ab_rows, ab_cols, nleft)
        if (nleft>0) then
            call handle_lhs_derivatives(knots, n_knots, deg, x(1), AB, ab_rows, ab_cols, kl, ku, &
            derivs_l_ord, derivs_known_l, 0)
        end if
        if (nright>0) then
            call handle_lhs_derivatives(knots, n_knots, deg, x(n), AB, ab_rows, ab_cols, kl, ku, &
            derivs_r_ord, derivs_known_r,  nt-nright)
        end if

        ! Set up right hand side
        ! extradim in scipy can be ignored since we are interpolating only one set of y values

        ALLOCATE(rhs(nt))
        if (nleft>0) then
            rhs(1:nleft) = derivs_l_val
        end if
        if (nright>0) then
            rhs(nt - nright +1: nt) = derivs_l_val
        end if


        !!!! Everything above does not need to be repeated if x and grid stay the same!!!
        rhs(nleft+1:nt-nright) = y ! assign y values

        ! Now we are ready to solve the matrix!
        ! Using LAPACKs' GBSV routine
        ! On exit, rhs is overwritten by the solution c
        ALLOCATE(ipiv(nt))
        call DGBSV(nt, kl, ku, 1, AB, 2*kl+ku+1, ipiv, rhs, nt, info)

        if ((info >0).or.(info<0)) then
            print*, "LAPACK ERROR"
        end if

        ! Nice, if you made it here, you have all coefficients!
        
        ! Now use them to evaluate at x_new positions

        call evaluate_bspline(knots, n_knots, x_new, y_new, n_new, deg, 0, rhs, nt)
        
        RETURN

    end subroutine interpolate_knots

    pure integer function get_n_knots(deg, n) result(n_knots)
        use, intrinsic :: iso_fortran_env
        implicit none
        integer, parameter :: dp = REAL64
        !in
        integer, INTENT(IN) :: deg,n
        
        n_knots = n + 2*deg
        RETURN
    end function get_n_knots

    subroutine get_natural_knots(x, n, deg, knots, n_knots)
        use, intrinsic :: iso_fortran_env
        implicit none
        integer, parameter :: dp = REAL64

        ! in
        integer,intent(in) :: n, deg
        real(dp), INTENT(IN) :: x(n)

        ! out
        real(dp), intent(out), ALLOCATABLE ::  knots(:)
        integer, INTENT(OUT) :: n_knots

        n_knots = get_n_knots(deg, n)
        ALLOCATE(knots(n_knots))
        knots(1       : deg     ) = x(1)
        knots(1+deg   : n+deg   ) = x(:)
        knots(n+deg+1 : n_knots ) = x(n)
        RETURN
    end subroutine get_natural_knots

    subroutine build_colloc_matrix(x, n, knots, n_knots, deg, AB, ab_rows, ab_cols, offset)
        ! -> bspl._colloc from scipy
        use, intrinsic :: iso_fortran_env
        implicit none
        integer, parameter :: dp = REAL64

        ! in
        integer, INTENT(IN) :: deg, n, n_knots, ab_cols, ab_rows
        real(dp), intent(IN) :: x(n), knots(n_knots)
        ! out
        real(dp), INTENT(INOUT) :: AB(ab_rows,ab_cols)
        ! internal
        integer :: left, j, a, kl, ku, clmn, nt, offset
        real(dp) :: x_val
        real(dp), ALLOCATABLE :: work(:)

        ku = deg
        kl = deg
        nt = n_knots - deg - 1

        allocate(work(2*deg+2))

        left = deg
        ! TODO fix indices
        do j=0, n-1
            x_val = x(j+1)
            left = find_interval(knots, n_knots, deg, x_val, left, 0) ! 0 is extrapolate = False

            call deBoor_D(knots, n_knots, x_val, deg, left, 0, work)

            do a = 0, deg
                clmn = left - deg + a
                AB(kl + ku + j + offset - clmn + 1, clmn +1 ) = work(a +1)
            end do
        end do

    end subroutine build_colloc_matrix

    subroutine deBoor_D(knots, n_knots, x, deg, ell, m, result)
        ! analogous to _deBoor_D from fitpack.h in scipy
        ! Comment from SciPy:

    
        ! On completion the result array stores
        ! the k+1 non-zero values of beta^(m)_i,k(x):  for i=ell, ell-1, ell-2, ell-k.
        ! Where t[ell] <= x < t[ell+1].
        !
        ! Implements a recursive algorithm similar to the original algorithm of
        ! deBoor.

        use, intrinsic :: iso_fortran_env
        implicit none
        integer, parameter :: dp = REAL64

        !in
        integer, INTENT(IN) :: n_knots, ell
        integer, INTENT(IN) :: deg !degree
        integer, INTENT(IN) :: m ! which derivatives to evaluate

        real(dp), INTENT(IN) :: x, knots(n_knots)

        real(dp), INTENT(INOUT) :: result(deg+1)

        !internal
        real(dp) :: xb, xa, w
        integer :: ind, j, n, i
        real(dp) :: hh(deg)

        ! Comment from SciPy:
        ! Perform k-m "standard" deBoor iterations
        ! so that h contains the k+1 non-zero values of beta_{ell,k-m}(x)
        ! needed to calculate the remaining derivatives.
        
        result(1) = 1.0d0
        j =1
        do while (j<= deg-m)
            do i=1,j
                hh(i) = result(i)
            end do
            result(1) = 0.0d0
            n = 1
            do while (n<=j)
                ind = ell +n
                xb = knots(ind +1)
                xa = knots(ind -j +1)
                if (xb == xa) then
                    result(n+1) = 0.0d0
                    continue
                end if
                w = hh(n-1 +1)/(xb-xa)
                result(n-1+1) = result(n-1+1) + w*(xb-x) ! result(n) += w*(xb-x)
                result(n+1) = w*(x-xa)
            n = n+1
            end do
        j = j+1
        end do

        ! Comment form SciPy:
        ! Now do m "derivative" recursions
        ! to convert the values of beta into the mth derivative
        j = deg-m+1
        do while(j<=deg)
            do i=1,j
                hh(i) = result(i)
            end do
            result(1) = 0.0d0
            n = 1
            do while(n<=j)
                ind = ell + n
                xb= knots(ind +1)
                xa = knots(ind -j +1)
                if (xb == xa) then
                    result(m +1) = 0.0d0
                    continue
                end if
                w = j*hh(n-1+1)/(xb-xa)
                result(n-1+1) = result(n-1+1) -w
                result(n+1) = w
            n= n+1
            end do
        j = j+1
        end do
    end subroutine deBoor_D

    subroutine handle_lhs_derivatives(knots, n_knots, deg, x_val, AB, ab_rows, ab_cols, kl, ku, deriv_ords, n_derivs, offset)
        use, intrinsic :: iso_fortran_env
        implicit none
        integer, parameter :: dp = REAL64

        ! in
        integer, INTENT(IN) :: deg, n_knots, ab_cols, ab_rows, offset, kl, ku, n_derivs, deriv_ords(n_derivs)
        real(dp), intent(IN) :: knots(n_knots), x_val
        ! out
        real(dp), INTENT(INOUT) :: AB(ab_rows,ab_cols)
        ! internal
        integer :: a, left, nu, row, clmn

        real(dp), ALLOCATABLE :: work(:)
        
        ! derivatives @ x_val
        left = find_interval(knots, n_knots, deg, x_val, deg, 0)

        ! allocate work array
        ALLOCATE(work(2*deg+2))
        do row=0, n_derivs-1
            nu = deriv_ords(row+1)
            call deBoor_D(knots, n_knots, x_val, deg, left, nu, work)

            do a = 0, deg
                clmn = left - deg + a
                AB(kl + ku + offset + row - clmn + 1, clmn +1 ) = work(a +1)
            end do
        end do
    end subroutine handle_lhs_derivatives

    integer function find_interval(knots, n_knots, deg, x_val, prev_l, extrapolate) result(interval)
        use, intrinsic :: iso_fortran_env
        implicit none
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

        if (((x_val < tb).or.(x_val > te)).and.(extrapolate==0)) then
            print*, "Some error occured"
            print*, knots
            print*, x_val
            interval = -1
            RETURN
        end if 

        if ((deg < prev_l).and.(prev_l<n)) then
            l = prev_l
        else
            l = deg
        end if

        do while ((x_val < knots(l+1)).and.(l.ne.deg)) !Fortran index
            l = l-1
        end do

        l = l+1

        do while ((x_val>=knots(l+1)).and.(l.ne.n)) !Fortran index
            l = l + 1
        end do

        interval = l-1

        return
    end function find_interval

    subroutine get_intervals(knots, n_knots, deg, x, n, extrapolate, intervals)

        ! in
        integer, INTENT(IN) :: deg, n, n_knots, extrapolate
        real(dp), INTENT(IN) :: knots(n_knots), x(n)

        integer, INTENT(out) :: intervals(n)

        integer :: i, prev_l

        prev_l = deg
        do i=1,n
            intervals(i) = find_interval(knots, n_knots, deg, x(i), prev_l, extrapolate)
        end do
        RETURN
    end subroutine get_intervals

    subroutine evaluate_bspline(knots, n_knots, x_new, y_new, n_new, deg, nu, coeffs, n)
        use, intrinsic :: iso_fortran_env
        implicit none
        integer, parameter :: dp = REAL64

        ! in
        integer, INTENT(IN) :: deg, n, n_knots, n_new
        integer, INTENT(IN) :: nu ! Order of derivative to evaluate (0 if interp value)
        real(dp), intent(IN) :: x_new(n_new), knots(n_knots), coeffs(n)
        ! out
        real(dp), INTENT(OUT) :: y_new(n_new)
        ! internal
        integer :: interval, a, i
        real(dp), ALLOCATABLE :: work(:)

        ! allocate work
        ALLOCATE(work(2*deg+2))

        interval = deg
        ! loop over new x_values
        do i = 1, n_new
            interval = find_interval(knots, n_knots, deg, x_new(i), interval, 0)
            ! Could check here if interval <0
            if (interval < 0) then
                print*, "Interval Error in evaluate. x_new monotonic?"
            end if
            ! Evaluate k+1 which are non-zero in the interval
            call deBoor_D(knots, n_knots, x_new(i), deg, interval, nu, work)

            y_new(i) = 0
            do a = 0,deg
                y_new(i) = y_new(i) + coeffs(interval + a - deg+1)*work(a+1)
            end do
        end do
        RETURN
    end subroutine evaluate_bspline

    subroutine build_LHS(x, n, knots, n_knots, deg, derivs_known_l, derivs_known_r, nleft, nright, &
                            derivs_l_ord, derivs_r_ord, AB, ab_rows, ab_cols, nt, kl, ku, ierr)
        implicit none
        ! in
        integer,intent(in) :: n, n_knots, deg
        real(dp), intent(in) :: x(n), knots(n_knots)
        integer, INTENT(IN) :: derivs_known_l, derivs_known_r, nleft, nright
        integer, INTENT(IN) :: derivs_l_ord(derivs_known_l)
        integer, INTENT(IN) :: derivs_r_ord(derivs_known_r)

        !out
        integer, INTENT(OUT) :: nt, ku, kl, ab_rows, ab_cols, ierr
        real(dp), INTENT(OUT), ALLOCATABLE :: AB(:,:)

        ! Set up left hand side
        nt = n_knots - deg - 1
        kl = deg
        ku = deg
        ab_rows =2*kl + ku +1
        ab_cols = nt
        allocate(AB(ab_rows, ab_cols))
        AB = 0.0d0

        call build_colloc_matrix(x, n, knots, n_knots, deg, AB, ab_rows, ab_cols, nleft)
        if (nleft>0) then
            call handle_lhs_derivatives(knots, n_knots, deg, x(1), AB, ab_rows, ab_cols, kl, ku, &
            derivs_l_ord, derivs_known_l, 0)
        end if
        if (nright>0) then
            call handle_lhs_derivatives(knots, n_knots, deg, x(n), AB, ab_rows, ab_cols, kl, ku, &
            derivs_r_ord, derivs_known_r, nt-nright)
        end if
    end subroutine build_LHS

    subroutine pre_factorize_LHS(AB, ab_rows, ab_cols, nt, kl, ku, ipiv, info)

        ! in
        integer, INTENT(IN) :: ab_rows, ab_cols
        integer, INTENT(IN) :: nt, kl, ku

        ! inplace
        real(dp), INTENT(INOUT) :: AB(ab_rows, ab_rows)

        ! out
        integer, Allocatable, INTENT(OUT) :: ipiv(:) ! TODO is this even needed afterwards?
        integer, INTENT(OUT) :: info

        ! Compute the LU factorization of the band matrix A.

        ALLOCATE(ipiv(nt))

        CALL dgbtrf(nt, nt, kl, ku, AB, 2*kl+ku+1, ipiv, info )
        if (info.ne.0) then
            !error?
            print*, "Error in pre_factorize"
        end if
    end subroutine pre_factorize_LHS

    subroutine solve_coeffcients(nt, kl, ku, ab_rows, ab_cols, LHS, RHS, ipiv, info)
        ! Now we are ready to solve the matrix!
        ! Using LAPACKs' GBSV routine
        ! On exit, rhs is overwritten by the solution coeffs

        ! in
        integer, INTENT(IN) :: nt, kl, ku, ab_rows, ab_cols
        real(dp), intent(in) :: LHS(ab_rows, ab_cols), RHS(nt)

        ! inout
        integer, intent(INOUT) :: ipiv(nt) ! must be inout, but not actually used
        integer, INTENT(INOUT) :: info
    
        ! result (= coeffs) get written into RHS -> thus copy first
        ! LHS will also be overwritten! TODO
        call dgbtrs('No transpose', nt, kl, ku, 1, LHS, 2*kl+ku+1, ipiv, RHS, nt, info)

        return
    end subroutine solve_coeffcients

    subroutine prepare_RHS(nt, nleft, nright, derivs_l_val, derivs_r_val, rhs)
        ! in
        integer, INTENT(IN) :: nt, nleft, nright
        real(dp), INTENT(IN) :: derivs_l_val(nleft), derivs_r_val(nright)

        ! out
        real(dp), ALLOCATABLE, INTENT(OUT) :: rhs(:)


        ALLOCATE(rhs(nt))

        if (nleft>0) then
            rhs(1:nleft) = derivs_l_val
        end if
        if (nright>0) then
            rhs(nt - nright +1: nt) = derivs_l_val
        end if
        ! Assign NaNs in all places where y needs to go, so that if the array is used without assignment, an error will be thrown.
        rhs(nleft + 1 : nt - nright) = ieee_value(real(dp), ieee_signaling_nan)
        RETURN
    end subroutine prepare_RHS

    subroutine assign_y_to_RHS(y, rhs, nt, nleft, nright)
        ! in
        integer, INTENT(IN) :: nt, nleft, nright
        real(dp), INTENT(IN) :: y(nt - nleft - nright)
        ! inout
        real(dp), INTENT(INOUT) :: rhs(nt)

        rhs(nleft+1:nt-nright) = y ! assign y values
        RETURN
    end subroutine assign_y_to_RHS

    subroutine pre_eval_bspline(knots, n_knots, x_new, n_new, deg, nu, intervals, deBoor_matrix)
        integer, INTENT(in) :: n_knots, n_new, deg, nu
        real(dp), INTENT(in) :: knots(n_knots), x_new(n_new)

        integer, INTENT(out) :: intervals(n_new)
        real(dp), intent(out), ALLOCATABLE :: deBoor_matrix(:,:)

        real(dp), ALLOCATABLE :: work(:)

        ! internal
        integer :: i


        ! calculate interval vector
        call get_intervals(knots, n_knots, deg, x_new, n_new, 0, intervals)

        ALLOCATE(deBoor_matrix(n_new, 2*deg+2))
        ALLOCATE(work(2*deg+2))
        do i = 1, n_new
            call deBoor_D(knots, n_knots, x_new(i), deg, intervals(i), nu, work)
            deBoor_matrix(i,:) = work
        end do
        RETURN
    end subroutine pre_eval_bspline

    subroutine eval_bspline_fast(deg, coeffs, intervals, deBoor_matrix, n_new, y_new)

        ! in
        integer, INTENT(in) :: deg, n_new
        integer, intent(in) :: intervals(n_new)
        real(dp), INTENT(in) :: coeffs(n_new), deBoor_matrix(n_new, 2*deg+2)

        ! out
        real(dp), intent(out) :: y_new(n_new)

        ! internal
        integer :: i, a

        y_new = 0.0d0
        do i = 1, n_new
            
            do a = 0, deg
                y_new(i) = y_new(i) + coeffs(intervals(i) +a -deg+1)*deBoor_matrix(i, a+1)
            end do
        end do
        RETURN
    end subroutine eval_bspline_fast

end module interpolation
