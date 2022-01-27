! Created by Alexander M. Imre on 04.12.21.
!
! Fast B-Spline Interpolation
! A part of ViPErLEED
!
! This module is designed to enable fast interpolation tailored to the requirements of TensorLEED calculations in ViPErLEED.
! At its heart, it performs a spline interpolation in the b-spline basis. The structure search of TensorLEED calculations requires frequent interpolation from
! a less dense origin grid to a more dense target grid. Unlike more general interpolation libraries, this module leverages that fact, by separating the interpolation
! into a pre-evaluation that only depends on the abcissaes of the grids, and an actual interpolation between the ordinate values.
! Large parts of the required calculations, e.g. the setup of the entire left hand side of the matrix equation, thus need only be performed once and can then be reused
! to interpolate an aribtrary number of datasets on the same grid.
! 
! This module is a part of ViPErLEED and is distributed under the same license.


! Error codes:
! ViPErLEED Interpolation Error codes start with 6
! 0 .. no error 
! 1-10 .. ?
! 11-19 .. error raised while performing checks on input data
! 11 ... origin grid not sorted
! 12 ... target grid not sorted


module interpolation
    use, intrinsic :: ieee_arithmetic
    implicit none
    
    contains

! *****************************************************************************************
! Routines to get array sizes

    pure subroutine get_array_sizes(n, deg, n_knots, nt, LHS_rows)
        ! Returns sizes of arrays used various parts of the interpolation. This is necessary since F2PY 
        ! wrapping requires explicit array shapes and does not allow for allocatable or assumed shape arrays.
        integer, INTENT(IN) :: n, deg

        integer, INTENT(OUT) :: n_knots
        integer, INTENT(OUT) :: nt
        integer, INTENT(OUT) :: LHS_rows

        integer :: kl, ku

        n_knots = n + 2*deg
        nt = n_knots - deg - 1
        kl = deg
        ku = deg
        LHS_rows = 2*kl + ku + 1
        RETURN
    end subroutine get_array_sizes

    pure subroutine de_Boor_size(n_new, deg, de_Boor_dim)
        ! Returns the size of the de Boor coefficient matrix
        integer,intent(in) :: n_new
        integer,intent(in) :: deg

        integer,INTENT(OUT) :: de_Boor_dim(2)

        de_Boor_dim(1) = 2*deg + 2
        de_Boor_dim(2) = n_new
        RETURN
    end subroutine de_Boor_size

! *****************************************************************************************
! Input grid pre-evaluation and spline coefficients

    subroutine single_calc_spline(&
        n, x, y, &
        deg, &
        n_knots, nt, LHS_rows, &
        knots, &
        coeffs, &
        ierr)

        ! IN
        integer, INTENT(IN) :: deg              !! order of spline
        integer, INTENT(IN) :: n        !! # points in
        real(8), INTENT(IN):: x(n)             !! origin grid values
        real(8), INTENT(IN):: y(n)
        integer, INTENT(IN) :: n_knots, nt, LHS_rows

        ! OUT
        real(8), INTENT(OUT)  :: knots(n_knots)                 !! knot points
        real(8), INTENT(OUT)  :: coeffs(nt)
        integer, INTENT(OUT)   :: ierr                    !! integer error code

        ! Internal
        integer  :: nleft, nright
        real(8) :: LHS(LHS_rows, nt)
        real(8) :: RHS (nt)
        integer  :: ipiv(nt)

        ierr = 0

        call pre_eval_input_grid(n, x, deg, &
            n_knots, nt, LHS_rows, &
            nleft, nright, &
            knots, LHS, RHS, ipiv, &
            ierr &
            )

        if (ierr .ne. 0) RETURN

        call calc_spline_with_pre_eval(&
            n, deg, &
            y, &
            n_knots, nt, LHS_rows, &
            knots, &
            nleft, nright, &
            LHS, RHS, ipiv, &
            coeffs, &
            ierr)

    end subroutine single_calc_spline


    subroutine pre_eval_input_grid( &
        n, x, deg, &                                ! actual input arguments
        n_knots, nt, LHS_rows, & ! array sizes of output - formal input arguments
        nleft, nright, &
        knots, &
        LHS, &
        RHS, &
        ipiv, &
        ierr &
        )
        ! prepares LHS and RHS of the matrix equation for the spline coefficiencts
        ! Used in Structure search where many theoretical intensities need to be interpolated
        ! from the same starting grid.
        ! If you use this version, you should first call subroutine get_array_sizes to get the sizes
        ! of the arrays to expect returned.

        ! IN
        integer, INTENT(IN) :: n
        real(8), INTENT(IN):: x(n)             !! origin grid values
        integer, INTENT(IN) :: deg

        integer, INTENT(IN) :: n_knots
        integer, INTENT(IN) :: nt ! Size of LHS_cols, RHS_cols, coeffs and ipiv
        integer, INTENT(IN) :: LHS_rows
        

        ! OUT
        integer, INTENT(OUT)  :: nleft, nright !! number of knots for bc left and right
        real(8), INTENT(OUT) :: knots(n_knots) !! knot points
        real(8), INTENT(OUT) :: LHS(LHS_rows, nt)
        real(8), INTENT(OUT) :: RHS (nt)
        integer, INTENT(OUT) :: ipiv(nt)
        
        integer, INTENT(OUT)   :: ierr                    !! integer error code

        ! Internal
        integer                             :: kl, ku           !! matrix shape specifiers - in our case always equal deg
        integer                             :: LAPACK_info
        integer                             :: derivs_known_l, derivs_known_r
        integer, ALLOCATABLE                :: derivs_l_ord(:), derivs_r_ord(:)
        real(8), ALLOCATABLE               :: derivs_l_val(:), derivs_r_val(:)

        ierr = 0

        ! Compute knot vector and size
        call get_natural_knots(x, n, deg, n_knots, knots)

        ! Get the correct derivatives at the edge for the degree givennatural boundary conditions
        call prepare_derivs_nat_bc(deg, derivs_known_l, derivs_known_r, nleft, nright,&
        derivs_l_ord, derivs_r_ord, derivs_l_val, derivs_r_val)

        ! Using the derivatives and the x_values, we build up the LHS of the equation
        call build_LHS(n, x, n_knots, knots, deg, derivs_known_l, derivs_known_r, nleft, nright, &
            derivs_l_ord, derivs_r_ord, nt, LHS_rows, LHS)

        kl = deg
        ku = deg

        ! with this done, we can now prefactorize the LHS
        call pre_factorize_LHS(LHS, LHS_rows, nt, kl, ku, ipiv, LAPACK_info)

        if (LAPACK_info.ne.0) then
            ! Error in pre-factorization
            ierr = 621
            RETURN
        end if

        ! Now set up the RHS
        call prepare_RHS(nt, nleft, nright, derivs_l_val, derivs_r_val, RHS)

    end subroutine pre_eval_input_grid


    subroutine calc_spline_with_pre_eval(&
        n, deg, &
        y, &
        n_knots, nt, LHS_rows, &
        knots, &
        nleft, nright, &
        LHS, RHS, ipiv, &
        coeffs, &
        ierr)
        ! Calculate spline coefficients given pre-calculated LHS and RHS from routine pre_eval_input_grid

        integer, INTENT(IN):: n
        real(8), INTENT(IN):: y(n)
        integer, INTENT(IN) :: deg

        integer, INTENT(IN)   :: n_knots, nt, LHS_rows
        real(8), INTENT(IN)  :: knots(n_knots)
        integer, INTENT(IN)   :: nleft, nright
        real(8), INTENT(IN)  :: LHS(LHS_rows, nt)
        real(8), INTENT(IN)  :: RHS(nt)
        integer, INTENT(IN)  :: ipiv(nt)

        real(8), INTENT(OUT) :: coeffs(nt)
        integer, INTENT(OUT) :: ierr

        ! Internal
        integer :: kl, ku, LAPACK_info

        ierr = 0

        call assign_y_to_RHS(y, RHS, nt, nleft, nright, coeffs)

        kl = deg
        ku = deg

        call solve_coeffcients(nt, kl, ku, LHS_rows, LHS, coeffs, ipiv, LAPACK_info)
        if ((LAPACK_info >0).or.(LAPACK_info<0)) then
            ierr = 631 ! LAPACK Error
        end if
        RETURN
    end subroutine calc_spline_with_pre_eval

! *****************************************************************************************
! Output grid pre-evaluation and interpolation

    subroutine single_interpolate_coeffs_to_grid(n_knots, knots, nt, coeffs, deg, &
        n_new, x_new, &
        nu, &
        y_new)

        ! IN
        integer, INTENT(IN) :: deg              !! order of spline
        integer, INTENT(IN) :: n_knots, nt, n_new         !! # points in, out
        real(8), INTENT(IN):: x_new(n_new)     !! result grid values
        real(8), INTENT(IN):: knots(n_knots)
        real(8), INTENT(IN):: coeffs(nt)
        integer, INTENT(IN) :: nu

        !OUT
        real(8), INTENT(OUT):: y_new(n_new)
        ! Internal
        integer              :: intervals(n_new)        !! 1D integer array containg information, which spline interval the new grid points fall into.
        real(8), ALLOCATABLE  :: deBoor_matrix(:,:)      !! pre-evaluated deBoor matrix, required for evaluation of values on the new grid.
        integer :: de_Boor_dim(2)

        ! calculate interval vector for x_new
        call get_intervals(n_knots, knots, n_new, x_new, deg, intervals)

        call de_Boor_size(n_new, deg, de_Boor_dim)

        ALLOCATE(deBoor_matrix(de_Boor_dim(1), de_Boor_dim(2)))

        ! Build up the deBoor Matrix for x_new
        call calc_deBoor(n_knots, knots, n_new, x_new, deg, nu, intervals, &
            de_Boor_dim(1), de_Boor_dim(2), deBoor_matrix)

        call eval_bspline_fast(deg, nt, n_new, coeffs, intervals, deBoor_matrix, y_new)
        
    end subroutine single_interpolate_coeffs_to_grid


    subroutine get_intervals(n_knots, knots, n_new, x_new, deg, intervals)

        ! in
        integer, INTENT(IN) :: deg, n_new, n_knots
        real(8), INTENT(IN) :: knots(n_knots), x_new(n_new)

        integer, INTENT(out) :: intervals(n_new)

        integer :: extrapolate
        integer :: i, prev_l

        extrapolate = 0
        prev_l = deg
        do i=1,n_new
            intervals(i) = find_interval(knots, n_knots, deg, x_new(i), prev_l, extrapolate)
        end do
        RETURN
    end subroutine get_intervals


    subroutine calc_deBoor(n_knots, knots, n_new, x_new, deg, nu, intervals, deBoor_rows, deBoor_cols, deBoor_matrix)
        ! Calculates the deBoor_matrix
        ! IN
        integer, INTENT(in) :: n_knots, n_new, deg, nu
        real(8), INTENT(in) :: knots(n_knots), x_new(n_new)
        integer, INTENT(in) :: intervals(n_new)
        integer, INTENT(IN) :: deBoor_rows, deBoor_cols
        
        ! OUT
        real(8), intent(out) :: deBoor_matrix(deBoor_rows, deBoor_cols)

        ! internal
        integer :: i

        do i = 1, n_new
            call deBoor_D(knots, n_knots, x_new(i), deg, intervals(i), nu, deBoor_matrix(:,i))
        end do
        RETURN
    end subroutine calc_deBoor

! *****************************************************************************************
! Fast Evaluation of the B-Spline using prepared values

    pure subroutine eval_bspline_fast(&
        deg, nt, n_new, &
        coeffs, &
        intervals, &
        deBoor_matrix, &
        y_new &
        )
        ! TODO: see if this can be sped up. This is the rate limiting procedure

        ! in
        integer, INTENT(in) :: deg, n_new
        integer, INTENT(IN) :: nt
        integer, intent(in) :: intervals(n_new)
        real(8), INTENT(in) :: coeffs(nt)
        real(8), INTENT(IN):: deBoor_matrix(2*deg+2, n_new)

        ! out
        real(8), intent(out) :: y_new(n_new)

        ! internal
        integer :: i ! loop variable

        y_new = 0.0d0
        do concurrent (i = 1:n_new)
            ! whole array operations like the one below are faster than manual loops
            y_new(i) = sum(coeffs(intervals(i)-deg+1:intervals(i)+1)*deBoor_matrix(1:deg+1,i))
        end do
        RETURN
    end subroutine eval_bspline_fast

    pure subroutine eval_y_pendry_fast(&
        deg, nt, n_new, &
        coeffs, &
        intervals, &
        v0i, &
        deBoor_matrix_deriv_0, deBoor_matrix_deriv_1, &
        y_func &
        )
        ! evaluates the B-spline in the same way as eval_bspline_fast, but does it at once for
        ! the intensity and the derivative and returns the Pendry Y function directly.

        ! in
        integer, INTENT(in) :: deg, n_new, nt
        integer, intent(in) :: intervals(n_new)
        real(8), INTENT(in) :: coeffs(nt)
        real(8), INTENT(in) :: deBoor_matrix_deriv_0(2*deg+2, n_new)
        real(8), INTENT(in) :: deBoor_matrix_deriv_1(2*deg+2, n_new)
        real(8), INTENT(IN) :: v0i
    
        ! OUT
        real(8), INTENT(OUT) :: y_func(n_new)

        ! Internal
        integer i
        real(8) :: intensity(n_new), derivative(n_new)
        
        y_func = 0.0d0
        do concurrent (i=1:n_new)
            intensity(i) = sum(coeffs(intervals(i)-deg+1:intervals(i)+1)*deBoor_matrix_deriv_0(1:deg+1,i))
            derivative(i) = sum(coeffs(intervals(i)-deg+1:intervals(i)+1)*deBoor_matrix_deriv_1(1:deg+1,i))
            y_func(i) = intensity(i)*derivative(i) /(intensity(i)*intensity(i) + v0i*v0i*derivative(i)*derivative(i))
        end do

    end subroutine eval_y_pendry_fast

! *****************************************************************************************
! Utility procedures

    subroutine perform_checks(n, x, n_new, x_new, ierr)
        implicit none
        integer, intent(in) :: n, n_new
        real(8), INTENT(IN) :: x(n), x_new(n_new)
        integer,intent(out) ::  ierr

        integer :: i ! loop var

        ierr = 0

        ! check if x values are sorted
        do concurrent (i = 2: n)
            if (x(i).le.x(i-1)) then
                ! origin grid not sorted
                ierr = 611
            end if
        end do

        do concurrent (i = 2: n_new)
            if (x_new(i).le.x_new(i-1)) then
                ! target grid not sorted
                ierr = 612
            end if
        end do

        RETURN
    end subroutine perform_checks


    pure subroutine equidist_to_array(n, x_start, x_step, x)
        integer, INTENT(IN)  :: n
        real(8), INTENT(IN) :: x_start, x_step

        real(8), INTENT(OUT) :: x(n)

        integer :: i

        do i = 1, n
            x(i) = x_start + x_step*(i-1)
        end do
        RETURN
    end subroutine equidist_to_array

! *****************************************************************************************
! Handeling of knots and boundary conditions
   
    subroutine get_natural_knots(x, n, deg, n_knots, knots)
        ! in
        integer,intent(in) :: n, deg
        real(8), INTENT(IN) :: x(n)
        integer, INTENT(IN) :: n_knots

        ! out
        real(8), intent(out) ::  knots(n_knots)


        knots(1       : deg     ) = x(1)
        knots(1+deg   : n+deg   ) = x(:)
        knots(n+deg+1 : n_knots ) = x(n)
        RETURN
    end subroutine get_natural_knots


    subroutine prepare_derivs_nat_bc(deg, derivs_known_l, derivs_known_r, nleft, nright,&
        derivs_l_ord, derivs_r_ord, derivs_l_val, derivs_r_val)
        integer, INTENT(IN) :: deg
        

        integer, INTENT(OUT) :: derivs_known_l, derivs_known_r, nleft, nright
        integer, ALLOCATABLE, INTENT(OUT) :: derivs_l_ord(:), derivs_r_ord(:)
        real(8), ALLOCATABLE, INTENT(OUT) :: derivs_l_val(:), derivs_r_val(:)

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

! *****************************************************************************************
! Internals for building LHS

    subroutine build_colloc_matrix(x, n, knots, n_knots, deg, AB, LHS_rows, LHS_cols, offset)
        ! -> bspl._colloc from scipy

        ! in
        integer, INTENT(IN) :: deg, n, n_knots, LHS_cols, LHS_rows
        real(8), intent(IN) :: x(n), knots(n_knots)
        ! out
        real(8), INTENT(INOUT) :: AB(LHS_rows,LHS_cols)
        ! internal
        integer :: left, j, a, kl, ku, clmn, nt, offset
        real(8) :: x_val
        real(8), ALLOCATABLE :: work(:)

        ku = deg
        kl = deg
        nt = n_knots - deg - 1

        allocate(work(2*deg+2))

        left = deg
        ! TODO fix indices
        do j=0, n-1
            x_val = x(j+1)
            left = find_interval(knots, n_knots, deg, x_val, left, 0) ! 0 means extrapolate = False

            call deBoor_D(knots, n_knots, x_val, deg, left, 0, work)

            do a = 0, deg
                clmn = left - deg + a
                AB(kl + ku + j + offset - clmn + 1, clmn +1 ) = work(a +1)
            end do
        end do

    end subroutine build_colloc_matrix


    subroutine handle_lhs_derivatives(knots, n_knots, deg, x_val, AB, LHS_rows, LHS_cols, kl, ku, deriv_ords, n_derivs, offset)
        implicit none

        ! in
        integer, INTENT(IN) :: deg, n_knots, LHS_cols, LHS_rows, offset, kl, ku, n_derivs, deriv_ords(n_derivs)
        real(8), intent(IN) :: knots(n_knots), x_val
        ! out
        real(8), INTENT(INOUT) :: AB(LHS_rows,LHS_cols)
        ! internal
        integer :: a, left, nu, row, clmn

        real(8), ALLOCATABLE :: work(:)
        
        ! derivatives @ x_val
        left = find_interval(knots, n_knots, deg, x_val, deg, 0)

        ! allocate work array
        ALLOCATE(work(2*deg+2))
        do concurrent (row=0: n_derivs-1)
            nu = deriv_ords(row+1)
            call deBoor_D(knots, n_knots, x_val, deg, left, nu, work)

            do a = 0, deg
                clmn = left - deg + a
                AB(kl + ku + offset + row - clmn + 1, clmn +1 ) = work(a +1)
            end do
        end do
    end subroutine handle_lhs_derivatives

! *****************************************************************************************
! de Boor B-spline evaluation

    pure subroutine deBoor_D(knots, n_knots, x, deg, ell, m, result)
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

        real(8), INTENT(IN) :: x, knots(n_knots)

        real(8), INTENT(INOUT) :: result(deg+1)

        !internal
        real(8) :: xb, xa, w
        integer :: ind, j, n, i
        real(8) :: hh(deg)

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

! *****************************************************************************************
! Assigning grid points to spline intervals

    pure function find_interval(knots, n_knots, deg, x_val, prev_l, extrapolate) result(interval)

        integer, INTENT(IN) :: n_knots, deg, extrapolate, prev_l
        real(8), INTENT(IN) :: knots(n_knots), x_val

        integer interval

        ! internal
        integer :: l, n
        real(8) :: tb, te

        ! TODO Fortran indices fix!
        n = n_knots - deg - 1
        tb = knots(deg+1)   !Fortran index
        te = knots(n+1)     !Fortran index

        ! reference code checks for NaN in x_val -> skipped here

        if (((x_val < tb).or.(x_val > te)).and.(extrapolate==0)) then
            ! TODO implement ierr here!
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


! *****************************************************************************************
! Preparing the left hand side (LHS) of the equation. The LHs holds the original grid (knots) and info on the boundary conditions.
! As long as the origin grid and boundary conditions stay the same, the LHS is unchanged and can be pre-evaluated.    

    subroutine build_LHS(n, x, n_knots, knots, deg, derivs_known_l, derivs_known_r, nleft, nright, &
                            derivs_l_ord, derivs_r_ord, nt, LHS_rows, LHS)
        implicit none
        ! in
        integer,intent(in) :: n, n_knots, deg
        real(8), intent(in) :: x(n), knots(n_knots)
        integer, INTENT(IN) :: derivs_known_l, derivs_known_r, nleft, nright
        integer, INTENT(IN) :: derivs_l_ord(derivs_known_l)
        integer, INTENT(IN) :: derivs_r_ord(derivs_known_r)
        integer, INTENT(IN) :: nt, LHS_rows

        !out
        real(8), INTENT(OUT) :: LHS(LHS_rows, nt)

        ! Internal
        integer :: kl, ku

        ! Set up left hand side
        LHS = 0.0d0

        ku = deg
        kl = deg

        call build_colloc_matrix(x, n, knots, n_knots, deg, LHS, LHS_rows, nt, nleft)
        if (nleft>0) then
            call handle_lhs_derivatives(knots, n_knots, deg, x(1), LHS, LHS_rows, nt, kl, ku, &
            derivs_l_ord, derivs_known_l, 0)
        end if
        if (nright>0) then
            call handle_lhs_derivatives(knots, n_knots, deg, x(n), LHS, LHS_rows, nt, kl, ku, &
            derivs_r_ord, derivs_known_r, nt-nright)
        end if
    end subroutine build_LHS


    subroutine pre_factorize_LHS(LHS, LHS_rows, nt, kl, ku, ipiv, info)

        ! in
        integer, INTENT(IN) :: LHS_rows
        integer, INTENT(IN) :: nt, kl, ku

        ! inplace
        real(8), INTENT(INOUT) :: LHS(LHS_rows, nt)

        ! out
        integer, INTENT(OUT) :: ipiv(nt) ! TODO is this even needed afterwards?
        integer, INTENT(OUT) :: info

        ! Compute the LU factorization of the band matrix A.

        CALL dgbtrf(nt, nt, kl, ku, LHS, 2*kl+ku+1, ipiv, info)

    end subroutine pre_factorize_LHS

! *****************************************************************************************
! Deal with right hand side (RHS) of the equation
! The RHS holds the values at the knot points, i.e. in our case the y values and derivatives at the edge

    subroutine prepare_RHS(nt, nleft, nright, derivs_l_val, derivs_r_val, rhs)
        ! Allocates the right size for the RHS and fills in derivatives at the edge.
        ! The other values are set to NaN and need to be filled in by the y values for actual interpolation
        ! in
        integer, INTENT(IN) :: nt, nleft, nright
        real(8), INTENT(IN) :: derivs_l_val(nleft), derivs_r_val(nright)

        ! out
        real(8), INTENT(OUT) :: rhs(nt)


        if (nleft>0) then
            rhs(1:nleft) = derivs_l_val
        end if
        if (nright>0) then
            rhs(nt - nright +1: nt) = derivs_r_val
        end if
        ! Assign NaNs in all places where y needs to go, so that if the array is used without assignment, an error will be thrown.
        rhs(nleft + 1 : nt - nright) = ieee_value(real(8), ieee_signaling_nan)
        RETURN
    end subroutine prepare_RHS


    subroutine assign_y_to_RHS(y, rhs_prep, nt, nleft, nright, rhs)
        ! in
        integer, INTENT(IN) :: nt, nleft, nright
        real(8), INTENT(IN) :: y(nt - nleft - nright)
        real(8), INTENT(IN) :: rhs_prep(nt)
        ! out
        real(8), INTENT(OUT) :: rhs(nt)

        ! Since rhs will be overwritten by LAPACK, it needs to be copied at some point prior
        rhs = rhs_prep
        rhs(nleft+1:nt-nright) = y ! assign y values
        RETURN
    end subroutine assign_y_to_RHS

! *****************************************************************************************
! Solve for the B-Spline coefficients given prepared LHS and RHS

    subroutine solve_coeffcients(nt, kl, ku, LHS_rows, LHS, RHS, ipiv, info)
        ! Now we are ready to solve the matrix!
        ! Using LAPACKs' GBSV routine
        ! On exit, rhs is overwritten by the solution coeffs

        ! in
        integer, INTENT(IN) :: nt, kl, ku, LHS_rows
        real(8), intent(in) :: LHS(LHS_rows, nt)
        real(8), INTENT(INOUT) :: RHS(nt)

        ! inout
        integer, intent(IN) :: ipiv(nt) ! must be inout, but not actually used
        integer, INTENT(INOUT) :: info
    
        ! result (= coeffs) get written into RHS -> thus copy first
        ! LHS will also be overwritten! TODO
        call dgbtrs('No transpose', nt, kl, ku, 1, LHS, 2*kl+ku+1, ipiv, RHS, nt, info)

        return
    end subroutine solve_coeffcients


end module interpolation
