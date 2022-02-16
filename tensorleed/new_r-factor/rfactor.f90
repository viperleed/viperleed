! Subroutines and functions for R factor calculation in ViPErLEED
!
! v0.1.0
!
! Author: Alexander M. Imre, 2021
! for license info see ViPErLEED Package

! Note: file extentsion .f90 required for f2py compatibility, but newer standard can be used.

module r_factor_new
    use interpolation
    implicit none

    contains

    

! f2py -m rfactor rfactor.f90 -h rfactor.pyf --overwrite-signature --debug-capi
! f2py -c rfactor.pyf rfactor.f90



! The equidistant final grid is defined by the array E_new with length n_E_new.  This is equivalently defined by the set of integers n_E_new, E_new_start, E_new_step.

! Experimental data may come in as arrays 



   
subroutine r_pendry_beam_y(n_E, E_step, y1, y2, id_start_y1, id_start_y2, n_y1, n_y2, V0r_shift, &
                         R_pendry, numerator, denominator, N_overlapping_points)
    ! Innermost function Rfactor_beam:
    ! INPUT: Y1, Y2, Estart1, Estart2, V0rshift -> Y1 and Y2 are already prepared (i.e. same grid steps, I(V) was previously interpolated)
    ! DOES: figure out overlap region, calculate RPe
    ! RETURNS: RPe, Number of overlapping points/ !!! -> not just R pendry for each beam, but numerator & denominator! for further calculation!!
    !

    !###############
    !INPUTS
    !###############

    !f2py intent(in) n_E
    integer, intent(in) :: n_E
    real(8), INTENT(IN) :: E_step
    integer, INTENT(IN) :: id_start_y1, id_start_y2
    integer, INTENT(IN) :: n_y1, n_y2
    !f2py intent(in) y1
    real(8), intent (in), dimension(n_E)   :: y1
    !f2py intent(in) y2
    real(8), intent (in), dimension(n_E)   :: y2
    integer, intent (in)                   :: V0r_shift
    ! outputs
    real(8), intent (out)                  :: R_pendry, numerator, denominator
    integer,  intent (out)              :: N_overlapping_points

    real(8), dimension (:), allocatable         :: y_diff, y_squared_sum

    integer :: id_min, id_max

    ! V0r_shift shifts y2 by some number of grid points -> id_start_y2 = id_start_y2 + V0r_shift
    id_min = max(id_start_y1, id_start_y2 + V0r_shift)
    id_max = min(id_start_y1 + n_y1 - 1, id_start_y2 + n_y2 - 1 + V0r_shift)
    
    N_overlapping_points = id_max - id_min + 1

    allocate(y_diff(N_overlapping_points), y_squared_sum(N_overlapping_points))

    ! Shift it around

    ! Remember that indices are shifted by V0r

    y_diff=y1(id_min : id_max) - y2(id_min - V0r_shift: id_max - V0r_shift) ! difference between Y functions
    y_squared_sum=y1(id_min: id_max)**2 + y2(id_min - V0r_shift: id_max - V0r_shift)**2
    
    ! caclulate numerator = integral (Y1-Y2)**2 dE
    numerator = trapez_integration_const_dx(y_diff**2, E_step)
    ! calculate denominator = integral (Y1**2+Y2**2) dE
    denominator = trapez_integration_const_dx(y_squared_sum, E_step)

    R_pendry = numerator/denominator

    return

end subroutine r_pendry_beam_y

subroutine r2_beam_intensity(n_E, E_step, int_1, int_2, id_start_1, id_start_2, n_1, n_2, &
    V0r_shift, R2, n_overlapping_points)

    integer, intent(in) :: n_E
    real(8), INTENT(IN) :: E_step
    integer, INTENT(IN) :: id_start_1, id_start_2
    integer, INTENT(IN) :: n_1, n_2
    !f2py intent(in) int_1
    real(8), intent (in), dimension(n_E)   :: int_1
    !f2py intent(in) int_2
    real(8), intent (in), dimension(n_E)   :: int_2
    integer, intent (in)                   :: V0r_shift
    ! outputs
    real(8), intent (out)                  :: R2
    integer,  intent (out)              :: n_overlapping_points

    real(8), dimension (:), allocatable         :: int_diff, y_squared_sum


    integer :: id_min, id_max
    real(8) :: A, c

    id_min = max(id_start_1, id_start_2 + V0r_shift)
    id_max = min(id_start_1 + n_1 - 1, id_start_2 + n_2 - 1 + V0r_shift)

    n_overlapping_points = id_max - id_min + 1

    
    A = 1/trapez_integration_const_dx(int_1(id_min: id_max)**2, E_step)

    c = trapez_integration_const_dx(int_1(id_min: id_max), E_step) / &
        trapez_integration_const_dx(int_2(id_min -V0r_shift : id_max - V0r_shift), E_step)

    R2 = A*trapez_integration_const_dx( ( int_1(id_min: id_max) - &
        c*int_2(id_min -V0r_shift : id_max - V0r_shift) )**2, &
        E_step ) / 2
    ! not entirely sure where this factor of 1/2 comes from, but it was implemented in TensErLEED


end subroutine r2_beam_intensity


subroutine r_pendry_beamset_V0r_opt_on_grid( &
        range, start_guess, fast_search, best_V0r_step, best_V0r, best_R, n_evaluated, &       
        n_E, E_step, n_beams, y1, y2, id_start_y1, id_start_y2, n_y1, n_y2, &
        r_pendry_beams, numerators, denominators, n_overlapping_points, &
        ierr, &
        tol_R, tol_R_2, max_fit_range)
    
    implicit none
    integer, parameter :: n_start_guesses = 3
    ! IN
    integer, intent(in) :: range(2)
    integer, INTENT(IN) :: start_guess(n_start_guesses)

    ! OUT
    logical, INTENT(INOUT) :: fast_search
    integer, INTENT(OUT) :: best_V0r_step
    real(8), INTENT(OUT) :: best_R, best_V0r ! best_V0r as real
    integer, INTENT(OUT) :: n_evaluated
    integer, INTENT(OUT) :: ierr

    ! Internal variables used for V0r optimization
    real(8), ALLOCATABLE :: R_V0r(:), weights(:)
    integer :: n_steps
    real(8) :: fit_V0r_id
    real(8) :: fit_x_min, fit_y_min, fit_curvature
    integer, ALLOCATABLE :: V0r_step(:)
    real(8), ALLOCATABLE :: V0r_step_real(:)
    real(8) :: para_coeff(3)
    logical, ALLOCATABLE :: evaluated(:)

    ! Optionals
    !f2py real(8), optional :: tol_R = 0.00001
    real(8), INTENT(IN) :: tol_R, tol_R_2
    integer, INTENT(IN) :: max_fit_range

    integer :: i, j

    integer :: fit_range, min_fit_range
    integer :: next_step, closest_step
    integer :: worst_step
    integer :: best_id
    integer :: left, right
    integer :: closest_step_id
    real(8) :: R2
    integer :: ierr_tmp

    ! holder_arrays for passout
    integer, ALLocatable :: N_overlapping_points_V0r(:,:)
    real(8), ALLocatable :: r_pendry_beams_V0r(:,:)
    real(8), ALLOCATABLE :: numerators_V0r(:,:), denominators_V0r(:,:)


    ! Pass-through I/O for r_pendry_beamset
    !f2py integer, hidden, intent(in), depend(E_start1, E_start2), check(len(E_start1)==len(E_start2)) :: nr_beams=len(E_start1)
    integer :: n_E ! number of energy steps
    integer n_beams
    !f2py real(8), intent(in) :: E_step
    real(8) :: E_step
    !f2py real intent(in) :: y1 (n_E, n_beams), y2 (n_E, n_beams)
    real(8), intent(in) :: y1 (n_E, n_beams), y2 (n_E, n_beams)
    !f2py integer, intent(in) id_start_y1
    !f2py integer, intent(in) id_start_y2
    integer, intent(in)    :: id_start_y1(n_beams) ,id_start_y2(n_beams)
    !f2py integer, intent(in) n_y1
    !f2py integer, intent(in) n_y2
    integer, intent(in)    :: n_y1(n_beams) ,n_y2(n_beams)
    !f2py integer, intent(out) :: n_overlapping_points
    integer, intent(out) :: n_overlapping_points(n_beams)
    !f2py real(8), intent(out) :: r_pendry_beams(n_beams), r_pendry_weighted
    real(8), intent(out) :: r_pendry_beams(n_beams)
    real(8), intent(out) :: numerators(n_beams), denominators(n_beams)


    ! *****************************************************************************************

    n_steps = range(2) - range(1) + 1
    
    if (n_steps .le. 5) then
        ierr = 851
        RETURN
    end if


    ALLOCATE(V0r_step(n_steps), V0r_step_real(n_steps))
    ALLOCATE(R_V0r(n_steps), weights(n_steps))
    ! Below arrays are allocated with beam index in first position, because step changes more slowly!
    ALLOCATE(r_pendry_beams_V0r(n_beams, n_steps), N_overlapping_points_V0r(n_beams, n_steps))
    ALLOCATE(numerators_V0r(n_beams, n_steps), denominators_V0r(n_beams, n_steps))
    ALLOCATE(evaluated(n_steps))

    do i = 1, n_steps
        V0r_step(i) = range(1) + i - 1
        V0r_step_real(i) = range(1) + i - 1
    end do

    ! Assign NaNs to unknown R_V0r and impossibly large value to best_R
    R_V0r = 5d3

    best_R = 5d3
    weights = 0
    evaluated = .False.

    n_evaluated = 0
    fit_range = max_fit_range
    min_fit_range = max(max_fit_range - 6, 5)

    next_step = start_guess(1) - range(1) + 1
    do while(fast_search)

        !print*, "Evaluated:", n_evaluated, "next step", next_step
        call r_pendry_beamset_y(n_E, E_step, n_beams, y1, y2, &
            id_start_y1, id_start_y2, n_y1, n_y2, &
            V0r_step(next_step), & 
            R_V0r(next_step), &
            r_pendry_beams_V0r(:, next_step), &
            numerators_V0r(:, next_step), denominators_V0r(:, next_step), &
            N_overlapping_points_V0r(:, next_step), &
            ierr)
        ! Pass out possible R factor error message
        if (ierr .ne. 0) RETURN
        print*, ""
        print*, "Smart step:", V0r_step(next_step), "R = ", R_V0r(next_step)


        ! Calculate next point
        call update_best_V0r(n_steps, n_beams, best_R, best_id, R_V0r(next_step), next_step, &
            N_overlapping_points_V0r, N_overlapping_points, &
            numerators_V0r, numerators, &
            denominators_V0r, denominators, &
            r_pendry_beams_V0r, r_pendry_beams)

        evaluated(next_step) = .True.
        n_evaluated = n_evaluated + 1
        weights(next_step) = 1

        ! Do 3 guesses initially
        if (n_evaluated < 3) then
            next_step = start_guess(n_evaluated + 1) - range(1) +1
            CYCLE
        else if (n_evaluated == 3) then
            call parabola_lsq_fit(n_steps, V0r_step_real, R_V0r, weights, para_coeff, ierr)
            if (ierr .ne. 0) RETURN ! Error
            fit_curvature = 2*para_coeff(1)
            fit_x_min = -para_coeff(2)/para_coeff(1)/2
            closest_step = NINT(fit_x_min) - (range(1)) + 1

            
        else if (n_evaluated == n_steps) then
            ierr = 852
            exit
        else if (n_evaluated > 3) then
            weights = 0
            weights(next_step) = 1
        end if


        ! Decide how to proceed:
        ! *****************************************************************************************

        ! Everything withing fit_range of closest_step is taken for parabola

        do i = 0, fit_range
            if ((closest_step - i < 1) .or. (closest_step + i > n_steps)) then
                ierr = 854
                fast_search = .False.
                print*, "outside range!"
                EXIT
            end if
            if (evaluated(closest_step - i)) then
                weights(closest_step - i) = 1
            end if
            if (evaluated(closest_step + i)) then
                weights(closest_step + i) = 1
            end if
        end do

        if (.not. fast_search) EXIT

        !print*, "SUM weights", sum(weights), n_evaluated

        if (sum(weights) < 3.9) then
            if (.not. evaluated(closest_step)) then
                next_step = closest_step
                CYCLE
            end if
            do i = 1, fit_range
                if (.not. evaluated(closest_step - i)) then
                    next_step = closest_step - i
                    EXIT
                end if
                if (.not. evaluated(closest_step + i)) then
                    next_step = closest_step + i
                    EXIT
                end if
            end do
            CYCLE
        else 
            !Fit parabola to values
            call parabola_lsq_fit(n_steps, V0r_step_real, R_V0r, weights, para_coeff, ierr)
            if (ierr .ne. 0) RETURN ! Error
    
            ! perform parabola fit -> fit_V0r_id, fit_R, curvature,
            fit_curvature = 2*para_coeff(1)
            fit_x_min = -para_coeff(2)/para_coeff(1)/2
            fit_y_min = para_coeff(3) - para_coeff(2)**2/para_coeff(1)/4
    
            print*, "Expecting minimum of R =", fit_y_min, "at V0r = ", fit_x_min*E_step
    
    
            ! decide if fast search possible
            closest_step = NINT(fit_x_min) - (range(1)) + 1
            
            if ((closest_step .le. 1) .or. (closest_step .ge. n_steps) .or. (fit_curvature .le. 0.005)) then
                ! Parabola fit doesnt work...
                fast_search = .False.
                print*, "Parabola is futile!"
                print*, "curvature", fit_curvature
                exit
            end if
    
            ! Have we calculated the point closest to the predicted minimum yet?
            if (.not. evaluated(closest_step)) then
                next_step = closest_step
                CYCLE
            end if
        end if

        R2 = parabola_R_squared(n_steps, V0r_step_real, R_V0r, weights, para_coeff(1), para_coeff(2), para_coeff(3))
        print*, "Getting R^2 = ", R2


        if (R2 > tol_R) then
            ! Statisfied
            best_V0r = fit_x_min*E_step
            best_R = fit_y_min
            best_V0r_step = V0r_step(best_id)
            print*, ""
            print*, "Smart search success"
            print*, "Energy overlapp: ", sum(N_overlapping_points_V0r(:,best_id) -1)*E_step
            RETURN
        
        else if (sum(weights) < 2*fit_range + 1) then
            ! add another point in the range
            do i = 1, fit_range
                if (.not. evaluated(closest_step - i)) then
                    next_step = closest_step - i
                    EXIT
                else if (.not. evaluated(closest_step + i)) then
                    next_step = closest_step + i
                    EXIT
                end if
            end do
            CYCLE

        else if (R2 > tol_R_2) then
            ! try to limit fit_range and see if it gets better

            if (fit_range > min_fit_range) then
                fit_range = fit_range - 1
            else
                ! Minimum present, but not well behaved. Return minimum value on grid instead of prediciton.
                ierr = 856
                best_V0r = V0r_step_real(best_id)*E_step
                best_V0r_step = V0r_step(best_id)
                best_R = R_V0r(best_id)
                RETURN
            end if

        else
            !Default to brute force
            ierr = 855
            fast_search = .False.
            exit
        end if


    end do

    ! *****************************************************************************************
    ! Brute force grid - not efficient, but will always give a result

    do i = 1, n_steps
        if (.not. evaluated(i)) then
            call r_pendry_beamset_y(n_E, E_step, n_beams, y1, y2, id_start_y1, id_start_y2, n_y1, n_y2, &
                V0r_step(i), & 
                R_V0r(i), &
                r_pendry_beams_V0r(:, i), &
                numerators_V0r(:, i), denominators_V0r(:, i), &
                N_overlapping_points_V0r(:, i), &
                ierr_tmp)
            evaluated(i) = .True.
            n_evaluated = n_evaluated + 1
            call update_best_V0r(n_steps, n_beams, best_R, best_id, R_V0r(i), i, &
                N_overlapping_points_V0r, N_overlapping_points, &
                numerators_V0r, numerators, &
                denominators_V0r, denominators, &
                r_pendry_beams_V0r, r_pendry_beams)
            !print*, "Brute force. step:", V0r_step(i), "R = ", R_V0r(i)
            if (ierr_tmp .ne. 0) ierr = ierr_tmp
            !print*, "V0r:", V0r_step(i)*E_step, "Overlapp:", sum(N_overlapping_points_V0r(:,i))*E_step
        else
            CYCLE
        end if
    end do

    ! If brute forced, report best V0r as the best step
    best_V0r_step = V0r_step(best_id)
    best_V0r = V0r_step_real(best_id)*E_step
    best_R = R_V0r(best_id)
    RETURN

end subroutine r_pendry_beamset_V0r_opt_on_grid

subroutine update_best_V0r(n_steps, n_beams, curr_best_R, curr_best_id, R, id, &
                N_overlapping_points_V0r, N_overlapping_points, &
                numerators_V0r, numerators, &
                denominators_V0r, denominators, &
                r_pendry_beams_V0r, r_pendry_beams)
    integer :: n_steps
    integer :: n_beams
    real(8), INTENT(INOUT) :: curr_best_R
    integer, INTENT(INOUT) :: curr_best_id

    real(8), INTENT(IN) :: R
    integer, INTENT(IN) :: id

    integer, INTENT(IN) :: N_overlapping_points_V0r(n_beams, n_steps)
    real(8), INTENT(IN) :: numerators_V0r(n_beams, n_steps)
    real(8), INTENT(IN) :: denominators_V0r(n_beams, n_steps)
    real(8), INTENT(IN) :: r_pendry_beams_V0r(n_beams, n_steps)

    integer, INTENT(OUT) :: N_overlapping_points(n_beams)
    real(8), INTENT(OUT) :: numerators(n_beams)
    real(8), INTENT(OUT) :: denominators(n_beams)
    real(8), INTENT(OUT) :: r_pendry_beams(n_beams)


    if (R < curr_best_R) then
        curr_best_R = R
        curr_best_id = id
        N_overlapping_points = N_overlapping_points_V0r(:, id)
        numerators = numerators_V0r(:, id)
        denominators = denominators_V0r(:, id)
        r_pendry_beams = r_pendry_beams_V0r(:, id)
    end if
    RETURN
end subroutine update_best_V0r

subroutine parabola_lsq_fit(n,x,y, w, B, ierr)
    integer, INTENT(IN) :: n 
    real(8), INTENT(IN) :: x(n), y(n), w(n)

    real(8), INTENT(OUT) :: B(3,1)
    integer,  INTENT(OUT):: ierr

    real(8) :: s_x, s_x2, s_y, s_x3, s_x4, s_xy, s_x2y
    real(8) :: A(3,3)
    integer :: LAPACK_info, i, ipiv(3)
    integer :: WORK(18)
    real(8) :: n_real

    n_real = n
    s_x     = 0
    s_x2    = 0
    s_y     = 0
    s_x3    = 0
    s_xy    = 0
    s_x4    = 0
    s_x2y   = 0

    A = 0
    B = 0

    do i = 1, n
        A(1,1) = A(1,1) + w(i)*x(i)**4
        A(1,2) = A(1,2) + w(i)*x(i)**3
        A(1,3) = A(1,3) + w(i)*x(i)**2
        A(2,3) = A(2,3) + w(i)*x(i)
        A(3,3) = A(3,3) + w(i)

        B(1,1) = B(1,1) + w(i)*y(i)*x(i)**2
        B(2,1) = B(2,1) + w(i)*y(i)*x(i)
        B(3,1) = B(3,1) + w(i)*y(i)

    end do
    ! Impose symmetry afterwards and save statements in loop
    A(2,1) = A(1,2)
    A(2,2) = A(1,3)
    A(3,1) = A(1,3)
    A(3,2) = A(2,3)

    ! print*, "A"
    ! print*, A(1,1:3)
    ! print*, A(2, 1:3)
    ! print*, A(3, 1:3)
    ! print*, "B"
    ! print*, B

    call DSYSV('U', 3, 1, A, 3, ipiv, b, 3, work, 18, LAPACK_info)
    if (LAPACK_info .ne. 0) ierr = 860




end subroutine parabola_lsq_fit

function parabola(x, a, b, c) result(y)
    real(8), INTENT(IN)  :: x
    real(8), INTENT(IN)  :: a, b, c
    real(8) :: y

    y = a*(x**2) + b*x + c
    RETURN
end function parabola

function parabola_R_squared(n, x, y, w, a, b, c) result(R_squared)
    integer, INTENT(IN) :: n
    real(8), INTENT(IN) :: x(n), y(n), w(n)
    real(8), INTENT(IN) :: a, b, c
    real(8) :: R_squared

    real(8) :: y_tmp, sum_weights, y_avg
    integer :: i
    real(8) :: SS_res, SS_tot

    y_avg = sum(y*w)/sum(w)

    R_squared = 0
    SS_res = 0
    SS_tot = 0

    do i = 1, n
        y_tmp = parabola(x(i), a, b, c)
        SS_res = SS_res + w(i)*((y(i) - y_tmp)**2)
        SS_tot = SS_tot + w(i)*((y(i) - y_avg)**2)
    end do
    
    R_squared = 1 - SS_res/SS_tot

end function parabola_R_squared


subroutine r_pendry_beamset_y(n_E, E_step, n_beams, y1, y2, id_start_y1, id_start_y2, n_y1, n_y2, V0r_shift, &
    r_pendry_weighted, r_pendry_beams, numerators, denominators, n_overlapping_points, ierr)
    !Rfactor_beamset:
    !INPUT: set of arrays Y1, Y2, Estart1, Estart2, V0rshift, beamtypes -> calls the above in a loop over beams
    !DOES: call above, average based on overlapping points, also split into types ("integer/fractional")
    !RETURNS: array RPe [overall and for each beamtype]
    implicit none

    !###############
    !INPUTS
    !###############
    !f2py integer, hidden, intent(in), depend(E_start1, E_start2), check(len(E_start1)==len(E_start2)) :: nr_beams=len(E_start1)
    integer :: n_E ! number of energy steps
    integer n_beams
    !f2py real(8), intent(in) :: E_step
    real(8) :: E_step
    !f2py real intent(in) :: y1 (n_E, n_beams), y2 (n_E, n_beams)
    real(8), intent(in) :: y1 (n_E, n_beams), y2 (n_E, n_beams)
    !f2py integer, intent(in) id_start_y1
    !f2py integer, intent(in) id_start_y2
    integer, intent(in)    :: id_start_y1(n_beams) ,id_start_y2(n_beams)
    !f2py integer, intent(in) n_y1
    !f2py integer, intent(in) n_y2
    integer, intent(in)    :: n_y1(n_beams) ,n_y2(n_beams)
    !f2py integer, intent(out) :: n_overlapping_points
    integer, intent(out) :: n_overlapping_points(n_beams)
    !f2py real(8), intent(out) :: r_pendry_beams(n_beams), r_pendry_weighted
    real(8), intent(out) :: r_pendry_beams(n_beams), r_pendry_weighted
    real(8), intent(out) :: numerators(n_beams), denominators(n_beams)
    !f2py integer, intent(in)::  V0r_shift
    integer, intent(in) :: V0r_shift

    integer, intent(out):: ierr

    ! temporay variables used in subroutine
    integer i
    real(8) temp_num, temp_denom

    ierr = 0

    ! iterate over beams and call subroutine r_factor_beam
    do i=1 , n_beams
        call r_pendry_beam_y(n_E, E_step, y1(:, i), y2(:, i), id_start_y1(i), id_start_y2(i), n_y1(i), n_y2(i), V0r_shift, &
            r_pendry_beams(i), numerators(i), denominators(i), n_overlapping_points(i))
    end do
    !total_points = sum(N_overlapping_points)

    if (ANY(ieee_is_nan(r_pendry_beams))) then
        r_pendry_weighted = ieee_value(real(8), ieee_signaling_nan)

        ! do i =1, n_beams
        !     if (ieee_is_nan(r_pendry_beams(i))) then
        !         print*,"NaN in beam ", i
        !     end if
        ! end do

        ierr = 811
        RETURN
    end if
    
    r_pendry_weighted = sum(numerators)/sum(denominators)

    return
end subroutine r_pendry_beamset_y

subroutine r2_beamset_intensity(n_E, E_step, n_beams, int_1, int_2, id_start_1, id_start_2, n_1, n_2, V0r_shift, &
    r2_weighted, r2_beams,  n_overlapping_points, ierr)

    !###############
    !INPUTS
    !###############
    !f2py integer, hidden, intent(in), depend(E_start1, E_start2), check(len(E_start1)==len(E_start2)) :: nr_beams=len(E_start1)
    integer :: n_E ! number of energy steps
    integer n_beams
    !f2py real(8), intent(in) :: E_step
    real(8) :: E_step
    real(8), intent(in) :: int_1 (n_E, n_beams), int_2 (n_E, n_beams)
    integer, intent(in)    :: id_start_1(n_beams) ,id_start_2(n_beams)
    integer, intent(in)    :: n_1(n_beams) ,n_2(n_beams)
    !f2py integer, intent(out) :: n_overlapping_points
    integer, intent(out) :: n_overlapping_points(n_beams)
    real(8), intent(out) :: r2_beams(n_beams), r2_weighted
    !f2py integer, intent(in)::  V0r_shift
    integer, intent(in) :: V0r_shift

    integer, intent(out):: ierr

    ! temporay variables used in subroutine
    integer i
    real(8) temp_num, temp_denom

    ierr = 0

    ! iterate over beams and call subroutine r_factor_beam
    do i=1 , n_beams
        call r2_beam_intensity(n_E, E_step, int_1(:, i), int_2(:, i), id_start_1(i), id_start_2(i), n_1(i), n_2(i), V0r_shift, &
            r2_beams(i), n_overlapping_points(i))
    end do

    if (ANY(ieee_is_nan(r2_beams))) then
        r2_weighted = ieee_value(real(8), ieee_signaling_nan)
        ierr = 811
        RETURN
    end if
    
    ! factro E_step in numerator and denominator cancels out
    ! we use (n_overlapping_points - 1) for the range since 1 overlapping point would mean energy interval of 0
    r2_weighted = sum(r2_beams*(n_overlapping_points - 1))/(sum(n_overlapping_points - 1))

    return
end subroutine r2_beamset_intensity


subroutine Rfactor_v0ropt(opt_type, min_steps, max_steps, nr_used_v0)
    !Rfactor_v0ropt:
    !INPUT: same as above, but with v0rmin/v0rmax, maybe v0rstep, optimization_type
    !DOES: V0r optimization with Rfactor_beamset, either as loop over grid, or parabola
    !RETURNS: array RPe [overall and for each beamtype], V0r_shift_optimum
    implicit none
    integer, intent(in) :: opt_type, min_steps, max_steps
    integer, intent(out) :: nr_used_v0

    !integer beamtypes(:)

    real(8) V0, R_Pe_V0
    ! Decide which optimization is used
    nr_used_v0 = 0 ! initialize to 0

    ! For grid search: V0 will be sampled on an equally spaced grid in range (V0_min, V0_max) with max_steps steps.
    ! For  example, with values V0_min = 2, V0_max = 5 and max_steps = 10, the sampled points would be:
    ! 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0
    ! The optional paramter min_steps is unused in this case.
    !
    ! For parabola search: V0 will initially be sampled in three points (V0_min, V0_max and (V0_min + V0_max)/2),
    ! then an attempt is made to fit a parbola of the form R(V0) = a*V0**2 + b*V0 + c. If the curvature of the parabola
    ! is negative, subsequent V0s are choosen close to the predicted minimum and such as to test the fit. The
    ! optimization will stop if the resultant R-factor has converged (criterion V0_conv) though at least
    ! min_steps values of V0 are evaluated. If the curvature of the parabola is positive, or the fit of the parabola
    ! becomes continually worse the optimization defaults to sampling a grid. (Similar as above, though already
    ! evaluated values will not be discarded).
    !
    if (opt_type == 0) then
    ! optimize grid
    elseif (opt_type == 1) then
    ! optimize with parabola, uses function parabolic_optimize
    end if

    ! Calculate next V0

    ! Calculate R factor, check if converged

    !call Rfactor_beamset(y1, y2, sizes_y1, sizes_y2, E_start1, E_start2, nr_beams, E_step, V0, &
    !                     beamtypes, R_Pe_V0, R_Pe_beams, N_overlapping_points)

    RETURN
end subroutine Rfactor_v0ropt


subroutine prepare_beams(n_beams, n_E_in, E_grid_in, intensities_in, E_start_beams, n_E_beams, &
                         skip_stages, n_beams_out, averaging_scheme, &
                         which_R, &
                         deg, &
                         n_E_out, E_grid_out, &
                         V0i, &
                         E_start_beams_out, n_E_beams_out, &
                         intpol_intensity, y_func, ierr)
    !Prepare_beams:
    !INPUT: array[I, E_min, E_step, NE], E_grid_step, averaging_scheme, smoothing?, E_min, E_max
    !DOES: (0) Limit_range, (1) Average/discard/reorder according to scheme; (2) smooth?; (3) interpolate on grid; (4) compute Y on new grid
    !Probably call several functions... see existing PREEXP, similar!
    !RETURNS: array of Y

    ! should this thing even return the Y function values ??

    !###############
    !INPUTS
    !###############
    ! Sizes
    integer, intent(in)                          :: n_beams          ! number of beams
    !f2py integer, intent(in) :: n_beams
    integer, INTENT(IN)                          :: n_E_in           ! number of energy values in the input
    !f2py integer, intent(hide), depend(E_grid_in) :: n_E_in=shape(E_grid_in,0)
    integer, INTENT(IN)                          :: n_E_out          ! number of energy values in the input
    !f2py integer, intent(hide), depend(E_grid_out) :: n_E_out=shape(E_grid_out,0)

    ! Operation information
    !f2py integer, intent(in), dimension(5) :: skip_stages
    integer, intent(in) :: skip_stages(5) ! which stages to execute
    !f2py integer, intent(in) :: n_beams_out
    integer, INTENT(IN) :: n_beams_out ! How many non-equivalent beams are there?
    !f2py integer, intent(in), dimension(n_beams) :: averaging_scheme
    integer, intent(inout) :: averaging_scheme(n_beams) ! How to average the beams
                                                     ! Integer between 1 and n_averaging_types - determines with which beams to average together
    !f2py intent(in) :: deg
    integer, INTENT(IN) :: deg ! Degree of interpolation to use
    integer, INTENT(IN) :: which_R ! Which R factor to prepare for (1: pendry, 2: R2)

    ! Energy grids
    real(8), INTENT(IN)                         :: E_grid_in(n_E_in) ! Input energy grid
    real(8), INTENT(IN)                         :: E_grid_out(n_E_out) ! Output energy grid

    ! Intensities is a rectangular array with values for each beam at each possible energy, though values will unused if the beam does not have data over the entire
    ! energy range. The array E_start_beams determines the index of the first used energy, the array n_E_beams determines the number of the number of
    ! energy steps to used by the beam. Gaps in the beam are not allowed here, but they can be implemented externally using the averaging functionalities.
    ! Input Data
    integer, intent(in)                          :: E_start_beams(n_beams) ! First energy step to be used for the beam
    integer, intent(in)                          :: n_E_beams(n_beams) ! Number of energy steps to use for the beam
    real(8), intent(in)                         :: intensities_in(n_E_in, n_beams) ! Beam intensities - input data
    
    real(8), INTENT(IN) :: V0i


    real(8), INTENT(OUT) :: intpol_intensity(n_E_out, n_beams)
    real(8) :: intpol_derivative(n_E_out, n_beams)
    real(8), INTENT(OUT) :: y_func(n_E_out, n_beams)


    !f2py integer, intent(out), dimension(n_beams) :: E_start_beams_out
    integer, INTENT(OUT) :: E_start_beams_out(n_beams) ! First energy step for each beam of the interpolated data
    !f2py integer, intent(out), dimension(n_beams):: n_E_beams_out
    integer, INTENT(OUT) :: n_E_beams_out(n_beams) ! Number of energy steps for each beam of the interpolated data

    integer, INTENT(OUT) :: ierr

    !###############
    ! Internal
    !###############
    integer :: new_min_index, new_max_index, cut_n_E_in
    integer                          :: cut_E_start_beams(n_beams) ! First energy step to be used for the beam
    integer                          :: cut_n_E_beams(n_beams) ! Number of energy steps to use for the beam

    integer :: beams_max_id_in(n_beams)
    integer :: cut_beams_max_id_in(n_beams)
    integer :: beams_max_id_out(n_beams_out)

    integer :: max_n_knots, max_nt, max_LHS_rows
    integer :: n_knots_beams(n_beams), nt_beams(n_beams), LHS_rows_beams(n_beams)
    real(8), ALLOCATABLE :: knots(:,:), LHS(:,:,:), coeffs(:,:)

    ! TYPE(grid_pre_evaluation) :: grid_pre_eval(n_beams)
    ! type(deriv_pre_evaluation):: pre_eval_derivs(n_beams)
    ! type(grid)                :: grid_origin(n_beams)
    ! type(grid)                :: grid_target(n_beams)

    real(8)    :: intensities(n_E_in, n_beams)


    integer :: ierrs(n_beams)


    integer :: i ! Loop indices

    ierr = 0
    ierrs = 0
    intensities = intensities_in

    beams_max_id_in = E_start_beams + n_E_beams -1

    !###############
    ! Limit range of beams
    !###############
    ! Assumes, experimental beams are recorded on the same energy grid already. This should always be the case.
    ! Thus only the information of E_min and n_E are required for each beam
    
    ! discard all data outside of final grid (except one data point on either side for interpolation, if available)
    ! This requires:
    ! 1) cutting the E_grid_in
    ! 2) recalculating the starting and end indices of all beams on the new grid

    if (skip_stages(1).EQ.0) then
        !print*, E_grid_in
        !print*, E_grid_out
        ! minimum index to keep:
        if (E_grid_in(1) .le. E_grid_out(1)) then
            new_min_index = find_grid_correspondence(E_grid_out(1), n_E_in, E_grid_in, 1)
        else
            new_min_index = 1
        end if

        if (E_grid_in(n_E_in) .ge. E_grid_out(n_E_out)) then
            new_max_index = find_grid_correspondence(E_grid_out(n_E_out), n_E_in, E_grid_in, new_min_index)
        else
            new_max_index = n_E_in
        end if
        

        cut_n_E_in = new_max_index - new_min_index + 1
        
        do concurrent( i=1: n_beams)
            ! new start and end indices
            cut_E_start_beams(i) = max(E_start_beams(i), new_min_index)
            cut_beams_max_id_in(i) = min(new_max_index, beams_max_id_in(i))
            cut_n_E_beams(i) = cut_beams_max_id_in(i) - cut_E_start_beams(i) + 1
            if (cut_n_E_beams(i) < 2*deg+1) then
                ! Beam cannot be used ...skip somehow...
                ierrs(i) = 12345 !TODO give error code... at least one beam did not contain usable information... - return averaging_scheme with corresponding beams?
                ! Change energy range or discard beam 
                averaging_scheme(i) = 0 ! (can I implement this via averaging group 0?)
            end if
            !print*, "cut_n_E_in", cut_n_E_in
            !print*, "cut_E_start_beams", cut_E_start_beams(i)
            
        end do
        ierr = MAXVAL(ierrs)
    else
        new_min_index = 1
        new_max_index = n_E_in
        cut_E_start_beams = E_start_beams
        cut_n_E_beams = n_E_beams
        cut_n_E_in = n_E_in
        cut_beams_max_id_in = beams_max_id_in

    end if

    !print*, "new_min_index, new_max_index", new_min_index, new_max_index
    !print*, "E_start_beams", cut_E_start_beams(15:20)
    !print*, "n_E_out", n_E_out(15:20)
    !print*, "max id beams", cut_beams_max_id_in(15:20)

    
    !###############
    ! Average/discard/reorder
    !###############

    ! TODO implement/ fix averaging
    if (skip_stages(2).EQ.0) then
        if (n_beams_out > n_beams) then
            ierr = 220
            RETURN
        end if
        call avg_reo_disc(n_beams, n_beams_out, n_E_in, averaging_scheme, 2*deg+1, &
        intensities, cut_E_start_beams, cut_n_E_beams, &
        ierr)

        if (ierr .ne. 0) then
            print*, "Oh no!", ierr
            RETURN
        end if
    else
        if (n_beams_out .ne. n_beams) then
            ierr = 223
            RETURN
        end if
    end if
    ! print*, "E_start_beams:"
    ! print*, cut_E_start_beams
    ! print*, "n_E_beams:"
    ! print*, cut_n_E_beams

    E_start_beams_out(:) = 0
    n_E_beams_out(:) = 0 

    do concurrent( i=1: n_beams_out)
        ! Set out beam indices & length
        E_start_beams_out(i) = find_grid_correspondence( &
                E_grid_in(cut_E_start_beams(i)), &
                n_E_out, &
                E_grid_out, &
                1 &
            )

        n_E_beams_out(i) = - E_start_beams_out(i) + &
            find_grid_correspondence( &
                E_grid_in(cut_beams_max_id_in(i)), &
                n_E_out, &
                E_grid_out, &
                E_start_beams_out(i) &
            ) + 1
        beams_max_id_out(i) = E_start_beams_out(i) + n_E_beams_out(i) - 1
    end do



    ! Debug below...
    do i = 1, n_beams
       ! print*, i
       ! print*, n_E_out
       ! print*, beams_max_id_out(i)
        if (E_grid_in(cut_beams_max_id_in(i)) < E_grid_out(beams_max_id_out(i))) then
            print*, "issue size out - beam", i
        end if
        if (E_grid_in(cut_E_start_beams(i))>E_grid_out(E_start_beams_out(i))) then
            print*, "issue size in - beam", i
        end if
    end do


    !###############
    ! Smoothing
    !###############

    ! Smoothing is not (yet) implemented.

    ! This is a tricky issue, because excessive smoothing strongly affects the R-factor, so you want to avoid over-smoothing at all cost.
    ! We expect the experimental data to come in from the beam-extraction already pre-smoothed. If we implement an additional later smoothing here,
    ! it only makes sense to do it equally on experimental and theoretical data. Thus, both should use the same smoothing kernel and probably not 
    ! use a smoothing parameter larger than ~6 eV.y

    if (skip_stages(3).EQ.0) then
        ierr = 23
        continue
    end if


    !###############
    ! Interpolation on new grid
    !###############

    ! How to deal with coeffs array? - size depends on deg:
    ! size of coeffs is nt= n_knots - deg -1 = n + 2*deg -deg -1 = n + deg -1 
    ! where max(n) is the number of datapoints in, i.e. n_E_in. Thus max(nt) = n_E_in + deg + 1.
    ! Allocating coeffs with size (n_E_in + deg) should therfore cover all possible beam sizes. Unfortunately, we need to keep track of the number of coeffs per beam in an array though.
    if (skip_stages(4).EQ.0) then

        intpol_derivative  = 0
        intpol_intensity = 0
        y_func = 0

        call perform_checks(cut_n_E_in, E_grid_in(new_min_index:new_max_index), n_E_out, E_grid_out, ierr)
        if (ierr .ne. 0) RETURN

        call get_array_sizes(cut_n_E_in, deg, max_n_knots, max_nt, max_LHS_rows)
        Allocate(knots(max_n_knots+1, n_beams), LHS(max_LHS_rows+1, max_nt+1, n_beams), coeffs(max_nt+1, n_beams))
        
        if (which_R == 1) then
            do i= 1,n_beams_out
                ! print*, "beam: ", i
                ! !DEBUG
                ! print*, cut_E_start_beams(i), cut_E_start_beams(i)+cut_n_E_beams(i) - 1
                ! print*, E_start_beams_out(i), E_start_beams_out(i)+n_E_beams_out(i) - 1
                ! print*, intensities(cut_E_start_beams(i), i)
                ! print*, intensities( cut_beams_max_id_in(i) , i)
                ! print*, E_grid_in(cut_E_start_beams(i))
                ! print*, E_grid_in(cut_beams_max_id_in(i))
                
                call get_array_sizes(cut_n_E_beams(i), deg, n_knots_beams(i), nt_beams(i), LHS_rows_beams(i))
                ! print*, cut_E_start_beams(i), cut_beams_max_id_in(i)
                ! print*, ""
                ! print*, n_knots_beams(i), nt_beams(i)


                call single_calc_spline(&
                    cut_n_E_beams(i), &
                    E_grid_in(cut_E_start_beams(i) : cut_beams_max_id_in(i)), &
                    intensities(cut_E_start_beams(i) : cut_beams_max_id_in(i), i), &
                    deg, &
                    n_knots_beams(i), nt_beams(i), LHS_rows_beams(i), &
                    knots(1:n_knots_beams(i), i), &
                    coeffs(1:nt_beams(i), i), &
                    ierrs(i))
                ! !Interpolate intensitites
                !     print*, "ierr", ierr
                !     print*, "First and last intensity in"
                !     print*, intensities_in(cut_E_start_beams(i),i)
                !     print*, intensities_in(cut_E_start_beams(i) + cut_n_E_beams(i) - 1, i)
                !     print*, "Any NAN?", ANY(ieee_is_nan(intensities_in))
                !     print*, "First and last knot"
                !     print*, knots(1,i)
                !     print*, knots(n_knots_beams(i), i)
                !     print*, "n_knots, nt, deg, LHS_rows", n_knots_beams(i), nt_beams(i), deg, LHS_rows_beams(i)

                call single_interpolate_coeffs_to_grid( &
                    n_knots_beams(i), knots(1:n_knots_beams(i), i), nt_beams(i), coeffs(1:nt_beams(i), i), &
                    deg, &
                    n_E_beams_out(i), &
                    E_grid_out(E_start_beams_out(i) : beams_max_id_out(i)), &
                    0, &
                    intpol_intensity(E_start_beams_out(i) : beams_max_id_out(i), i))
                !print*, "Intpol"
                !print*, intpol_intensity(E_start_beams_out(i) : E_start_beams_out(i) + n_E_beams_out(i),i)
                    ! Interpolate derivatives

                ! Intensity must be > 0, so set to 0 if not:
                intpol_intensity = max(intpol_intensity, 0d0)

                call single_interpolate_coeffs_to_grid( &
                    n_knots_beams(i), knots(1:n_knots_beams(i), i), nt_beams(i), coeffs(1:nt_beams(i), i), &
                    deg, &
                    n_E_beams_out(i), &
                    E_grid_out(E_start_beams_out(i) : beams_max_id_out(i)), &
                    1, &
                    intpol_derivative(E_start_beams_out(i) : beams_max_id_out(i), i))

                ! Calculate pendry Y function
                call pendry_y( &
                    n_E_beams_out(i), &
                    intpol_intensity(E_start_beams_out(i):beams_max_id_out(i), i), &
                    intpol_derivative(E_start_beams_out(i):beams_max_id_out(i), i), &
                    v0i, &
                    y_func(E_start_beams_out(i):beams_max_id_out(i), i) &
                    )
            end do
            ierr = MAXVAL(ierrs)
            !end Pendry
        else if (which_R == 2) then
            ! Beginn R2
            do i= 1,n_beams_out
                
                call get_array_sizes(cut_n_E_beams(i), deg, n_knots_beams(i), nt_beams(i), LHS_rows_beams(i))

                call single_calc_spline(&
                    cut_n_E_beams(i), &
                    E_grid_in(cut_E_start_beams(i) : cut_beams_max_id_in(i)), &
                    intensities(cut_E_start_beams(i) : cut_beams_max_id_in(i), i), &
                    deg, &
                    n_knots_beams(i), nt_beams(i), LHS_rows_beams(i), &
                    knots(1:n_knots_beams(i), i), &
                    coeffs(1:nt_beams(i), i), &
                    ierrs(i))

                call single_interpolate_coeffs_to_grid( &
                    n_knots_beams(i), knots(1:n_knots_beams(i), i), nt_beams(i), coeffs(1:nt_beams(i), i), &
                    deg, &
                    n_E_beams_out(i), &
                    E_grid_out(E_start_beams_out(i) : beams_max_id_out(i)), &
                    0, &
                    intpol_intensity(E_start_beams_out(i) : beams_max_id_out(i), i))
                
                ! Done! R2 only needs intensities
                y_func = ieee_value(real(8), ieee_signaling_nan) ! Assign NaNs so it throws if used by accident
            end do
        else
            ierr = 701
        end if
    end if

    !###############
    ! Y function calculation on new grid
    !###############

    ! TODO: clean up implementation

    ! output is the y functions

    ! i = 4

    ! print* , n_E_beams(i), cut_n_E_beams(i), n_E_beams_out(i)

    ! print*, "Intensities old index l:", intensities_in(E_start_beams(i), i), &
    ! intensities(E_start_beams(i) +1, i)
    ! print*, "Intensities old index u:", intensities_in(E_start_beams(i) + n_E_beams(i) - 2, i), &
    ! intensities(E_start_beams(i) + n_E_beams(i) -1, i)

    ! print*, "Intensities new id:", intensities(cut_E_start_beams(i), i), &
    ! intensities(cut_E_start_beams(i) +1, i)
    ! print*, "Intensities new id:", intensities(cut_beams_max_id_in(i) -1, i), &
    ! intensities(cut_beams_max_id_in(i), i)

    ! print*, "Intensities out l:", intpol_intensity(E_start_beams_out(i), i), &
    ! intpol_intensity(E_start_beams_out(i) +1, i)
    ! print*, "Intensities out u:", intpol_intensity(beams_max_id_out(i) -1, i), &
    ! intpol_intensity(beams_max_id_out(i), i)


    RETURN

end subroutine prepare_beams

pure subroutine find_index_on_new_grid(n_E_in, E_in, n_E_out, E_out, id_start_in, n_steps_in,&
                                  id_start_out, n_steps_out)
    integer, INTENT(IN) :: n_E_in, n_E_out
    real(8), INTENT(IN) :: E_in(n_E_in), E_out(n_E_out)
    integer, INTENT(IN) :: id_start_in, n_steps_in

    integer, INTENT(OUT) :: id_start_out, n_steps_out

    integer :: id_end_in, id_end_out

    id_end_in = id_start_in + n_steps_in
    id_start_out = min(1, find_interval(E_out, n_E_out, 0, E_in(id_start_in), 1, 0))
    id_end_out = max(n_E_out, find_interval(E_out, n_E_out, 0, E_in(id_start_out), 1, 0))
    n_steps_out = id_end_out - id_start_out
    RETURN
end subroutine find_index_on_new_grid

pure integer function find_grid_correspondence(value, n_steps, steps, min_bound) result(interval)
    ! returns the index of the step lower than value
    real(8), INTENT(IN) :: value
    integer, INTENT(IN) :: n_steps
    real(8), INTENT(IN) :: steps(n_steps)
    integer, INTENT(IN) :: min_bound

    integer :: i

    interval = min_bound

    do i = min_bound, n_steps
        if (steps(i) < value) then
            interval = interval + 1
        else
            RETURN
        end if
    end do
    RETURN
end function find_grid_correspondence


subroutine avg_reo_disc(n_beams_in, n_beams_out, n_E, scheme, min_len, &
    intensities, beam_start, beam_n, &
    ierr)
    ! Average, reorder and discard beams according to scheme
    ! scheme is an array of integers between 0 and n_beams_out
    ! beams are averaged into the 
    ! For any beam with scheme == 0, the beam is discarded.
    integer, INTENT(IN) :: n_beams_in
    integer, INTENT(IN) :: n_beams_out
    integer, INTENT(IN) :: n_E
    integer, INTENT(IN) :: min_len
    real(8), INTENT(INOUT) :: intensities(n_E, n_beams_in)
    integer, INTENT(INOUT) :: scheme(n_beams_in)
    integer, INTENT(INOUT) :: beam_start(n_beams_in)
    integer, INTENT(INOUT) :: beam_n(n_beams_in)


    real(8) :: intensities_in(n_E,  n_beams_in)
    integer :: beam_start_in(n_beams_in)
    integer:: beam_n_in(n_beams_in)
    integer, INTENT(OUT) :: ierr

    integer :: i,j
    integer :: avg_counter(n_beams_out)

    intensities_in = intensities
    beam_start_in = beam_start
    beam_n_in = beam_n

    beam_start(:) = 0
    beam_n(:) = n_E

    ! Adjust beam length if necessary...
    do i = 1, n_beams_out
        do j = 1, n_beams_in
            if(scheme(j)==i) then
                if (beam_start_in(j) > beam_start(i)) then
                    beam_start(i) = beam_start_in(j)
                end if 
                if (beam_start(i) + beam_n(i) > beam_start_in(j) + beam_n_in(j)) then
                    beam_n(i) = beam_start_in(j) + beam_n_in(j) - beam_start(i)
                end if
            end if
        end do
    end do

    ! Decreasing the beam lengths may have made some beams unviable for calculations...
    do j = 1, n_beams_in
        if (beam_n(j) < min_len) then
            scheme(j) = 0
            ierr = 211
        end if
    end do

    intensities(:,:) = 0.d0

    avg_counter(:) = 0
    ! Average per scheme
    do i = 1, n_beams_out
        do j = 1, n_beams_in
            if(scheme(j)==i) then
                intensities(:,i) = intensities(:,i) + intensities_in(:,j)
                avg_counter(i) = avg_counter(i) + 1
            end if
        end do
        ! Divide by nr beams if more than one
        if (avg_counter(i) > 0) then
            intensities(:,i) = intensities(:,i)/avg_counter(i)
        end if
    end do
    RETURN
end subroutine avg_reo_disc

subroutine range_index_from_Energy(E_min_current, NE_in, E_step, E_min_cut, E_max_cut, new_start_stop_step)
    integer, intent(in) :: NE_in
    real(8), intent(in) :: E_min_current, E_min_cut, E_max_cut, E_step

    integer , intent(out) :: new_start_stop_step(3)

    real(8) :: old_max, new_max, new_min

    old_max = E_min_current + E_step*NE_in
    new_max = min(E_max_cut, old_max)
    new_min = max(E_min_current, E_min_cut)

    new_start_stop_step(3) = floor((E_max_cut-E_min_cut)/E_step) ! typecast to integer; will cut off float...

    new_start_stop_step(1) = (E_min_cut - E_min_current)/E_step
    new_start_stop_step(2) = new_start_stop_step(1) + new_start_stop_step(3)
    return
end subroutine range_index_from_Energy




! Routines for working with rfactors

subroutine r_pendry_beamtype_grouping(n_beams, numerators, denominators, n_overlapping_points, &
        n_groups, grouping, r_factor_groups, n_overlapping_points_groups, ierr)

    integer, INTENT(IN) :: n_beams
    real(8), INTENT(IN) :: numerators(n_beams)
    real(8), INTENT(IN) :: denominators (n_beams)
    integer, INTENT(IN) :: n_overlapping_points(n_beams)

    integer, INTENT(IN) :: n_groups
    integer, INTENT(IN) :: grouping(n_beams)

    real(8), INTENT(OUT) :: r_factor_groups(n_groups)
    integer, INTENT(OUT) :: ierr

    integer, INTENT(OUT) :: n_overlapping_points_groups(n_groups)

    integer :: group, i
    real(8) :: num(n_groups), denom(n_groups)

    ierr = 0

    ! initialize
    num = 0
    denom = 0
    n_overlapping_points_groups = 0
    
    if (ANY(grouping < 0) .or. ANY(grouping > n_groups)) then
        ierr = 902
        RETURN
    end if

    do i = 1, n_beams
        group = grouping(i)
        if (group == 0) CYCLE
        num(group) = num(group) + numerators(i)
        denom(group) = denom(group) + denominators(i)
        n_overlapping_points_groups(group) = n_overlapping_points_groups(group) + n_overlapping_points(i)
    end do

    do group = 1, n_groups
        if (abs(denom(group)) .ge. 1e-10) then
            r_factor_groups(group) = num(group)/denom(group)
        else
            r_factor_groups(group) = ieee_value(real(8), ieee_signaling_nan)
            ierr = 903
        end if 
    end do

end subroutine r_pendry_beamtype_grouping 





! Library functions    


subroutine pendry_y(n_data, intensity, derivative, v0i, y_func)
    integer, INTENT(IN) :: n_data
    real(8), INTENT(IN) :: intensity(n_data)
    real(8), INTENT(IN) :: derivative(n_data)
    real(8), INTENT(IN) :: v0i

    real(8), INTENT(OUT) :: y_func(n_data)

    y_func = intensity*derivative /(intensity*intensity + v0i*v0i*derivative*derivative)

end subroutine pendry_y


pure function trapez_integration_const_dx(f, dx) result(integral)
    implicit none
    real(8), intent(in)    :: f(:)
    real(8), intent(in)    :: dx
    real(8) integral
    integer n

    n = size(f)
    integral = 0
    integral = sum(f(1:n-1)+f(2:n))*dx/2
    return
end function trapez_integration_const_dx



! Compatibility functions
subroutine translate_evaluation_grid(e_min, e_max, EINCR, VINCR, V0RR, V01, V02, &
    n_E, energies, E_step)

    ! e_min != EMIN, e_max != EMAX
    real(4), INTENT(IN) :: e_min, e_max, EINCR, VINCR
    real(4), INTENT(IN) :: V0RR, V01, V02

    integer, INTENT(OUT) :: n_E
    real(8), INTENT(OUT), ALLOCATABLE :: energies(:)
    real(8), INTENT(OUT) :: E_step

    integer i
    real(8) :: min_tmp, max_tmp

    E_step = real(ABS(MIN(EINCR, VINCR)), 8)
    ! V01 and V02 need to be considered, because we want to keep extra data for V0r shifts
    min_tmp = e_min + real(-V0RR + V01, 8)
    max_tmp = e_max + real(-V0RR + V02, 8) + E_step

    n_E = int((max_tmp - min_tmp)/E_step)

    allocate(energies(n_E))

    !print*, size(energies)

    do i = 1, n_E
        energies(i) = min_tmp + (i-1)*E_step
    end do
    RETURN
end subroutine translate_evaluation_grid


subroutine translate_exp(n_beams, MNGP, NEE, EE, AE, &
    n_E, exp_grid, E_start_beams, n_E_beams, intensities)
    
    integer, INTENT(IN) :: n_beams
    integer, INTENT(IN) :: MNGP
    integer, INTENT(IN) :: NEE(n_beams)
    real(4), INTENT(IN) :: EE(n_beams, MNGP), AE(n_beams, MNGP)
    
    integer, INTENT(OUT) :: n_E
    integer, INTENT(OUT) :: n_E_beams(n_beams), E_start_beams(n_beams)
    real(8), INTENT(OUT), ALLOCATABLE :: exp_grid(:), intensities(:,:)
  
    
    real(8) :: min_energy, max_energy, E_step
    integer :: i, beam
    
    n_E_beams = NEE

    ! get E_step
    E_step = real(EE(1,2) - EE(1,1), 8)
    
    min_energy = MINVAL(EE(:,1))
    max_energy = MAXVAL(EE)
    print*, min_energy, max_energy, E_step

    n_E = NINT((max_energy - min_energy) / E_step) + 1
    allocate(exp_grid(n_E))

    do i = 1, n_E
      exp_grid(i) = min_energy + E_step*(i-1)
    end do
    
    allocate(intensities(n_E, n_beams))
    intensities = ieee_value(real(8), ieee_signaling_nan)
        ! row major order used in EE for some reason
        
    do beam = 1, n_beams
        do i = 1, n_E
            if (ABS(EE(beam,1) - exp_grid(i)) < 0.01) then
                E_start_beams(beam) = i
                exit
            end if
        end do
        
        intensities(E_start_beams(beam) : E_start_beams(beam) + n_E_beams(beam) -1, beam) = &
            AE(beam, 1 : NEE(beam))
    end do
    RETURN
end subroutine translate_exp


subroutine translate_theo(n_beams,  MNGP, MNBTD, ES, ATS, &
    n_E, theo_grid, E_start_beams, n_E_beams,  intensities)

    integer, INTENT(IN) :: n_beams
    integer, INTENT(IN) :: MNGP, MNBTD
    real(4), INTENT(IN) :: ES(MNGP), ATS(1, MNBTD, MNGP)

    integer, INTENT(OUT) :: n_E
    real(8), INTENT(OUT), ALLOCATABLE :: theo_grid(:)

    integer, INTENT(OUT) :: n_E_beams(n_beams), E_start_beams(n_beams)
    real(8), INTENT(OUT), ALLOCATABLE :: intensities(:,:)

    integer :: i, beam, en
    integer :: tmp_max_id_in_grid, tmp_max_id_out_grid


    ! Theory data is available for all energies - find which ones
    do i = 1, MNGP
        if ((ES(i) < 0.01)) then
            n_E = i - 1
            EXIT
        end if
    end do

    allocate(theo_grid(n_E))
    theo_grid = ES(1:n_E)
    
    ! is all data always available for 
    n_E_beams = n_E

    ALLOCATE(intensities(n_E, n_beams))
    intensities = ieee_value(real(8), ieee_signaling_nan)

    do beam = 1,MNBTD
        do en = 1,n_E
            if (ATS(1, beam, en) > 1e-8) then
                E_start_beams(beam) = en
                EXIT
            end if
        end do
        n_E_beams(beam) = n_E - E_start_beams(beam) + 1
        do en = E_start_beams(beam), n_E
            if (ATS(1, beam, en) < 1e-8) then
                ! max id reached
                n_E_beams = en - E_start_beams(beam) + 1
                EXIT
            end if
        end do

        intensities(E_start_beams(beam) : E_start_beams(beam) + n_E_beams(beam) - 1, beam) =  &
            ATS(1, beam, E_start_beams(beam) : E_start_beams(beam) + n_E_beams(beam) - 1)
    end do

    RETURN
end subroutine translate_theo

end module r_factor_new