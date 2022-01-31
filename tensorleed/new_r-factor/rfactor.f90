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

    integer :: id_start_y1_tmp, id_start_y2_tmp, id_min, id_max

    ! V0r_shift shifts y2 by some number of grid points -> id_start_y2 = id_start_y2 + V0r_shift
    id_start_y1_tmp = id_start_y1
    id_start_y2_tmp = id_start_y2 + V0r_shift

    id_min = max(id_start_y1_tmp, id_start_y2_tmp)
    id_max = min(id_start_y1_tmp + n_y1, id_start_y2_tmp+ n_y2)
    N_overlapping_points = id_max - id_min
    allocate(y_diff(N_overlapping_points), y_squared_sum(N_overlapping_points))

    ! Remember that indices are shifted by V0r
    y_diff=y1(id_min: id_max) - y2(id_min + V0r_shift: id_max + V0r_shift) ! difference between Y functions
    y_squared_sum=y1(id_min: id_max)**2+y2(id_min + V0r_shift: id_max +V0r_shift)**2

    ! caclulate numerator = integral (Y1-Y2)**2 dE
    numerator = trapez_integration_const_dx(y_diff**2, E_step)
    ! calculate denominator = integral (Y1**2+Y2**2) dE
    denominator = trapez_integration_const_dx(y_squared_sum, E_step)

    R_pendry = numerator/denominator
    return

end subroutine r_pendry_beam_y



!V0rshift not yet implemented
    ! beamtypes not yet implemented -> does it really make sense to do that here even?
    ! maybe extra subroutine actually
subroutine r_pendry_beamset_y(n_E, E_step, n_beams, y1, y2, id_start_y1, id_start_y2, n_y1, n_y2, V0r_shift, &
    r_pendry_weighted, r_pendry_beams, n_overlapping_points, ierr)
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
    !f2py integer, intent(in)::  V0r_shift
    integer, intent(in) :: V0r_shift

    integer, intent(out):: ierr

    real(8) numerators(n_beams), denominators(n_beams)

    ! temporay variables used in subroutine
    integer i
    integer total_points
    real(8) temp_num, temp_denom

    ierr = 0

    ! iterate over beams and call subroutine r_factor_beam
    do i=1 , n_beams
        call r_pendry_beam_y(n_E, E_step, y1(:, i), y2(:, i), id_start_y1(i), id_start_y2(i), n_y1(i), n_y2(i), V0r_shift, &
            r_pendry_beams(i), numerators(i), denominators(i), n_overlapping_points(i))
    end do
    total_points = sum(N_overlapping_points)

    if (ANY(ieee_is_nan(r_pendry_beams))) then
        r_pendry_weighted = ieee_value(real(8), ieee_signaling_nan)
        ierr = 811
    end if
    r_pendry_weighted = sum(numerators/denominators*N_overlapping_points)/total_points

    return
end subroutine r_pendry_beamset_y

subroutine Rfactor_beamtypes(y1, sizes_y1, y2, sizes_y2, E_start1, E_start2, nr_beams, E_step, V0rshift, &
                           beamtypes, nr_beamtypes, R_Pe_weighted, R_Pe_beams, N_overlapping_points)
    !###############
    !INPUTS
    !###############
    !f2py integer, hidden, intent(in), depend(E_start1, E_start2), check(len(E_start1)==len(E_start2)) :: nr_beams=len(E_start1)
    integer nr_beams
    !f2py real(8), intent(in) :: E_step
    real(8) E_step
    !f2py integer, intent(in) :: sizes_y1
    integer, intent(in) :: sizes_y1(nr_beams)
    !f2py integer, intent(in) :: sizes_y2
    integer, intent(in) :: sizes_y2(nr_beams)
    !f2py real intent(in) :: y1(:,:), y2(:,:)
    real(8), intent(in) :: y1 (:,:), y2 (:,:)
    !f2py intent(in) E_start_1
    !f2py intent(in) E_start_2
    real(8), intent(in)    :: E_start1(nr_beams), E_start2(nr_beams)
    integer, intent(in) :: beamtypes(nr_beams)
    integer, intent(in) :: nr_beamtypes
    !f2py integer, intent(out) :: N_overlapping_points(:)
    integer, intent(out) :: N_overlapping_points(nr_beams)
    !f2py real(8), intent(out) :: R_Pe_beams(nr_beams), R_Pe_weighted(nr_beamtypes)
    real(8), intent(out) :: R_Pe_beams(nr_beams), R_Pe_weighted(nr_beamtypes)
    !f2py real(8), intent(in)::  V0rshift
    real(8), intent(in) :: V0rshift

    ! variables used internally
    integer beam, group, points_group, tmp_points(nr_beams)
    !
    real(8) tmp_num(nr_beams), tmp_denom(nr_beams), tmp_pendry


    ! beamtypes is array that contains assignment o beamtype group for each beam
    ! contains integers from 1 to nr_beamtypes

    ! Loop over beamtype groups
    do group=1, nr_beamtypes
        tmp_points = 0
        tmp_num = 0d0
        tmp_denom = 0d0
        do beam = 1, nr_beams
            if (beamtypes(beam)==group) then
                !call r_factor_beam(y1(beam,:), sizes_y1(beam), y2(beam,:), sizes_y2(beam), E_start1(beam), &
                !        E_start2(beam), E_step, V0rshift, tmp_pendry, tmp_num(beam), tmp_denom(beam), tmp_points(beam))
            end if
        end do
        points_group = sum(tmp_points)
        R_Pe_weighted(group) = sum(tmp_num/tmp_denom*points_group)/points_group

    end do

    return
end subroutine Rfactor_beamtypes

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
                         deg, &
                         n_E_out, E_grid_out, &
                         E_start_beams_out, n_E_beams_out, &
                         intpol_intensity, intpol_derivative, y_func, ierr)
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

    integer :: max_n_knots, max_nt, max_LHS_rows
    integer :: n_knots_beams(n_beams), nt_beams(n_beams), LHS_rows_beams(n_beams)
    real(8), ALLOCATABLE :: knots(:,:), LHS(:,:,:), coeffs(:,:)

    ! TYPE(grid_pre_evaluation) :: grid_pre_eval(n_beams)
    ! type(deriv_pre_evaluation):: pre_eval_derivs(n_beams)
    ! type(grid)                :: grid_origin(n_beams)
    ! type(grid)                :: grid_target(n_beams)

    real(8)    :: intensities(n_E_in, n_beams)

    real(8), INTENT(OUT) :: intpol_intensity(n_E_out, n_beams)
    real(8), INTENT(OUT) :: intpol_derivative(n_E_out, n_beams)
    real(8), INTENT(OUT) :: y_func(n_E_out, n_beams)

    integer :: ierrs(n_beams)
    real(8) :: v0i

    integer :: i ! Loop indices

    ierr = 0
    ierrs = 0
    intensities = intensities_in

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
            print* , "This should not happen"
            ierrs(:) = 211
            RETURN
        end if

        if (E_grid_in(n_E_in) .ge. E_grid_out(n_E_out)) then
            new_max_index = find_grid_correspondence(E_grid_out(n_E_out), n_E_in, E_grid_in, new_min_index)
        else
            print* , "This should not happen"
            ierrs(:) = 211
            RETURN
            !id_end = n_E_in
        end if
         
        cut_n_E_in =  new_max_index - new_min_index +1
        
        do concurrent( i=1: n_beams)
            ! new start and end indices
            !cut_E_start_beams(i) = max(E_start_beams(i)-new_min_index+1, new_min_index)
            cut_E_start_beams(i) = max(E_start_beams(i), new_min_index)
            cut_n_E_beams(i) = min(new_max_index, E_start_beams(i)+n_E_beams(i)) - cut_E_start_beams(i)
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

    end if

    ! print*, "new_min_index, new_max_index", new_min_index, new_max_index
    ! print*, "E_start_beams", cut_E_start_beams(15:20)
    ! print*, "n_E_out", n_E_out
    


    if (ANY(ieee_is_nan(intensities))) print*,"NaN before average"

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
            E_grid_in(cut_E_start_beams(i) + cut_n_E_beams(i) - 1), &
            n_E_out, &
            E_grid_out, &
            E_start_beams_out(i) &
            )
    end do

    ! Debug below...
    do i = 1, n_beams
        if (E_grid_in(cut_E_start_beams(i) + cut_n_E_beams(i)-1) < E_grid_out(E_start_beams_out(i) + n_E_beams_out(i)-1)) then
            print*, "issue size out - beam", i
        end if
        if (E_grid_in(cut_E_start_beams(i))>E_grid_out(E_start_beams_out(i))) then
            print*, "issue size in - beam", i
        end if
    end do

    v0i = 2.0d0


    !###############
    ! Smoothing
    !###############

    if (skip_stages(3).EQ.0) then
        ierr = 23
        write(*,*) "Smoothing not yet implemented!" ! TODO: implement smoothing in prepare_beams
    end if

    !if (ANY(ieee_is_nan(intensities))) print*,"NaN before interpolation"

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
        Allocate(knots(max_n_knots, n_beams), LHS(max_LHS_rows, max_nt, n_beams), coeffs(max_nt, n_beams))
        

        do i= 1,n_beams_out
            !DEBUG
            !print*, cut_E_start_beams(i), cut_E_start_beams(i)+cut_n_E_beams(i)
            !print*, E_start_beams_out(i), E_start_beams_out(i)+n_E_beams_out(i)
            !print*, intensities(cut_E_start_beams(i), i)
            !print*, intensities(cut_E_start_beams(i)+cut_n_E_beams(i), i)
            
            call get_array_sizes(cut_n_E_beams(i), deg, n_knots_beams(i), nt_beams(i), LHS_rows_beams(i))
            call single_calc_spline(&
                cut_n_E_beams(i), &
                E_grid_in(cut_E_start_beams(i) : cut_E_start_beams(i) + cut_n_E_beams(i)), &
                intensities_in(cut_E_start_beams(i) : cut_E_start_beams(i) + cut_n_E_beams(i), i), &
                deg, &
                n_knots_beams(i), nt_beams(i), LHS_rows_beams(i), &
                knots(1:n_knots_beams(i), i), &
                coeffs(1:nt_beams(i), i), &
                ierrs(i))
            ! Interpolate intensitites
                ! print*, "ierr", ierr
                ! print*, "First and last intensity in"
                ! print*, intensities_in(cut_E_start_beams(i),i)
                ! print*, intensities_in(cut_E_start_beams(i)+ cut_n_E_beams(i), i)
                ! print*, "Any NAN?", ANY(ieee_is_nan(intensities_in))
                ! print*, "First and last knot"
                ! print*, knots(1,i)
                ! print*, knots(n_knots_beams(i), i)
                ! print*, "n_knots, nt, deg, LHS_rows", n_knots_beams(i), nt_beams(i), deg, LHS_rows_beams(i)
            call single_interpolate_coeffs_to_grid( &
                n_knots_beams(i), knots(1:n_knots_beams(i), i), nt_beams(i), coeffs(1:nt_beams(i), i), &
                deg, &
                n_E_beams_out(i), &
                E_grid_out(E_start_beams_out(i) : E_start_beams_out(i) + n_E_beams_out(i)), &
                0, &
                intpol_intensity(E_start_beams_out(i) : E_start_beams_out(i) + n_E_beams_out(i), i))
            !print*, "Intpol"
            !print*, intpol_intensity(E_start_beams_out(i) : E_start_beams_out(i) + n_E_beams_out(i),i)
                ! Interpolate derivatives
            call single_interpolate_coeffs_to_grid( &
                n_knots_beams(i), knots(1:n_knots_beams(i), i), nt_beams(i), coeffs(1:nt_beams(i), i), &
                deg, &
                n_E_beams_out(i), &
                E_grid_out(E_start_beams_out(i) : E_start_beams_out(i) + n_E_beams_out(i)), &
                1, &
                intpol_derivative(E_start_beams_out(i) : E_start_beams_out(i) + n_E_beams_out(i), i))

            ! Calculate pendry Y function
            call pendry_y( &
                n_E_beams_out(i), &
                intpol_intensity(E_start_beams_out(i):E_start_beams_out(i)+n_E_beams_out(i), i), &
                intpol_derivative(E_start_beams_out(i):E_start_beams_out(i)+n_E_beams_out(i), i), &
                v0i, &
                y_func(E_start_beams_out(i):E_start_beams_out(i)+n_E_beams_out(i), i) &
                )
        end do
        ierr = MAXVAL(ierrs)
    end if

    !###############
    ! Y function calculation on new grid
    !###############

    ! TODO: clean up implementation

    ! output is the y functions


    print*, "New verison successfull"

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


! subroutine limit_range_index(in_beams, n_beams, n_energies, E_min_current, E_step, beam_starts, NE_beams, &
!                        NE_out_start, NE_out_stop, out_beams, new_NE, new_NE_beams, new_start_id)
!     !Limit_range:
!     !INPUT: E_min, E_max, array[I, E_min, E_step, NE]
!     !RETURNS: Beam cut to only within [E_min, E_max]

!     ! Takes beams recorded on the same energy grid and cut off any parts below or above a certain threshold Energy.
!     ! Returns beams on uniform size array, cut to right size. Also returns E_min of each beam + new min_index and cut NE
!     ! Works for one or for multiple beams (though the latter may defeat the purpose)

!     integer, intent(in) :: n_beams, n_energies ! number of beams, number of energies in in_beams
!     real(8), intent(in) :: in_beams(n_energies, n_beams)
!     real(8), intent(out) :: out_beams(NE_out_stop-NE_out_start+1, n_beams)

!     real(8), intent(in) :: E_min_current, E_step
!     integer, intent(in) :: NE_out_start, NE_out_stop
!     integer, intent(in)  :: beam_starts(n_beams), NE_beams(n_beams)
!     integer, intent(out) :: new_NE, new_NE_beams(n_beams)


!     ! Internal
!     integer steps, to_new_min, i, new_end_id(n_beams)
!     real new_max, new_min, new_start_id(n_beams)

!     new_NE = NE_out_stop-NE_out_start+1
    

!     ! allocate out_beams to right size
!     !allocate(out_beams(new_NE, n_beams))
!     out_beams = in_beams(NE_out_start : NE_out_stop, :)

!     ! Determine new boundary indices E_min& NE:
!     do i=0, n_beams
!         new_start_id(i) = max(beam_starts(i) - NE_out_start, 1)
!         new_end_id(i) = min(beam_starts(i) + NE_beams(i)-NE_out_start, NE_out_stop)
!         new_NE_beams(i) = new_end_id(i) - new_start_id(i) +1
!     end do

! end subroutine limit_range_index


! subroutine pre_evaluate_grid_beam(n, x, deg, grid_origin, grid_target, grid_pre_eval, ierr)
!     integer, INTENT(IN) :: deg
!     type(grid), INTENT(IN) :: grid_origin
!     type(grid), INTENT(IN) :: grid_target
    

!     type(grid_pre_evaluation), INTENT(OUT) :: grid_pre_eval
!     integer, INTENT(OUT) :: ierr

!     ! Internal
!     integer :: max_n_knots, max_nt, max_LHS_rows
!     integer :: n_knots_beams(n_beams), nt_beams(n_beams), LHS_rows_beams(n_beams)
!     real(8) :: knots(max_n_knots, n_beams), LHS(max_LHS_rows, max_nt, n_beams), coeffs(max_nt, n_beams)

!     call get_array_sizes(n_E, deg, max_n_knots, max_nt, max_LHS_rows)

!     call perform_checks(n_E, E_grid, n_E_out, E_grid_out, ierr)



!     do i=1,n_beams
!         call get_array_sizes(n_E_beams(i), deg, n_knots_beams(i), nt_beams(i), LHS_rows_beams(i))
!     end do

!     call single_calc_spline(n_E_beams(i), E_beams(:), ...)

!     ALLOCATE(knots(n_knots, n_beams), LHS(LHS_rows, nt, n_beams), coeffs(nt))

!     call pre_evaluate_grid(grid_origin%x, grid_origin%n, grid_target%x, grid_target%n, &
!                             deg, 1, grid_pre_eval, ierr) ! the 1 means do_checks ON
!     RETURN
! end subroutine pre_evaluate_grid_beam


! subroutine pre_evaluate_beam_deriv(grid_pre_eval, nu, grid_deriv_pre_eval, ierr)
!     integer, INTENT(IN) :: nu
!     type(grid_pre_evaluation), INTENT(IN) :: grid_pre_eval

!     type(deriv_pre_evaluation), INTENT(OUT) :: grid_deriv_pre_eval
!     integer, INTENT(OUT) :: ierr

!     call pre_evaluate_deriv(grid_pre_eval, nu, grid_deriv_pre_eval, ierr)
! end subroutine pre_evaluate_beam_deriv


! subroutine interpolate_beam(intensity, grid_pre_eval, interpolated_intensity, coeffs, ierr)
!     type(grid_pre_evaluation), INTENT(IN) :: grid_pre_eval
!     real(8), INTENT(IN) :: intensity(grid_pre_eval%grid_origin%n)

!     real(8), INTENT(OUT) :: interpolated_intensity(grid_pre_eval%grid_target%n)
!     integer, INTENT(OUT) :: ierr
!     real(8), INTENT(INOUT) :: coeffs(grid_pre_eval%infos%nt)

!     real(8), ALLOCATABLE :: coeffs_tmp(:)

!     call interpolate_fast(grid_pre_eval%grid_origin%n, intensity, grid_pre_eval, interpolated_intensity, coeffs_tmp, ierr)
!     coeffs = coeffs_tmp
! end subroutine interpolate_beam

! subroutine interpolate_deriv_beam(grid_pre_eval, deriv_pre_eval, coeffs, intensity_deriv, ierr)

!     type(grid_pre_evaluation), INTENT(IN) :: grid_pre_eval
!     TYPE(deriv_pre_evaluation), INTENT(IN):: deriv_pre_eval
!     real(8), INTENT(IN) :: coeffs(grid_pre_eval%infos%nt)

!     real(8), INTENT(OUT) :: intensity_deriv(grid_pre_eval%grid_target%n)
!     integer, INTENT(OUT) :: ierr

!     call interpolate_deriv_fast(grid_pre_eval%infos, grid_pre_eval%grid_target%n, grid_pre_eval%intervals, coeffs, &
!                                 deriv_pre_eval%deriv_deBoor_matrix, intensity_deriv, ierr)
!     RETURN
! end subroutine interpolate_deriv_beam



! Library functions    
! get rid of below - not required due to interpolation!


subroutine pendry_y(n_data, intensity, derivative, v0i, y_func)
    integer, INTENT(IN) :: n_data
    real(8), INTENT(IN) :: intensity(n_data)
    real(8), INTENT(IN) :: derivative(n_data)
    real(8), INTENT(IN) :: v0i

    real(8), INTENT(OUT) :: y_func(n_data)

    y_func = intensity*derivative /(intensity*intensity + v0i*v0i*derivative*derivative)

end subroutine pendry_y


pure function parabolic_optimize(values) result(minimum_pair)
    implicit none

    real(8), intent(in) :: values(:)
    real(8) :: minimum_pair(2) ! Don't declare as intent(out)
    integer known_values

    real(8) minium_pair(2)

    known_values = size(values)

    if (known_values < 3) then
        ! Not allowed; Problem since Parabola fit is impossible - return something?
    else
        ! Should be the case!
    end if

    minimum_pair = (2d0, 2d0)
    return
end function parabolic_optimize

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


end module r_factor_new