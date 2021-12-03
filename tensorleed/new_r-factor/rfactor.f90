! Subroutines and functions for R factor calculation in ViPErLEED
!
! v0.1.0
!
! Author: Alexander M. Imre, 2021
! for license info see ViPErLEED Package

! Note: file extentsion .f95 required for f2py compatibility, but newer standard can be used.

module r_factor_new
    implicit none
    contains


! f2py -m rfactor rfactor.f90 -h rfactor.pyf --overwrite-signature --debug-capi
! f2py -c rfactor.pyf rfactor.f90

subroutine r_factor_beam(y1, size_y1, y2, size_y2, E_start1, E_start2, E_step, V0r_shift, R_pendry, &
                         numerator, denominator, N_overlapping_points)
    ! Innermost function Rfactor_beam:
    ! INPUT: Y1, Y2, Estart1, Estart2, V0rshift -> Y1 and Y2 are already prepared (i.e. same grid steps, I(V) was previously interpolated)
    ! DOES: figure out overlap region, calculate RPe
    ! RETURNS: RPe, Number of overlapping points/ !!! -> not just R pendry for each beam, but numerator & denominator! for further calculation!!
    !

    implicit none
    !###############
    !INPUTS
    !###############

    !f2py intent(in) size_y1, size_y2
    !f2py depend(y1) size_y1
    !f2py depend(y2) size_y2
    integer, intent (in)                    :: size_y1, size_y2
    !f2py intent(in) y1
    real, intent (in), dimension(size_y1)   :: y1
    !f2py intent(in) y2
    real, intent (in), dimension(size_y2)   :: y2
    real, intent (in)                       ::  E_start1, E_start2, E_step, V0r_shift
    ! outputs
    real, intent (out)                      ::  R_pendry, numerator, denominator
    integer,  intent (out)                  ::  N_overlapping_points

    real E_end1, E_end2, E_start, E_end, E_range
    real diff_start
    real, dimension (:), allocatable         :: y_diff, y_squared_sum
    integer y1_start_id, y2_start_id
    integer i


    ! end of energy ranges of Y functions
    E_end1 = E_start1+size_y1*E_step
    E_end2 = E_start2+size_y2*E_step

    E_start = max(E_start1, E_start2)
    E_end = min(E_end1, E_end2)


    ! Caclulate number of overlapping points - can be used as weight of R_factor outside
    E_range = E_end - E_start
    write(*,*) "E_range: ", E_range
    write(*,*) "E_start: ", E_start
    write(*,*) "E_step: ", E_step
    N_overlapping_points = int(E_range/E_step)
    write(*,*) "N_ovl: ", N_overlapping_points

    y1_start_id = int(E_start1-E_start)+1 !Fortan counting starts at 1
    y2_start_id = int(E_start2-E_start)+1

    ! Allocate array sizes; not technically required in Fortran 2003 standard
    allocate(y_diff(N_overlapping_points), y_squared_sum(N_overlapping_points))

    y_diff=y1(y1_start_id: y1_start_id + N_overlapping_points-1) - y2(y2_start_id: y2_start_id + N_overlapping_points-1) ! difference between Y functions
    y_squared_sum=y1(y1_start_id:y1_start_id+N_overlapping_points-1)**2+y2(y2_start_id:y2_start_id+N_overlapping_points-1)**2

    write(*,*) "y_1: ", y1
    write(*,*) "y_2: ", y2
    write(*,*) "y_diff: ", y_diff
    write(*,*) "y_squared_sum: ", y_squared_sum

    ! caclulate numerator = integral (Y1-Y2)**2 dE
    numerator = trapez_integration_const_dx(y_diff**2, E_step)
    ! calculate denominator = integral (Y1**2+Y2**2) dE
    denominator = trapez_integration_const_dx(y_squared_sum, E_step)

    R_pendry = numerator/denominator
    return

end subroutine r_factor_beam


!V0rshift not yet implemented
    ! beamtypes not yet implemented -> does it really make sense to do that here even?
    ! maybe extra subroutine actually
subroutine Rfactor_beamset(y1, sizes_y1, y2, sizes_y2, E_start1, E_start2, nr_beams, E_step, V0rshift, &
                           R_Pe_weighted, R_Pe_beams, N_overlapping_points)
    !Rfactor_beamset:
    !INPUT: set of arrays Y1, Y2, Estart1, Estart2, V0rshift, beamtypes -> calls the above in a loop over beams
    !DOES: call above, average based on overlapping points, also split into types ("integer/fractional")
    !RETURNS: array RPe [overall and for each beamtype]
    implicit none

    !###############
    !INPUTS
    !###############
!f2py integer, hidden, intent(in), depend(E_start1, E_start2), check(len(E_start1)==len(E_start2)) :: nr_beams=len(E_start1)
    integer nr_beams
!f2py real, intent(in) :: E_step
    real E_step
!f2py integer, intent(in) :: sizes_y1
    integer, intent(in) :: sizes_y1(nr_beams)
!f2py integer, intent(in) :: sizes_y2
    integer, intent(in) :: sizes_y2(nr_beams)
!f2py real intent(in) :: y1(:,:), y2(:,:)
    real, intent(in) :: y1 (:,:), y2 (:,:)
!f2py intent(in) E_start_1
!f2py intent(in) E_start_2
    real, intent(in)    :: E_start1(nr_beams), E_start2(nr_beams)
!f2py integer, intent(out) :: N_overlapping_points
    integer, intent(out) :: N_overlapping_points(nr_beams)
!f2py real, intent(out) :: R_Pe_beams(nr_beams), R_Pe_weighted
    real, intent(out) :: R_Pe_beams(nr_beams), R_Pe_weighted
!f2py real, intent(in)::  V0rshift
    real, intent(in) :: V0rshift

    real numerators(nr_beams), denominators(nr_beams)

    ! temporay variables used in subroutine
    integer i
    integer total_points
    real temp_num, temp_denom

    ! iterate over beams and call subroutine r_factor_beam
    do i=1, nr_beams
        call r_factor_beam(y1(i,:), sizes_y1(i), y2(i,:), sizes_y2(i), E_start1(i), E_start2(i), E_step, &
                           V0rshift, R_Pe_beams(i), numerators(i), denominators(i), N_overlapping_points(i))
    end do
    total_points = sum(N_overlapping_points)
    R_Pe_weighted = sum(numerators/denominators*N_overlapping_points)/total_points

    return
end subroutine Rfactor_beamset

subroutine Rfactor_beamtypes(y1, sizes_y1, y2, sizes_y2, E_start1, E_start2, nr_beams, E_step, V0rshift, &
                           beamtypes, nr_beamtypes, R_Pe_weighted, R_Pe_beams, N_overlapping_points)
    !###############
    !INPUTS
    !###############
!f2py integer, hidden, intent(in), depend(E_start1, E_start2), check(len(E_start1)==len(E_start2)) :: nr_beams=len(E_start1)
    integer nr_beams
!f2py real, intent(in) :: E_step
    real E_step
!f2py integer, intent(in) :: sizes_y1
    integer, intent(in) :: sizes_y1(nr_beams)
!f2py integer, intent(in) :: sizes_y2
    integer, intent(in) :: sizes_y2(nr_beams)
!f2py real intent(in) :: y1(:,:), y2(:,:)
    real, intent(in) :: y1 (:,:), y2 (:,:)
!f2py intent(in) E_start_1
!f2py intent(in) E_start_2
    real, intent(in)    :: E_start1(nr_beams), E_start2(nr_beams)
    integer, intent(in) :: beamtypes(nr_beams)
    integer, intent(in) :: nr_beamtypes
!f2py integer, intent(out) :: N_overlapping_points(:)
    integer, intent(out) :: N_overlapping_points(nr_beams)
!f2py real, intent(out) :: R_Pe_beams(nr_beams), R_Pe_weighted(nr_beamtypes)
    real, intent(out) :: R_Pe_beams(nr_beams), R_Pe_weighted(nr_beamtypes)
!f2py real, intent(in)::  V0rshift
    real, intent(in) :: V0rshift

    ! variables used internally
    integer beam, group, points_group, tmp_points(nr_beams)
    !
    real tmp_num(nr_beams), tmp_denom(nr_beams), tmp_pendry


    ! beamtypes is array that contains assignment o beamtype group for each beam
    ! contains integers from 1 to nr_beamtypes

    ! Loop over beamtype groups
    do group=1, nr_beamtypes
        tmp_points = 0
        tmp_num = 0d0
        tmp_denom = 0d0
        do beam = 1, nr_beams
            if (beamtypes(beam)==group) then
                call r_factor_beam(y1(beam,:), sizes_y1(beam), y2(beam,:), sizes_y2(beam), E_start1(beam), &
                        E_start2(beam), E_step, V0rshift, tmp_pendry, tmp_num(beam), tmp_denom(beam), tmp_points(beam))
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

    real V0, R_Pe_V0
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


end subroutine Rfactor_v0ropt

subroutine prepare_beams(intensity, n_beams, NE_in, E_min_global, E_min, E_step, NE_beams, &
                         skip_stages, averaging_scheme, n_averaging_types, &
                         NE_after, E_min_after, E_step_after)
    !Prepare_beams:
    !INPUT: array[I, E_min, E_step, NE], E_grid_step, averaging_scheme, smoothing?, E_min, E_max
    !DOES: (0) Limit_range, (1) Average/discard/reorder according to scheme; (2) smooth?; (3) interpolate on grid; (4) compute Y on new grid
    !Probably call several functions... see existing PREEXP, similar!
    !RETURNS: array of Y

    !###############
    !INPUTS
    !###############
    integer, intent(in) :: n_beams, NE_in ! number of beams, number of passed energies
    integer, intent(in) :: stages(5) ! which stages to execute
    real, intent(in) :: intensity(NE_in,n_beams), E_min_global ! beam intensities
    real, intent(in) :: E_min(n_beams), E_step(n_beams), NE_beams(n_beams) ! minimum energies, nr of steps
    integer, intent(in) :: averaging_scheme(n_beams)
    integer, intent(in) :: smoothing ! smooth or not?
    integer, intent(in) :: NE_after ! number of steps (size E_grid_step) after interpolation
    real, intent(in) :: E_min_after, E_step_after ! desired grid step after interpolation

    !###############
    ! Internal
    !###############
    real :: E_max_after
    integer :: new_start_stop_step(3)
    real, allocatable :: cut_intensity(:,:)

    !###############
    ! Limit range of beams
    !###############
    ! Everything assumes, experimental beams are recorded on the same energy grid already. This should always be the case.
    ! Thus only the information of E_min and NE are required for each beam
    if (skip_stages(1).EQ.0) then
    E_max_after = E_min_after + NE_after*E_step_after
        ! range_index_from_Energy will calculate size of intensity array after cutting to desired energy range
        call range_index_from_Energy(E_min_global, NE_in, E_step, E_min_afer, E_max_after, new_start_stop_step)

        allocate(cut_intensity(new_start_stop_step(3),n_beams))

        ! limit_range_index cuts range of intensity array
        call limit_range_index(intensity, n_beams, n_energies, E_min_current, E_step, beam_starts, NE_beams, &
                           new_start_stop_step(1), new_start_stop_step(2), cut_intensity, new_NE, new_NE_beams, new_start_id)
    else
        new_start_stop_step(3) = NE_in
    end if

    !###############
    ! Average/discard/reorder
    !###############

    if (skip_stages(2).EQ.0) then

        allocate(averaged_intensity(new_start_stop_step(3), n_averaging_types))
        call avg_scheme(cut_intensity, n_beams, NE, averaging_scheme, n_averaging_types, averaged_intensity)
        n_beams = n_averaging_types
    end if


    !###############
    ! Smoothing
    !###############

    if (skip_stages(3).EQ.0) then
        write(*,*) "Smoothing not yet implemented!" ! TODO: implement smoothing in prepare_beams
    end if

    !###############
    ! Interpolation on new grid
    !###############

    if (skip_stages(4).EQ.0) then
    ! TODO: implement
    end if

    !###############
    ! Y function calculation on new grid
    !###############


    ! TODO: implement


end subroutine prepare_beams


subroutine range_index_from_Energy(E_min_current, NE_in, E_step, E_min_cut, E_max_cut, new_start_stop_step)
    integer, intent(in) :: NE_in
    real, intent(in) :: E_min_current, E_min_cut, E_max_cut, E_step

    integer , intent(out) :: new_start_stop_step(3)

    real :: old_max, new_max, new_min

    old_max = E_min_current + E_step*NE_in
    new_max = min(E_max_cut, old_max)
    new_min = max(E_min_current, E_min_cut)

    new_start_stop_step(3) = floor((E_max_cut-E_min_cut)/E_step) ! typecast to integer; will cut off float...

    new_start_stop_step(1) = (E_min_cut - E_min_current)/E_step
    new_start_stop_step(2) = new_start_stop_step(1) + new_start_stop_step(3)
    return
end subroutine range_index_from_Energy


subroutine limit_range_index(in_beams, n_beams, n_energies, E_min_current, E_step, beam_starts, NE_beams, &
                       NE_out_start, NE_out_stop, out_beams, new_NE, new_NE_beams, new_start_id)
    !Limit_range:
    !INPUT: E_min, E_max, array[I, E_min, E_step, NE]
    !RETURNS: Beam cut to only within [E_min, E_max]

    ! Takes beams recorded on the same energy grid and cut off any parts below or above a certain threshold Energy.
    ! Returns beams on uniform size array, cut to right size. Also returns E_min of each beam + new min_index and cut NE
    ! Works for one or for multiple beams (though the latter may defeat the purpose)

    integer, intent(in) :: n_beams, n_energies ! number of beams, number of energies in in_beams
    real, intent(in) :: in_beams(n_energies, n_beams)
    real, intent(out) :: out_beams(NE_out_stop-NE_out_start+1, n_beams)

    real, intent(in) :: E_min_current, E_step
    integer, intent(in) :: NE_out_start, NE_out_stop
    integer, intent(in)  :: beam_starts(n_beams), NE_beams(n_beams)
    integer, intent(out) :: new_NE, new_NE_beams(n_beams)


    ! Internal
    integer steps, to_new_min, i, new_end_id(n_beams)
    real new_max, new_min, new_start_id(n_beams)

    new_NE = NE_out_stop-NE_out_start+1
    

    ! allocate out_beams to right size
    !allocate(out_beams(new_NE, n_beams))
    out_beams = in_beams(NE_out_start : NE_out_stop, :)

    ! Determine new boundary indices E_min& NE:
    do i=0, n_beams
        new_start_id(i) = max(beam_starts(i) - NE_out_start, 1)
        new_end_id(i) = min(beam_starts(i) + NE_beams(i)-NE_out_start, NE_out_stop)
        new_NE_beams(i) = new_end_id(i) - new_start_id(i) +1
    end do

end subroutine limit_range_index

subroutine avg_scheme(in_beams, n_beams, NE, averaging_scheme, avg_types averaged_beams)
    real, intent(in), dimension(NE, n_beams) :: in_beams
    integer, intent(in) :: n_beams, NE, avg_types
! avg_types = max(averging_scheme)
    real, intent(out), dimension(NE, avg_types) :: averaged_beams
    ! any beam with scheme 0 is discarded (intentionally)

    ! Internal
    integer, dimension(avg_types) :: avg_counter
    ! initialize
    avg_counter(:) = 0
    averaged_beams(:,:) = 0

    do concurrent (i = 1:avg_types)
        do j = 1,n_beams+1
            if (avg_scheme(j) == i) then
                averaged_beams(:,i) = averaged_beams(:,i) + in_beams(:,j)
                avg_counter(i) = avg_counter(i) +1
            end if
        end do
        if (avg_counter(i)>0) then
        averaged_beams(:,i) = averaged_beams(:,i)/avg_counter(i)
        end if
    end do

end subroutine avg_scheme


    ! Library functions

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
    !external :: derivative

    call derivative(intensity, de, i_prime, n_data)
    y_func = intensity * i_prime / (intensity**2 + v0i**2 * i_prime**2)  ! TODO: Maybe can be better optimized?

end subroutine pendry_y

pure function parabolic_optimize(values) result(minimum_pair)
    implicit none

    real, intent(in) :: values(:)
    real :: minimum_pair(2) ! Don't declare as intent(out)
    integer known_values

    real minium_pair(2)

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
    real, intent(in)    :: f(:)
    real, intent(in)    :: dx
    real integral
    integer n

    n = size(f)
    integral = 0
    integral = sum(f(1:n-1)+f(2:n))*dx/2
    return
end function trapez_integration_const_dx

pure function test_exec(input) result(output)
    implicit none
    real, intent(in) :: input
    real output

    output = input

    return
end function test_exec

end module r_factor_new