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

module interpolation
    use, intrinsic :: iso_fortran_env
    implicit none
    
    integer, parameter :: sp = REAL32
    integer, parameter :: dp = REAL64
    integer, parameter :: qp = REAL128

    integer, PARAMETER :: info_size=10
    

    type :: info_values
        integer :: deg
        integer :: n_knots
        integer :: nt
        integer :: kl, ku
        integer :: nleft, nright
        integer :: LHS_rows, LHS_cols
        integer :: RHS_cols
    end type info_values
    contains



    subroutine pack_values(deg, n_knots, nt, kl, ku, nleft, nright, LHS_rows, LHS_cols, RHS_cols, info)
        integer, INTENT(IN)                :: deg                  !! degree of spline
        integer, INTENT(IN)                :: n_knots              !! # knots
        integer, INTENT(IN)                :: nt, kl, ku           !! matrix shape specifiers
        integer, INTENT(IN)                :: LHS_rows, LHS_cols   !! shape LHS
        integer, INTENT(IN)                :: RHS_cols             !! shape RHS
        integer, INTENT(IN) :: nleft, nright
        
        type(info_values), INTENT(OUT)               :: info(info_size)
        info%deg       = deg
        info%n_knots   = n_knots
        info%nt        = nt
        info%kl        = kl
        info%ku        = ku
    info%nleft     = nleft 
    info%nright    = nright
    info%LHS_rows  = LHS_rows
    info%LHS_cols  = LHS_cols
    info%RHS_cols  = RHS_cols
    RETURN
end subroutine pack_values


subroutine unpack_values(info, deg, n_knots, nt, kl, ku, nleft, nright, LHS_rows, LHS_cols, RHS_cols)
    type(info_values), INTENT(IN)       :: info
    integer, INTENT(OUT)                :: deg                  !! degree of spline
    integer, INTENT(OUT)                :: n_knots              !! # knots
    integer, INTENT(OUT)                :: LHS_rows, LHS_cols   !! shape LHS
    integer, INTENT(OUT)                :: RHS_cols             !! shape RHS
    integer, INTENT(out) :: nleft, nright, nt, kl, ku
    deg         = info%deg
    n_knots     = info%n_knots
    nt          = info%nt
    kl          = info%kl
    ku          = info%ku
    nleft       = info%nleft
    nright      = info%nright
    LHS_rows    = info%LHS_rows
    LHS_cols    = info%LHS_cols
    RHS_cols    = info%RHS_cols
    RETURN
end subroutine unpack_values

end module interpolation
