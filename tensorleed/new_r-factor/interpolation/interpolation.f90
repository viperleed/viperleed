! Created by Alexander M. Imre on 04.12.21.

module interpolation
    implicit none
    contains

    subroutine equidist_quint_spline_coeffs(in_data, n_in, coeffs)
        integer, (intent in) :: n_in
        real(dp), (intent in) :: in_data(n_in)

        real(dp), (intent out) :: coeffs(5*n_in)




    end subroutine equidist_spline_intpol_coeffs


end module interpolation