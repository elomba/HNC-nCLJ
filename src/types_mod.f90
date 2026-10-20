!===============================================================================
! Module: types_mod
! Description: Standard precision definitions and physical/mathematical constants.
!===============================================================================
module types_mod
  use iso_fortran_env, only: dp => real64, ip => int32
  implicit none
  public

  ! Mathematical constants in double precision
  real(dp), parameter :: PI     = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: TWOPI  = 2.0_dp * PI
  real(dp), parameter :: FOURPI = 4.0_dp * PI

end module types_mod
