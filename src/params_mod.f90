!===============================================================================
! Module: params_mod
! Description: Global system parameters, dimensions, and component counts.
!===============================================================================
module params_mod
  use types_mod
  implicit none
  public

  ! Dimension parameters set at runtime from input file
  integer(ip) :: NSP       ! Number of species (supported: 2 to 5+)
  integer(ip) :: NUM       ! Number of grid points for Newton-Raphson scheme
  integer(ip) :: NIT       ! Number of distinct interaction pairs = (NSP*NSP + NSP)/2
  integer(ip) :: ND        ! Dimension of Newton-Raphson Jacobian system = NIT * NUM
  integer(ip), parameter :: NTEMP = 2000 ! Maximum number of temperatures/densities

end module params_mod
