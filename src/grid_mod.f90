!===============================================================================
! Module: grid_mod
! Description: Radial and reciprocal grid vectors, spacings, and Fourier factors.
!===============================================================================
module grid_mod
  use types_mod
  use params_mod
  implicit none
  public

  ! Real-space and reciprocal-space grid arrays
  real(dp), dimension(:), allocatable :: R  ! R(i) = i * DR
  real(dp), dimension(:), allocatable :: Q  ! Q(i) = i * DQ

  ! Grid parameters
  real(dp) :: DR   ! Real-space step size
  real(dp) :: DQ   ! Reciprocal-space step size = pi / (NR * DR)
  real(dp) :: CR1  ! Forward transform normalization = pi * DR
  real(dp) :: CQ1  ! Backward transform normalization = DQ / (8 * pi^2)

contains

  subroutine init_grid(nr_in, dr_in)
    integer(ip), intent(in) :: nr_in
    real(dp), intent(in) :: dr_in
    integer(ip) :: i

    DR  = dr_in
    DQ  = PI / (real(nr_in, dp) * DR)
    CR1 = PI * DR
    CQ1 = DQ / (8.0_dp * PI**2)

    if (allocated(R)) deallocate(R)
    if (allocated(Q)) deallocate(Q)
    allocate(R(nr_in), Q(nr_in))

    do i = 1, nr_in
       R(i) = real(i, dp) * DR
       Q(i) = real(i, dp) * DQ
    end do
  end subroutine init_grid

  subroutine cleanup_grid()
    if (allocated(R)) deallocate(R)
    if (allocated(Q)) deallocate(Q)
  end subroutine cleanup_grid

end module grid_mod
