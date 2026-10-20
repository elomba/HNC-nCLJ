!===============================================================================
! Module: fft_mod
! Description: Fast Sine and Cosine Fourier Transforms with precomputed
!              trigonometric tables and bit-reversal lookup arrays.
! Reference: E.O. Brigham, "The Fast Fourier Transform", Chap. 10; F. Lado.
!===============================================================================
module fft_mod
  use types_mod
  use params_mod
  implicit none
  public

  ! FFT dimensions
  integer(ip) :: N, NR, NU

  ! Precomputed tables
  real(dp), dimension(:), allocatable :: SINe, COSe
  real(dp), dimension(:), allocatable :: FR, FI
  integer(ip), dimension(:), allocatable :: ibit_rev

contains

  subroutine init_fft(nr_in)
    integer(ip), intent(in) :: nr_in
    integer(ip) :: i
    real(dp) :: a1

    NU = int(log(real(nr_in, dp)) / log(2.0_dp) + 0.540_dp, ip)
    NR = 2**NU
    N  = NR - 1

    if (allocated(SINe))     deallocate(SINe)
    if (allocated(COSe))     deallocate(COSe)
    if (allocated(FR))       deallocate(FR)
    if (allocated(FI))       deallocate(FI)
    if (allocated(ibit_rev)) deallocate(ibit_rev)

    allocate(SINe(NR), COSe(NR), FR(NR), FI(NR), ibit_rev(0:NR-1))

    a1 = PI / real(NR, dp)
    do i = 1, NR
       SINe(i) = sin(real(i, dp) * a1)
       COSe(i) = cos(real(i, dp) * a1)
    end do

    do i = 0, NR - 1
       ibit_rev(i) = ibitr_scalar(i, NU)
    end do
  end subroutine init_fft

  pure function ibitr_scalar(j_val, nu_val) result(res)
    integer(ip), intent(in) :: j_val, nu_val
    integer(ip) :: res
    integer(ip) :: i, j1, j2

    j1 = j_val
    res = 0
    do i = 1, nu_val
       j2 = j1 / 2
       res = res + res + (j1 - j2 - j2)
       j1 = j2
    end do
  end function ibitr_scalar

  subroutine fftsin(F, C, Tf)
    real(dp), dimension(NR), intent(inout) :: F
    real(dp), intent(in) :: C
    real(dp), dimension(NR), intent(out)   :: Tf

    integer(ip) :: i, ii, k, k0, k1, k2, l, nr1
    real(dp) :: ti, tr

    k0 = NR / 2
    nr1 = NR + 1
    FR(1) = 0.0_dp
    F(NR) = 0.0_dp

    do i = 1, k0
       FR(i + 1)   =  F(2 * i)
       FI(i)       =  F(2 * i - 1)
       FR(nr1 - i) = -FR(i + 1)
       FI(nr1 - i) = -FI(i)
    end do

    k = 0
    do l = 1, NU
       do i = 1, k0
          k1 = k + 1
          k2 = k1 + k0
          tr = FR(k2)
          ti = FI(k2)
          FR(k2) = FR(k1) - tr
          FI(k2) = FI(k1) - ti
          FR(k1) = FR(k1) + tr
          FI(k1) = FI(k1) + ti
          k = k + 1
       end do

       do while (.true.)
          k = k + k0
          if (k < NR) then
             ii = 2 * ibit_rev(k / k0)
             do i = 1, k0
                k1 = k + 1
                k2 = k1 + k0
                tr = FR(k2) * COSe(ii) + FI(k2) * SINe(ii)
                ti = FI(k2) * COSe(ii) - FR(k2) * SINe(ii)
                FR(k2) = FR(k1) - tr
                FI(k2) = FI(k1) - ti
                FR(k1) = FR(k1) + tr
                FI(k1) = FI(k1) + ti
                k = k + 1
             end do
          else
             k = 0
             k0 = k0 / 2
             exit
          end if
       end do
    end do

    do k = 1, NR
       i = ibit_rev(k - 1) + 1
       if (i > k) then
          tr = FR(k)
          ti = FI(k)
          FR(k) = FR(i)
          FI(k) = FI(i)
          FR(i) = tr
          FI(i) = ti
       end if
    end do

    do i = 1, N
       Tf(i) = C * ( -(FI(i + 1) - FI(nr1 - i)) &
                   +  (FR(i + 1) - FR(nr1 - i)) * COSe(i) &
                   +  (FI(i + 1) + FI(nr1 - i)) * SINe(i) )
    end do
  end subroutine fftsin

  subroutine fftcos(F0, F, Tf)
    real(dp), intent(in) :: F0
    real(dp), dimension(NR), intent(inout) :: F
    real(dp), dimension(NR), intent(out)   :: Tf

    integer(ip) :: i, ii, k, k0, k1, k2, l, nr1
    real(dp) :: ti, tr

    k0 = NR / 2
    nr1 = NR + 1
    FR(1) = F0
    F(NR) = 0.0_dp

    do i = 1, k0
       FR(i + 1)   = F(2 * i)
       FI(i)       = F(2 * i - 1)
       FR(nr1 - i) = FR(i + 1)
       FI(nr1 - i) = FI(i)
    end do

    k = 0
    do l = 1, NU
       do i = 1, k0
          k1 = k + 1
          k2 = k1 + k0
          tr = FR(k2)
          ti = FI(k2)
          FR(k2) = FR(k1) - tr
          FI(k2) = FI(k1) - ti
          FR(k1) = FR(k1) + tr
          FI(k1) = FI(k1) + ti
          k = k + 1
       end do

       do while (.true.)
          k = k + k0
          if (k < NR) then
             ii = 2 * ibit_rev(k / k0)
             do i = 1, k0
                k1 = k + 1
                k2 = k1 + k0
                tr = FR(k2) * COSe(ii) + FI(k2) * SINe(ii)
                ti = FI(k2) * COSe(ii) - FR(k2) * SINe(ii)
                FR(k2) = FR(k1) - tr
                FI(k2) = FI(k1) - ti
                FR(k1) = FR(k1) + tr
                FI(k1) = FI(k1) + ti
                k = k + 1
             end do
          else
             k = 0
             k0 = k0 / 2
             exit
          end if
       end do
    end do

    do k = 1, NR
       i = ibit_rev(k - 1) + 1
       if (i > k) then
          tr = FR(k)
          ti = FI(k)
          FR(k) = FR(i)
          FI(k) = FI(i)
          FR(i) = tr
          FI(i) = ti
       end if
    end do

    do i = 1, N
       Tf(i) = FR(i + 1) + FR(nr1 - i) &
             + (FI(i + 1) + FI(nr1 - i)) * COSe(i) &
             - (FR(i + 1) - FR(nr1 - i)) * SINe(i)
    end do
  end subroutine fftcos

  subroutine cleanup_fft()
    if (allocated(SINe))     deallocate(SINe)
    if (allocated(COSe))     deallocate(COSe)
    if (allocated(FR))       deallocate(FR)
    if (allocated(FI))       deallocate(FI)
    if (allocated(ibit_rev)) deallocate(ibit_rev)
  end subroutine cleanup_fft

end module fft_mod
