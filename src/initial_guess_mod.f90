!===============================================================================
! Module: initial_guess_mod
! Description: Analytical approximations (MSA for charged spheres,
!              Percus-Yevick for hard spheres) to provide initial solutions.
!===============================================================================
module initial_guess_mod
  use types_mod
  use params_mod
  use grid_mod
  use fft_mod
  implicit none
  public

contains

  !=============================================================================
  ! Subroutine: msa
  ! Computes the Mean Spherical Approximation (MSA) solution for equal sized
  ! charged hard spheres. Sout(i, j) contains the series function h - c.
  ! Programmed by E. Lomba, Trondheim 23.02.90.
  !=============================================================================
  subroutine msa(sout, nr_pts, rhoe, qq_ions, tr, sigma_val)
    integer(ip), intent(in) :: nr_pts
    real(dp), dimension(2), intent(in) :: rhoe, qq_ions
    real(dp), intent(in) :: tr, sigma_val
    real(dp), dimension(nr_pts, NIT), intent(out) :: sout

    integer(ip) :: i, i1, i2, mx, n_pts
    real(dp) :: b, ck, cq, cr, cst, ct, e, rho_eff, rk, rlk, sk_val, v, vk_val, x, xe, k_val
    real(dp), dimension(nr_pts) :: ts_arr, s_arr, ss_arr

    cr = CR1
    cq = CQ1
    rk = -1.0_dp / tr
    rho_eff = sum(rhoe(1:2) * qq_ions(1:2)**2)

    e = rho_eff * PI / 6.0_dp
    x = sqrt(24.0_dp * e / tr)
    if (x > 1.0_dp) then
       xe = 4.0_dp * x
    else
       xe = 2.0_dp
    end if
    b = (1.0_dp + x - sqrt(1.0_dp + 2.0_dp * x)) / x
    v = 2.0_dp * b / tr
    n_pts = nr_pts - 1

    do i = 1, nr_pts
       k_val  = Q(i)
       sk_val = sin(k_val)
       ck     = cos(k_val)
       ct = ( v * (-k_val * ck + sk_val) + v**2 * (-k_val**2 * ck + 2.0_dp * k_val * sk_val &
            + 2.0_dp * (ck - 1.0_dp)) / (4.0_dp * rk * k_val) ) / k_val**2 - rk * ck / k_val
       ct = -FOURPI * ct

       rlk = FOURPI * rk * (1.0_dp / k_val - k_val / (k_val**2 + xe**2))
       cst = ct - rlk

       vk_val    = k_val * rho_eff / (k_val - rho_eff * rlk)
       ts_arr(i) = k_val * cst * (vk_val / rho_eff)**2 / (k_val - cst * vk_val) - cst
    end do

    call fftsin(ts_arr, cq, s_arr)
    call pyhs(ss_arr, nr_pts, rho_eff, 1.0_dp)

    mx = 1
    do i1 = 1, 2
       do i2 = 1, i1
          do i = 1, n_pts
             sout(i, mx) = ss_arr(i) + qq_ions(i1) * qq_ions(i2) * s_arr(i)
          end do
          mx = mx + 1
       end do
    end do
  end subroutine msa

  !=============================================================================
  ! Subroutine: pyhs
  ! Computes the hard sphere series function r * s(r) in Percus-Yevick approx.
  ! Programmed by E. Lomba, Trondheim 20.09.90.
  !=============================================================================
  subroutine pyhs(s_out, nr_pts, rho_val, dhs)
    integer(ip), intent(in) :: nr_pts
    real(dp), intent(in) :: rho_val, dhs
    real(dp), dimension(nr_pts), intent(out) :: s_out

    integer(ip) :: k, n_pts
    real(dp) :: ck, cpi, cq, cr, et, f1, f2, f3
    real(dp) :: qi, qi2, qk, qr, qr2, qr3, qr4, rl1, rl2, rtc, sk_val, tdcf_val
    real(dp), dimension(nr_pts) :: ts_arr

    cr  = CR1
    cq  = CQ1
    et  = PI * rho_val * dhs**3 / 6.0_dp
    cpi = -FOURPI * dhs * dhs
    rl1 =  (1.0_dp + 2.0_dp * et)**2 / (1.0_dp - et)**4
    rl2 = -(1.0_dp + 0.5_dp * et)**2 / (1.0_dp - et)**4
    n_pts = nr_pts - 1

    do k = 1, n_pts
       qk    = Q(k)
       qi    = qk * dhs
       qi2   = qi * qi
       sk_val= sin(qi)
       ck    = cos(qi)
       qr    = 1.0_dp / qi
       qr2   = qr * qr
       qr3   = qr2 * qr
       qr4   = qr2 * qr2
       f1    = (sk_val - qi * ck) * qr2
       f2    = (-2.0_dp + 2.0_dp * qi * sk_val - qi2 * ck + 2.0_dp * ck) * qr3
       f3    = ( 4.0_dp * (6.0_dp + qi * (qi2 - 6.0_dp) * sk_val) &
               + (12.0_dp * qi2 - qi2 * qi2 - 24.0_dp) * ck ) * qr4 * qr
       tdcf_val = cpi * (rl1 * f1 + 6.0_dp * et * rl2 * f2 + 0.5_dp * et * rl1 * f3)
       rtc   = rho_val * tdcf_val
       ts_arr(k) = tdcf_val * rtc / (qk - rtc)
    end do

    call fftsin(ts_arr, cq, s_out)
  end subroutine pyhs

end module initial_guess_mod
