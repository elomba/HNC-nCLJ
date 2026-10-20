!===============================================================================
! Module: potentials_mod
! Description: Computes truncated and shifted Lennard-Jones interaction potentials,
!              virial functions, and long-range Allnatt Coulomb screening components.
!===============================================================================
module potentials_mod
  use types_mod
  use params_mod
  use grid_mod
  use fft_mod
  use arrays_mod
  implicit none
  public

contains

  !=============================================================================
  ! Subroutine: compute_lj_potential
  ! Evaluates cut-and-shifted Lennard-Jones potential u_LJ(r) and virial w_LJ(r).
  !=============================================================================
  subroutine compute_lj_potential(tr)
    real(dp), intent(in) :: tr
    integer(ip) :: i, j
    real(dp) :: ep4, sig, r6, r12, r_cut_j

    ulj(:, :) = 0.0_dp
    wlj(:, :) = 0.0_dp

    do j = 1, NIT
       ep4 = eps(j) * 4.0_dp / tr
       sig = sigma(j)
       r_cut_j = rcut(j)

       ucut(j) = ep4 * (sig / r_cut_j)**6 * ((sig / r_cut_j)**6 - 1.0_dp)

       do i = 1, NR
          if (R(i) <= r_cut_j) then
             r6  = (sig / R(i))**6
             r12 = r6 * r6
             ulj(i, j) = ep4 * (r12 - r6) - ucut(j)
             wlj(i, j) = -ep4 * (12.0_dp * r12 - 6.0_dp * r6)
          end if
       end do
    end do
  end subroutine compute_lj_potential

  !=============================================================================
  ! Subroutine: compute_screening_and_lr
  ! Computes Debye inverse screening length, singularity suppression factors,
  ! long-range direct correlation tails clr(r), sclr(r), and chain sum v(k).
  !=============================================================================
  subroutine compute_screening_and_lr(tr, rkx, gama, rkd, xp, xm)
    real(dp), intent(in)  :: tr, rkx, gama
    real(dp), intent(out) :: rkd, xp, xm

    integer(ip) :: id1, id2, mx, i
    real(dp) :: xe, xp2, xm2, dnt, cv, csv, rkdx2, qi, ri

    ! Debye inverse screening length
    rkd = sqrt(FOURPI * rkx / tr)

    ! Singularity suppression parameter
    if (rkd > 1.0_dp) then
       xe = 4.0_dp * rkd
    else
       xe = 2.0_dp
    end if

    xp2 = (xe**2 - sqrt(xe**4 - 4.0_dp * (xe * rkd)**2)) / 2.0_dp
    xm2 = (xe**2 + sqrt(xe**4 - 4.0_dp * (xe * rkd)**2)) / 2.0_dp
    xp  = sqrt(xp2)
    xm  = sqrt(xm2)

    mx = 1
    do id1 = 1, NSP
       do id2 = 1, id1
          cr(mx)   = CR1 * sqrt(rho(id1) * rho(id2))
          cq(mx)   = CQ1 / sqrt(rho(id1) * rho(id2))
          dnt      = cr(mx) / CR1
          gamm(mx) = -qq(id1) * qq(id2) * gama / tr
          cv       = FOURPI * dnt * gamm(mx) * xe**2
          csv      = gamm(mx) * xe**2 / (xm2 - xp2)
          rkdx2    = (rkd * xe)**2

          do i = 1, NR
             qi = Q(i)
             ri = R(i)

             ! Long-range tails of c(r) and s(r)
             clr(i, mx)  = gamm(mx) * (1.0_dp - exp(-xe * ri))
             sclr(i, mx) = csv * (exp(-xp * ri) - exp(-xm * ri))

             ! Hypervertex chain-sum matrix v(k) and its inverse
             vk(i, mx)   =  cv / (qi**4 + (qi * xe)**2 + rkdx2)
             vki(i, mx)  = -cv / (qi**4 + (qi * xe)**2)
             if (ipg(mx)) then
                vk(i, mx)  = 1.0_dp + vk(i, mx)
                vki(i, mx) = 1.0_dp + vki(i, mx)
             end if
          end do
          mx = mx + 1
       end do
    end do
  end subroutine compute_screening_and_lr

end module potentials_mod
