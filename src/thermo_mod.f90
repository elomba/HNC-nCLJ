!===============================================================================
! Module: thermo_mod
! Description: Evaluation of thermodynamic properties: internal energy,
!              Helmholtz free energy, virial and free-energy pressure (compressibility
!              factors), excess chemical potentials, and Bhatia-Thornton structure
!              factors (S_NN, S_CC, S_NC).
!===============================================================================
module thermo_mod
  use types_mod
  use params_mod
  use grid_mod
  use fft_mod
  use arrays_mod
  use linalg_mod
  implicit none
  public

contains

  !=============================================================================
  ! Subroutine: skcc
  ! Evaluates Bhatia-Thornton concentration-concentration structure factor S_cc(q).
  !=============================================================================
  subroutine skcc(sk_in, rho_in, rhot_val, itcc, nsp_val, scc_out, nit_val)
    integer(ip), intent(in) :: itcc, nsp_val, nit_val
    real(dp), dimension(nit_val), intent(in) :: sk_in
    real(dp), dimension(nsp_val), intent(in) :: rho_in
    real(dp), intent(in)  :: rhot_val
    real(dp), intent(out) :: scc_out

    integer(ip) :: i, k, l, idx
    real(dp), dimension(nsp_val) :: x_frac
    real(dp) :: s11, s12, s22, x1, x2

    do i = 1, nsp_val
       x_frac(i) = rho_in(i) / rhot_val
    end do

    if (nsp_val == 2 .or. itcc == 2) then
       x1  = x_frac(1)
       x2  = x_frac(2)
       s11 = sk_in(1)
       s22 = sk_in(3)
       s12 = sk_in(2)
    else
       x1  = x_frac(1)
       x2  = sum(x_frac(2:nsp_val))
       s11 = sk_in(1)
       s12 = 0.0_dp
       s22 = 0.0_dp
       do k = 2, nsp_val
          idx = (k * k - k) / 2 + 1
          s12 = s12 + sk_in(idx)
          do l = 2, nsp_val
             if (k >= l) then
                idx = (k * k - k) / 2 + l
             else
                idx = (l * l - l) / 2 + k
             end if
             s22 = s22 + sk_in(idx)
          end do
       end do
    end if
    scc_out = x2**2 * s11 + x1**2 * s22 - 2.0_dp * x1 * x2 * s12
  end subroutine skcc

  !=============================================================================
  ! Subroutine: compute_thermodynamics
  ! Calculates complete thermodynamics after convergence of the HNC equations.
  !=============================================================================
  subroutine compute_thermodynamics(tr, rhot, rkd, pres, sek, sekc, ankt, sex, stot, &
                                    za, xc, lamb1, lamb2, vrr, vcc)
    real(dp), intent(in)  :: tr, rhot, rkd
    real(dp), intent(out) :: pres, sek, sekc, ankt, sex, stot, za, xc
    real(dp), intent(out) :: lamb1, lamb2, vrr, vcc

    integer(ip) :: i, j, k, l, ncj, ic, ics, is, js, nps
    real(dp) :: dnt2, thlr, sncj, rms1, rms2, esum, hij2, a1sum, rmusum, alrsum
    real(dp) :: ri, hij, sij, psum, eljsum, a2kt, qi, trachk, fklr, fklrii, thij, dethk
    real(dp) :: seklj, akcor, ipk, fcnr, xd, xion, cionion, Mrr0, Mcc0, Mrc0
    real(dp) :: sm1, sm11, sm3, chss, chsum, sigx, ssm1, eik, alr, asr, elr, xcsum
    real(dp) :: tds12_eff, sid

    xcsum = 0.0_dp

    do j = 1, NIT
       ncj  = ncore(j)
       dnt2 = (cr(j) / CR1)**2
       ic   = 0
       if (ipg(j)) ic = 1

       ! Restore direct correlation function (subtracting Coulomb tail)
       do i = 1, N
          thlr     = (vk(i, j) - real(ic, dp)) * Q(i)
          th(i, j)  = ts(i, j) + tdcf(i, j) + thlr
          dcf(i, j) = dcf(i, j) + clr(i, j) - gamm(j)
       end do

       ! Contributions to energy, a/nkt, and chemical potentials
       sncj = s(ncj, j) + sclr(ncj, j) - clr(ncj, j) + gamm(j)
       rms1 = -R(ncj) - sncj

       ! Backward Simpson integration for r <= sigma_ij
       do i = ncj - 1, 1, -1
          ics  = 2 * (mod(ncj - i, 2) + 1)
          rms1 = rms1 + ics * dcf(i, j) * R(i)
       end do

       esum   = -R(ncj) * gamm(j) * h(ncj, j)
       hij2   = (h(ncj, j) * R(ncj))**2
       a1sum  = -(h(ncj, j) + 1.0_dp) * dcf(ncj, j) * R(ncj)
       rmusum = R(ncj) * h(ncj, j) * sncj - 2.0_dp * R(ncj) * dcf(ncj, j)
       rms2   = dcf(ncj, j) * R(ncj)
       alrsum = hij2

       ! Forward Simpson integration for r > sigma_ij
       do i = ncj + 1, N
          ri     = R(i)
          hij    = h(i, j) * ri
          sij    = hij - dcf(i, j)
          ics    = 2 * (mod(i - ncj, 2) + 1)
          esum   = esum   - ics * gamm(j) * hij
          hij2   = hij * hij
          a1sum  = a1sum  - ics * (hij + ri) * dcf(i, j)
          alrsum = alrsum + ics * hij2
          rmusum = rmusum + ics * (hij * sij - 2.0_dp * ri * dcf(i, j))
          rms2   = rms2   + ics * dcf(i, j) * ri
       end do

       ! Virial pressure and energy sums for Lennard-Jones
       psum   = 0.0_dp
       eljsum = 0.0_dp
       do i = ncore(j), N - 2, 2
          l = i + 1
          k = i + 2
          psum = R(i)**2 * wlj(i, j) * (h(i, j) + 1.0_dp) &
               + 4.0_dp * R(l)**2 * wlj(l, j) * (h(l, j) + 1.0_dp) &
               + R(k)**2 * wlj(k, j) * (h(k, j) + 1.0_dp) + psum

          eljsum = R(i)**2 * ulj(i, j) * (h(i, j) + 1.0_dp) &
                 + 4.0_dp * R(l)**2 * ulj(l, j) * (h(l, j) + 1.0_dp) &
                 + R(k)**2 * ulj(k, j) * (h(k, j) + 1.0_dp) + eljsum
       end do

       pr(j)    = -2.0_dp * PI * psum * DR / 9.0_dp
       ekt(j)   =  2.0_dp * PI * esum * DR / 3.0_dp + PI * gamm(j) * sigma(j)**2
       ektlj(j) =  2.0_dp * PI * eljsum * DR / 3.0_dp
       a1kt(j)  =  2.0_dp * PI * a1sum * DR / 3.0_dp
       chm(j)   =  2.0_dp * PI * (rmusum - rms1) * DR / 3.0_dp
       xcsum    = (2.0_dp - real(ic, dp)) * (rms1 + rms2) * dnt2 + xcsum
       aklr(j)  =  2.0_dp * PI * alrsum * DR / 3.0_dp
    end do

    ! Integral in k-space for Helmholtz free energy a2kt
    a2kt = 0.0_dp
    do i = 1, N
       qi  = Q(i)
       ics = 2 * (mod(i, 2) + 1)
       trachk = 0.0_dp
       fklr   = 0.0_dp
       fklrii = 0.0_dp
       do j = 1, NIT
          thij = th(i, j) / qi
          ic = 0
          if (ipg(j)) then
             ic = 1
             trachk = trachk + thij
             fklrii = fklrii + thij * thij
          end if
          xg(j) = real(ic, dp) + thij
          fklr  = fklr + thij * thij
       end do
       fklr = -(fklr - fklrii) - 0.5_dp * fklrii
       call sdet(xg, NSP, dethk)
       a2kt = a2kt + ics * (log(dethk) - trachk - fklr) * qi * qi
    end do
    a2kt = -a2kt * DQ / (12.0_dp * rhot * PI**2)

    sek   = 0.0_dp
    seklj = 0.0_dp
    psum  = 0.0_dp
    ankt  = 0.0_dp
    akcor = 0.0_dp

    ! Sum contributions across interaction pairs
    l = 1
    do i = 1, NSP
       do j = 1, i
          ipk     = real(j / i, dp)
          ic      = 2 - int(ipk)
          fcnr    = real(ic, dp)
          a1kt(l) = a1kt(l) * rho(i) * rho(j) / rhot
          aklr(l) = fcnr * aklr(l) * rho(i) * rho(j) / rhot
          pr(l)   = rho(i) * rho(j) * pr(l) / rhot
          sek     = sek + real(ic, dp) * ekt(l) * rho(i) * rho(j) / rhot
          seklj   = seklj + real(ic, dp) * ektlj(l) * rho(i) * rho(j) / rhot
          ankt    = ankt + a1kt(l) * real(ic, dp)
          psum    = psum + pr(l) * real(ic, dp)
          akcor   = akcor + aklr(l)

          tds(i, j) = 0.0_dp
          do k = 1, NR
             tds(i, j) = tds(i, j) + dcf(k, l) * R(k)
          end do
          tds(i, j) = FOURPI * tds(i, j) * DR
          if (i /= j) tds(j, i) = tds(i, j)
          l = l + 1
       end do
    end do

    ! Concentration fluctuations and Bhatia-Thornton eigenvalues
    xd   = XXc(1)
    xion = sum(XXc(2:NSP))
    if (NSP == 2) then
       cionion = tds(2, 2)
       tds12_eff = tds(1, 2)
    else
       cionion = 0.0_dp
       tds12_eff = 0.0_dp
       do k = 2, NSP
          tds12_eff = tds12_eff + (XXc(k) / xion) * tds(1, k)
          do l = 2, NSP
             cionion = cionion + (XXc(k) * XXc(l) / (xion**2)) * tds(k, l)
          end do
       end do
    end if

    Mrr0 = 1.0_dp - rhot * (xion**2 * cionion + xd**2 * tds(1, 1) + 2.0_dp * xion * xd * tds12_eff)
    Mcc0 = 1.0_dp - rhot * xion * xd * (cionion + tds(1, 1) - 2.0_dp * tds12_eff)
    Mrc0 = rhot * sqrt(xion * xd) * (xd * tds12_eff - xion * cionion - (xd - xion) * tds12_eff)

    lamb1 = (Mrr0 + Mcc0 - sqrt((Mrr0 - Mcc0)**2 + 4.0_dp * Mrc0**2)) / 2.0_dp
    lamb2 = (Mrr0 + Mcc0 + sqrt((Mrr0 - Mcc0)**2 + 4.0_dp * Mrc0**2)) / 2.0_dp
    vrr   = (lamb1 - Mcc0) / sqrt((lamb1 - Mcc0)**2 + Mrc0**2)
    vcc   = Mrc0 / sqrt((lamb1 - Mcc0)**2 + Mrc0**2)

    ! Chemical potentials and thermodynamic properties
    sm1  = 0.0_dp
    sm11 = 0.0_dp
    sm3  = 0.0_dp
    chss = 0.0_dp

    do is = 1, NSP
       ssm1  = 0.0_dp
       chsum = 0.0_dp
       do js = 1, NSP
          nps   = npos(is, js)
          sigx  = sigma(nps)
          sm1   = sm1 + rho(is) * rho(js) * sigx**3
          ssm1  = ssm1 + rho(js) * sigx**3
          chsum = rho(js) * (chm(nps) + ekt(nps)) + chsum
       end do
       sm3  = sm3 + (rho(is) * qq(is)**2)**2
       sm11 = sm11 + rho(is)**2 * (sigma(npos(is, is)))**3
       eik  = exp(-rkd * R(N))

       ! 1. Excess chemical potentials
       chemp(is) = 0.5_dp * eik * (0.5_dp * eik - 1.0_dp) * rkd * qq(is)**2 / tr &
                 + chsum + TWOPI * ssm1 / 3.0_dp
       chss = chss + chemp(is) * rho(is)
    end do

    alr = (2.0_dp * rkd**3 / (8.0_dp * PI * rhot) + PI * sm3 / (rkd * rhot * tr**2)) &
        * exp(-2.0_dp * rkd * R(N))
    asr = PI * (2.0_dp * (sm1 - sm11) + 2.0_dp * sm11) / (3.0_dp * rhot)
    elr = -rkd**3 * exp(-rkd * R(N)) / (8.0_dp * PI * rhot)

    ! 2. Reduced inverse isothermal compressibility
    xc = 1.0_dp - FOURPI * xcsum * DR / (3.0_dp * rhot)

    ! 3. Excess internal energy (Coulomb + LJ)
    sekc = sek + elr
    sek  = sekc + seklj

    ! 4. Excess Helmholtz free energy
    ankt = ankt + a2kt + asr + sekc + akcor + alr

    ! 5. Excess and total entropy
    sex = sek - ankt

    ! Ideal gas mixture entropy: s_id / (N * k_B) = 2.5 - ln(rho) - sum(x_i * ln(x_i))
    sid = 2.5_dp - log(rhot)
    do is = 1, NSP
       if (XXc(is) > 1.0d-12) then
          sid = sid - XXc(is) * log(XXc(is))
       end if
    end do
    stot = sid + sex

    ! 6. Compressibility factor from virial equation
    pres = 1.0_dp + psum + sekc / 3.0_dp

    ! 7. Compressibility factor from free energy and chemical potential
    za = 1.0_dp + chss / rhot - ankt
  end subroutine compute_thermodynamics

end module thermo_mod
