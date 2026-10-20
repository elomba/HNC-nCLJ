!===============================================================================
! Program: HNCR_nCLJ
! Description: Solution of the Ornstein-Zernike equation with HNC closure for
!              multicomponent (up to 4) charged Lennard-Jones fluids (soft-core
!              ions) using Allnatt's renormalized OZ equation and the Labik-
!              Malijevsky-Vonka Newton-Raphson scheme with Kinoshita-Harada
!              continuation.
!
! Authors:
!   - E. Lomba and J.S. Hoye (Original charged mixture code, Trondheim 1991)
!   - E. Lomba and C. Martin (Fortran 90 port, Madrid 2007)
!   - E. Lomba (Lennard-Jones soft potential version, Ferrol 2020)
!   - Antigravity Refactoring & Optimization (Modular F2008+, BLAS/LAPACK, 2026)
!===============================================================================
program main
  use types_mod
  use params_mod
  use grid_mod
  use fft_mod
  use arrays_mod
  use linalg_mod
  use potentials_mod
  use initial_guess_mod
  use thermo_mod
  use hnc_solver_mod
  implicit none

  character(len=256) :: input_file
  character(len=80)  :: title
  character(len=1)   :: tit
  character(len=10)  :: mu_header

  integer(ip) :: nr_input, num_input
  real(dp)    :: dr_input, tdiff, pmixt, tmax, amx, diffnr, gama
  logical     :: compsi, compj
  integer(ip) :: nrho, ntr, idn, itr, i, j, k, l, lx, inr, nrmax, ni
  integer(ip) :: mx, nx, in1, in2, ncmin, ntt, kwit, total_broyles
  real(dp)    :: rhot, rkx, tr, rkd, xp, xm, dri, tri, rr, final_err
  real(dp)    :: pres, sek, sekc, ankt, sex, stot, s2, za, xc, lamb1, lamb2, vrr, vcc
  real(dp)    :: sx, sq
  real(dp)    :: t_start, t_end, t_step, sksum, scccon

  ! Command-line argument handling
  if (command_argument_count() >= 1) then
     call get_command_argument(1, input_file)
  else
     input_file = 'datosLJ.hncr'
  end if

  write(*, '(/, A)') '=================================================================='
  write(*, '(A)')    '   HNCR-nCLJ: Hypernetted Chain Solver for Charged LJ Mixtures   '
  write(*, '(A, /)') '=================================================================='
  write(*, '(A, A)') ' Reading input file: ', trim(input_file)

  open(9, file=trim(input_file), status='old', action='read')
  read(9, '(A80)') title
  read(9, *) nr_input, dr_input, tdiff, pmixt, tmax, amx
  read(9, *) num_input, diffnr
  read(9, *) NSP

  if (NSP > 10) then
     write(*, '(A, I4)') ' *** Error: Number of species exceeds maximum supported (10): ', NSP
     stop
  end if

  NUM = num_input
  NIT = (NSP * NSP + NSP) / 2
  ND  = NIT * NUM

  ! Initialize grid and FFT modules
  call init_fft(nr_input)
  call init_grid(NR, dr_input)

  ! Allocate dynamic arrays
  call allocate_arrays()
  ulj(:, :) = 0.0_dp
  wlj(:, :) = 0.0_dp

  ! Read species parameters
  do i = 1, NSP
     read(9, *) ssig(i, 1:NSP)
  end do
  do i = 1, NSP
     read(9, *) eeps(i, 1:NSP)
  end do
  read(9, *) qq(1:NSP)
  read(9, *) nrho, ntr
  read(9, *) (densid(j), j = 1, nrho)
  read(9, *) (tk(j), j = 1, ntr)
  read(9, *) compsi

  if (NSP == 2) then
     read(9, *) XXc(1)
  else
     read(9, *) XXc(1:NSP-2)
  end if

  do i = 1, NSP
     read(9, *) rct(i, 1:NSP)
  end do
  read(9, *) gama
  close(9)

  ! Map cutoff radii to symmetric NIT storage
  lx = 1
  do i = 1, NSP
     do j = 1, i
        rcut(lx) = rct(i, j)
        lx = lx + 1
     end do
  end do

  ! Determine remaining mole fractions to guarantee electroneutrality
  if (NSP == 2) then
     XXc(1) = 1.0_dp / (1.0_dp - qq(1) / qq(2))
     XXc(2) = 1.0_dp - XXc(1)
  else
     if (abs(qq(NSP) - qq(NSP-1)) < 1.0d-12) then
        write(*, '(A)') ' *** Error: Charges of species NSP-1 and NSP must differ for electroneutrality ***'
        stop
     end if
     sx = sum(XXc(1:NSP-2))
     sq = sum(XXc(1:NSP-2) * qq(1:NSP-2))
     XXc(NSP)   = (qq(NSP-1) * (sx - 1.0_dp) - sq) / (qq(NSP) - qq(NSP-1))
     XXc(NSP-1) = 1.0_dp - sx - XXc(NSP)
  end if

  ! Initialize correlation functions
  dcf(:, :)   = 0.0_dp
  s(:, :)     = 0.0_dp
  sint(:, :)  = 0.0_dp
  ipg(1:NIT)  = .false.

  ! Set up index mappings and interaction pair arrays
  mx = 1
  do i = 1, NSP
     ipg((i * i + i) / 2) = .true.
     do j = 1, i
        npos(i, j) = mx
        npos(j, i) = mx
        sigma(mx) = ssig(i, j)
        eps(mx)   = eeps(i, j)
        nx = 1
        do k = 1, NSP
           do l = 1, k
              if (k > i) then
                 in2 = (k * k - k) / 2 + i
              else
                 in2 = (i * i - i) / 2 + k
              end if
              if (l > j) then
                 in1 = (l * l - l) / 2 + j
              else
                 in1 = (j * j - j) / 2 + l
              end if
              ic1(mx, nx) = in1
              ic2(mx, nx) = in2
              nx = nx + 1
           end do
        end do
        mx = mx + 1
     end do
  end do

  ! Cutoff indices
  do j = 1, NIT
     ncore(j) = nint(0.6_dp * sigma(j) / DR)
  end do

  ! Open output files
  open(3, file='Salida.hncr',    status='unknown')
  open(8, file='Thermohncr.dat', status='unknown')

  write(3, 99003) trim(title), (j, qq(j), j, sigma((j*j+j)/2), j, XXc(j), j = 1, NSP)
  write(8, 99003) trim(title), (j, qq(j), j, sigma((j*j+j)/2), j, XXc(j), j = 1, NSP)
  write(*, 99003) trim(title), (j, qq(j), j, sigma((j*j+j)/2), j, XXc(j), j = 1, NSP)

  write(3, 99004) tdiff, pmixt, NR, DR, NUM, diffnr
  write(8, 99004) tdiff, pmixt, NR, DR, NUM, diffnr
  write(*, 99004) tdiff, pmixt, NR, DR, NUM, diffnr

  write(8, '("#  Thermodynamics:")')
  write(8, '("# -----------------")')
  write(8, '("# ", 4A8, 8A10)', advance='no') &
       'x(1)', 'rho', 'beta*', 'z', 'u/nkt', 'u_c/nkt', 'a/nkt', 's_ex/nk', 's/nk', 's2/nk', 'z(a)', 'xc'
  do j = 1, NSP
     write(mu_header, '("mu/kt(", I0, ")")') j
     write(8, '(A10)', advance='no') adjustr(mu_header)
  end do
  write(8, '(4A10)') 'lamb1', 'lamb2', 'vrr', 'vcc'
  write(8, '("#", 190("*"))')

  ! Initial solution
  if (compsi) then
     rhot = densid(1)
     rho(1:NSP) = rhot * XXc(1:NSP)
     if (NSP == 2) then
        call msa(s, NR, rho, qq, 1.0_dp / tk(1), sigma(1))
     else
        s(1:NR, 1:NIT) = 0.0_dp
     end if
  else
     open(10, file='solin.hncr', status='old')
     read(10, *) rhot, dri, ni, tri, (gammi(j), j = 1, NIT)
     write(3, 99016)
     write(*, 99016)
     if (abs(DR - dri) > 1.0d-6) then
        write(3, 99017)
        write(*, 99017)
        nrmax = min(ni, NR - 1)
        do i = 1, nrmax
           read(10, *) rr, (sint(i, j), j = 1, NIT)
        end do
        nrmax = ni
        if (R(NR - 1) < rr) nrmax = NR - 1
        do i = 1, nrmax
           inr = int(R(i) / dri)
           do j = 1, NIT
              s(i, j) = sint(inr, j) + (sint(inr + 1, j) - sint(inr, j)) &
                      * (real(i, dp) * DR / dri - real(inr, dp))
           end do
        end do
     else
        do i = 1, min(ni, NR)
           read(10, *) rr, (s(i, j), j = 1, NIT)
        end do
     end if
     close(10)
  end if

  call cpu_time(t_start)
  kwit  = 0
  ntt   = 0
  compj = .true.

  ! Density and temperature loops
  density_loop: do idn = 1, nrho
     rhot = densid(idn)
     rho(1:NSP) = rhot * XXc(1:NSP)

     write(*, '(/, A)') ' Calculating for :'
     do j = 1, NSP
        write(*, '(" rho(", I0, ") = ", F10.6)') j, rho(j)
     end do
     write(*, 99022) PI * sum(rho(1:NSP) * sigma(1:NSP)**3)

     rkx = sum(rho(1:NSP) * qq(1:NSP)**2)

     temperature_loop: do itr = 1, ntr
        tr = 1.0_dp / tk(itr)

        ! Compute LJ potential and virial
        call compute_lj_potential(tr)

        ! Compute Debye screening length and long-range corrections
        call compute_screening_and_lr(tr, rkx, gama, rkd, xp, xm)

        write(3, 99012) 1.0_dp / tr, rhot
        write(*, 99012) 1.0_dp / tr, rhot
        write(3, 99001)
        write(*, 99001)

        ncmin = NR
        do j = 1, NIT
           ncmin = min(ncmin, ncore(j))
           call fftsin(s(1, j), cr(j), ts(1, j))
        end do
        ncmin = ncmin - 10

        ! Solve the HNC integral equations
        call solve_hnc(tdiff, pmixt, tmax, amx, diffnr, compj, ntt, kwit, t_start, &
                       total_broyles, final_err)

        if (kwit == 1) then
           write(*, '(A)') ' *** Aborting run due to error/timeout in solver ***'
           exit density_loop
        end if

        ! Calculate thermodynamics
        call compute_thermodynamics(tr, rhot, rkd, pres, sek, sekc, ankt, sex, stot, s2, za, xc, &
                                    lamb1, lamb2, vrr, vcc)

        call cpu_time(t_step)
        write(3, 99009) ntt, total_broyles, t_step - t_start
        write(8, '(2x, 4(f8.4), *(f10.4))') &
             XXc(1), rhot, 1.0_dp / tr, pres, &
             sek, sekc, ankt, sex, stot, s2, za, xc, &
             (chemp(j), j = 1, NSP), lamb1, lamb2, vrr, vcc

     end do temperature_loop
  end do density_loop

  call cpu_time(t_end)

  ! Store short-range s(r) solution for restart (lsol.hncr)
  open(2, file='lsol.hncr', status='unknown')
  write(2, 99031) rhot, DR, N, tr, (gamm(j), j = 1, NIT)
  do i = 1, N
     write(2, 99021) R(i), (s(i, j), j = 1, NIT)
  end do
  close(2)

  ! Output radial distribution functions g_ij(r) (RDFhncr.dat)
  open(4, file='RDFhncr.dat', status='unknown')
  do i = ncmin, N
     write(4, 99020) R(i), (h(i, j) + 1.0_dp, j = 1, NIT)
  end do
  close(4)

  ! Output structure factors S_ij(q) and Bhatia-Thornton S_NN, S_cc (skhnc.dat)
  open(99, file='skhnc.dat', status='unknown')
  do i = 1, N
     mx = 1
     sksum = 0.0_dp
     do k = 1, NSP
        do l = 1, k
           if (k == l) then
              sk(i, mx) = rho(l) * (1.0_dp + th(i, mx) / Q(i)) / rhot
              sksum = sksum + sk(i, mx)
           else
              sk(i, mx) = sqrt(rho(k) * rho(l)) * (th(i, mx) / Q(i)) / rhot
              sksum = sksum + 2.0_dp * sk(i, mx)
           end if
           mx = mx + 1
        end do
     end do
     call skcc(sk(i, 1:NIT), rho(1:NSP), rhot, 3, NSP, scccon, NIT)
     write(99, 99020) Q(i), sk(i, 1:NIT), sksum, scccon
  end do
  close(99)

  close(3)
  close(8)

  write(*, '(/, A, F9.3, A)') ' >>> Total computation CPU time: ', t_end - t_start, ' seconds.'
  write(*, '(A, /)')          ' >>> Output files generated: Salida.hncr, Thermohncr.dat, RDFhncr.dat, skhnc.dat, lsol.hncr'

  ! Clean up memory
  call deallocate_arrays()
  call cleanup_grid()
  call cleanup_fft()

  ! Format strings
99001 format(/2x, 'broyles it.   n-r it      error      time'/1x, 50('='))
99003 format(/, '#      program HNCR_nCLJ: '/"#   ", a, /"#  ", 60('*') &
             / *( '#   qq(', i1, ') = ', f9.4, ' sigma(', i1, ') = ', f7.4, &
                ' Xc(', i1, ') = ', f8.5 /: ))
99004 format('# tdiff = ', f10.7, ' pmixt = ', f6.3, ' nr = ', i5, ' dr = ', f7.4, &
             ' num = ', i3, ' diffnr = ', f10.7//)
99009 format(/1x, ' * end of computation after ', i4, ' n-r iter.', &
             ' and ', i3, ' broyles iter.'/4x, 'cpu time :', f7.2, ' sec.'/)
99012 format(/2x, 'computing for beta* = ', f9.4, 5x, 'rho (tot) = ', f8.5/1x, 60('-'))
99016 format(/' *** Reading initial solution from disk ***.'/)
99017 format(/'  ** Interpolating S(I,J) to match input dr **'/)
99020 format(1x, f9.5, 3x, *(1x, f11.5))
99021 format(f10.5, 4x, *(f12.5, 1x))
99022 format(" eta    = ", f10.6/)
99031 format(2(f8.5, 1x), i5, 1x, *(f14.9, 1x))

end program main
