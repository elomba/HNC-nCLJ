!===============================================================================
! Module: hnc_solver_mod
! Description: Core iterative solver for the Ornstein-Zernike equation with HNC
!              closure, combining the Newton-Raphson scheme of Labik, Malijevsky,
!              and Vonka with Broyles' adaptive mixing and the Kinoshita-Harada
!              continuation strategy.
! References:
!   - S. Labik, A. Malijevsky, P. Vonka, Mol. Phys. 56, 709 (1985).
!   - M. Kinoshita, M. Harada, Mol. Phys. 65, 599 (1988).
!   - E. Lomba, Mol. Phys. 68, 87 (1989).
!===============================================================================
module hnc_solver_mod
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
  ! Subroutine: solve_hnc
  ! Executes the coupled Broyles + Newton-Raphson iteration for the HNC closure.
  !=============================================================================
  subroutine solve_hnc(tdiff, pmixt, tmax, amx, diffnr, compj, ntt, kwit, t_start, &
                       total_broyles_iter, final_error)
    real(dp), intent(in)    :: tdiff, pmixt, tmax, amx, diffnr, t_start
    logical,  intent(inout) :: compj
    integer(ip), intent(inout) :: ntt, kwit
    integer(ip), intent(out)   :: total_broyles_iter
    real(dp), intent(out)   :: final_error

    integer(ip) :: iter, newton, i, j, k, ik, l1, l2, in1, in2, it1, it2, i1, i2
    integer(ip) :: kcut
    real(dp) :: rin, suma, tcc, rj, ssum, delnr, timer, current_time
    real(dp) :: aa1, diffm, alfa, dcflow, dcfc, hc, hlow

    rin = 1.0_dp / real(4 * NR, dp)
    iter = 1

    !---------------------------------------------------------------------------
    ! Outer Broyles iteration loop
    !---------------------------------------------------------------------------
    broyles_loop: do while (.true.)

       ! Step 3: Compute direct correlation c(r) and total correlation h(r)
       do j = 1, NIT
          do i = 1, ncore(j) - 1
             dcf(i, j) = -R(i) - sclr(i, j) - s(i, j)
          end do
          do i = ncore(j), N
             dcf(i, j) = R(i) * ( exp( -ulj(i, j) + (s(i, j) + gamm(j) &
                                - clr(i, j) + sclr(i, j)) / R(i) ) - 1.0_dp ) &
                         - sclr(i, j) - s(i, j)
          end do

          h(1:N, j) = (s(1:N, j) + dcf(1:N, j) + sclr(1:N, j)) / R(1:N)

          ! Step 4: Transform c(r) with core discontinuity correction
          dcflow = -R(ncore(j)) - sclr(ncore(j), j) - s(ncore(j), j)
          dcfc = dcf(ncore(j), j)
          dcf(ncore(j), j) = (dcf(ncore(j), j) + dcflow) / 2.0_dp
          call fftsin(dcf(1, j), cr(j), tdcf(1, j))
          dcf(ncore(j), j) = dcfc

          ! Step 5: Functions d(l) for Newton-Raphson convolution
          hc = h(ncore(j), j)
          hlow = -1.0_dp - sclr(ncore(j), j)
          h(ncore(j), j) = (h(ncore(j), j) + hlow) / 2.0_dp
          call fftcos(0.0_dp, h(1, j), d(1, j))

          d(0, j) = 4.0_dp * sum(h(1:N, j))
          h(ncore(j), j) = hc

          ! Initial guess for Newton-Raphson
          tsn(1:NUM, j) = ts(1:NUM, j)
       end do

       !------------------------------------------------------------------------
       ! Step 6: Newton-Raphson inner iterations over NUM points
       !------------------------------------------------------------------------
       newton = 0
       nr_loop: do while (.true.)

          ! 6.1 Check CPU time limit
          call cpu_time(current_time)
          timer = current_time - t_start
          if (timer > tmax) then
             kwit = 1
             write(*, '(A, I4, A, I4)') ' *** Convergence not achieved: time limit exceeded ***'
             return
          end if

          ! 6.2 Compute residuals dels and Jacobian matrix rja
          do i = 1, NUM
             do j = 1, NIT
                suma = 0.0_dp
                do k = 1, NUM
                   ik = abs(i - k)
                   suma = suma + (d(ik, j) - d(i + k, j)) * (tsn(k, j) - ts(k, j))
                end do
                tcc = tdcf(i, j) + suma * rin
                xg(j)  = Q(i) * vki(i, j) - tcc
                tsr(j) = tcc + Q(i) * vk(i, j)
             end do

             call ssinv(xg, NSP)

             do j = 1, NIT
                tsr(j) = xg(j) * Q(i)**2 - tsr(j)
                dels(NUM * (j - 1) + i) = tsr(j) - tsn(i, j)
             end do

             ! Only construct Jacobian if requested (Kinoshita-Harada strategy)
             if (compj) then
                do l1 = 1, NIT
                   do l2 = l1, NIT
                      in1 = ic1(l1, l2)
                      in2 = ic2(l1, l2)
                      rj  = -Q(i)**2 * xg(in1) * xg(in2)
                      if (l1 == l2) rj = rj + 1.0_dp
                      it1 = NUM * (l1 - 1)
                      i1  = it1 + i
                      it2 = NUM * (l2 - 1)
                      do k = 1, NUM
                         i2 = it2 + k
                         ik = abs(i - k)
                         rja(i1, i2) = rj * (d(ik, l2) - d(i + k, l2)) * rin
                         if (l1 /= l2) rja(i + it2, k + it1) = rj * (d(ik, l1) - d(i + k, l1)) * rin
                         if (i1 == i2) rja(i1, i1) = rja(i1, i1) + 1.0_dp
                      end do
                   end do
                end do
             end if
          end do

          ! Invert Jacobian matrix using LAPACK
          if (compj) call invmat(rja, ND, NIT * NUM, rjai, work, ipvt)

          ! Compute Newton correction vector using optimized BLAS DGEMV
          call matvec_mult(rjai, ND, NIT * NUM, dels, delg)

          ssum = sum(delg(1:NIT * NUM)**2)
          delnr = sqrt(ssum)
          newton = newton + 1
          ntt = ntt + 1

          if (delnr > 400.0_dp .and. ntt > 1000) then
             kwit = 1
             write(*, '(A)') ' *** Divergence detected in Newton-Raphson loop ***'
             return
          end if

          ! Check convergence of Newton-Raphson iteration
          if (delnr > diffnr) then
             do j = 1, NIT
                do i = 1, NUM
                   tsn(i, j) = tsn(i, j) + delg(i + NUM * (j - 1))
                end do
             end do
             cycle nr_loop
          else
             exit nr_loop
          end if
       end do nr_loop

       !------------------------------------------------------------------------
       ! Step 7: Direct OZ equation solution in k-space for i = NUM + 1 to N
       !------------------------------------------------------------------------
       aa1 = 0.0_dp
       do kcut = N, NUM + 1, -1
          do j = 1, NIT
             aa1 = max(aa1, abs(tdcf(kcut, j)))
          end do
          if (aa1 > 5.0d-4) exit
       end do

       do i = NUM + 1, N
          do j = 1, NIT
             tcc = tdcf(i, j)
             xg(j)  = Q(i) * vki(i, j) - tcc
             tsr(j) = tcc + Q(i) * vk(i, j)
          end do
          call ssinv(xg, NSP)
          tsn(i, 1:NIT) = xg(1:NIT) * Q(i)**2 - tsr(1:NIT)
       end do

       ! Step 8: Transform back to real-space s(r)
       do j = 1, NIT
          call fftsin(tsn(1, j), cq(j), snew(1, j))
       end do

       ! Step 9: Check convergence in real space
       diffm = 0.0_dp
       do j = 1, NIT
          diff(j) = sqrt(DR * sum((snew(1:N, j) - s(1:N, j))**2))
          diffm = max(diffm, diff(j))
       end do

       ! Step 10: Broyles mixing with adaptive mixing parameter
       do j = 1, NIT
          if (diff(j) > amx) then
             alfa = pmixt
          else
             alfa = 1.0_dp - (1.0_dp - pmixt) * (diff(j) / amx)**2
          end if
          s(1:N, j)  = s(1:N, j)  + alfa * (snew(1:N, j) - s(1:N, j))
          ts(1:N, j) = ts(1:N, j) + alfa * (tsn(1:N, j)  - ts(1:N, j))
       end do

       call cpu_time(current_time)
       timer = current_time - t_start
       write(3, '(3x, i3, 10x, i3, 6x, f9.6, 2x, f8.2)') iter, newton, diffm, timer
       write(*, '(3x, i3, 10x, i3, 6x, f9.6, 2x, f8.2)') iter, newton, diffm, timer

       if (diffm < tdiff) then
          compj = .false.

          ! Recompute Jacobian at final solution to ensure stability
          do i = 1, NUM
             do j = 1, NIT
                suma = 0.0_dp
                do k = 1, NUM
                   ik = abs(i - k)
                   suma = suma + (d(ik, j) - d(i + k, j)) * (tsn(k, j) - ts(k, j))
                end do
                tcc = tdcf(i, j) + suma * rin
                xg(j)  = Q(i) * vki(i, j) - tcc
                tsr(j) = tcc + Q(i) * vk(i, j)
             end do
             call ssinv(xg, NSP)
             do j = 1, NIT
                tsr(j) = xg(j) * Q(i)**2 - tsr(j)
                dels(NUM * (j - 1) + i) = tsr(j) - tsn(i, j)
             end do
             do l1 = 1, NIT
                do l2 = l1, NIT
                   in1 = ic1(l1, l2)
                   in2 = ic2(l1, l2)
                   rj  = -Q(i)**2 * xg(in1) * xg(in2)
                   if (l1 == l2) rj = rj + 1.0_dp
                   it1 = NUM * (l1 - 1)
                   i1  = it1 + i
                   it2 = NUM * (l2 - 1)
                   do k = 1, NUM
                      i2 = it2 + k
                      ik = abs(i - k)
                      rja(i1, i2) = rj * (d(ik, l2) - d(i + k, l2)) * rin
                      if (l1 /= l2) rja(i + it2, k + it1) = rj * (d(ik, l1) - d(i + k, l1)) * rin
                      if (i1 == i2) rja(i1, i1) = rja(i1, i1) + 1.0_dp
                   end do
                end do
             end do
          end do

          call invmat(rja, ND, NIT * NUM, rjai, work, ipvt)
          total_broyles_iter = iter
          final_error = diffm
          exit broyles_loop
       else
          iter = iter + 1
       end if

    end do broyles_loop

  end subroutine solve_hnc

end module hnc_solver_mod
