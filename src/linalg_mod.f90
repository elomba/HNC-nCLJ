!===============================================================================
! Module: linalg_mod
! Description: Linear algebra routines: LAPACK/BLAS wrappers for matrix inversion
!              and matrix-vector multiplication, plus specialized routines for
!              symmetric matrices (SSINV, SDET, SPROD).
!===============================================================================
module linalg_mod
  use types_mod
  implicit none
  public

  ! External LAPACK and BLAS routines
  interface
     subroutine dgetrf(m, n, a, lda, ipiv, info)
       import :: dp, ip
       integer(ip), intent(in)    :: m, n, lda
       real(dp), dimension(lda,*) :: a
       integer(ip), dimension(*)  :: ipiv
       integer(ip), intent(out)   :: info
     end subroutine dgetrf

     subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
       import :: dp, ip
       integer(ip), intent(in)    :: n, lda, lwork
       real(dp), dimension(lda,*) :: a
       integer(ip), dimension(*)  :: ipiv
       real(dp), dimension(*)     :: work
       integer(ip), intent(out)   :: info
     end subroutine dgetri

     subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
       import :: dp, ip
       character(len=1), intent(in) :: trans
       integer(ip), intent(in)      :: m, n, lda, incx, incy
       real(dp), intent(in)         :: alpha, beta
       real(dp), dimension(lda,*)   :: a
       real(dp), dimension(*)       :: x, y
     end subroutine dgemv
  end interface

contains

  !=============================================================================
  ! Subroutine: invmat
  ! Computes the inverse of a square general matrix using LAPACK DGETRF + DGETRI.
  !=============================================================================
  subroutine invmat(a, lda, n, b, work, ipvt)
    integer(ip), intent(in)  :: lda, n
    real(dp), dimension(lda, n), intent(in)  :: a
    real(dp), dimension(lda, n), intent(out) :: b
    real(dp), dimension(n),      intent(out) :: work
    integer(ip), dimension(n),   intent(out) :: ipvt

    integer(ip) :: info

    b(1:n, 1:n) = a(1:n, 1:n)

    call dgetrf(n, n, b, lda, ipvt, info)
    if (info /= 0) then
       write(*, '(A, I6)') ' *** Error from invmat: DGETRF failed, info = ', info
       stop ' *** Matrix is singular in invmat ***'
    end if

    call dgetri(n, b, lda, ipvt, work, n, info)
    if (info /= 0) then
       write(*, '(A, I6)') ' *** Error from invmat: DGETRI failed, info = ', info
       stop ' *** Matrix inversion failed in invmat ***'
    end if
  end subroutine invmat

  !=============================================================================
  ! Subroutine: matvec_mult
  ! Optimized matrix-vector multiplication y = A * x using BLAS DGEMV.
  !=============================================================================
  subroutine matvec_mult(a, lda, n, x, y)
    integer(ip), intent(in) :: lda, n
    real(dp), dimension(lda, n), intent(in)  :: a
    real(dp), dimension(n),      intent(in)  :: x
    real(dp), dimension(n),      intent(out) :: y

    call dgemv('N', n, n, 1.0_dp, a, lda, x, 1, 0.0_dp, y, 1)
  end subroutine matvec_mult

  !=============================================================================
  ! Subroutine: sdet
  ! Returns the determinant of a symmetric N x N matrix stored in packed mode.
  ! Programmed by F. Lado.
  !=============================================================================
  subroutine sdet(a, n, det)
    integer(ip), intent(in) :: n
    real(dp), dimension(*), intent(inout) :: a
    real(dp), intent(out) :: det

    integer(ip) :: i, j, j1, j2, jm, jp, k, k1
    real(dp) :: a1, sum_val

    if (a(1) >= 0.0_dp) then
       det = a(1)
       if (n == 1) return

       a(1) = sqrt(a(1))
       k1 = 1
       do k = 2, n
          a(k1 + 1) = a(k1 + 1) / a(1)
          k1 = k1 + k
       end do

       j1 = 1
       do j = 2, n
          sum_val = 0.0_dp
          jm = j - 1
          do i = 1, jm
             sum_val = sum_val + a(j1 + i)**2
          end do
          j2 = j1 + j
          a1 = a(j2) - sum_val
          if (a1 < 0.0_dp) goto 100

          det = det * a1
          a(j2) = sqrt(a1)

          if (j /= n) then
             jp = j + 1
             k1 = j2
             do k = jp, n
                sum_val = 0.0_dp
                do i = 1, jm
                   sum_val = sum_val + a(j1 + i) * a(k1 + i)
                end do
                a(k1 + j) = (a(k1 + j) - sum_val) / a(j2)
                k1 = k1 + k
             end do
          end if
          j1 = j2
       end do
       return
    end if

100 write(*,*) ' *** Error in routine sdet: non-positive definite matrix ***'
    det = 0.0_dp
  end subroutine sdet

  !=============================================================================
  ! Subroutine: sprod
  ! Returns matrix product C = A * B for symmetric N x N matrices in packed mode.
  ! Adapted from IBM Scientific Subroutine Package by F. Lado.
  !=============================================================================
  subroutine sprod(a, b, c, n)
    integer(ip), intent(in) :: n
    real(dp), dimension(*), intent(in)  :: a, b
    real(dp), dimension(*), intent(out) :: c

    integer(ip) :: i, i1, is, j, j1, js, k, ks
    real(dp) :: sum_val

    js = 0
    do j = 1, n
       is = 0
       do i = 1, j
          sum_val = 0.0_dp
          do k = 1, i
             sum_val = sum_val + a(is + k) * b(js + k)
          end do
          ks = is + i
          if (i /= j) then
             i1 = i + 1
             do k = i1, j
                sum_val = sum_val + a(ks + i) * b(js + k)
                ks = ks + k
             end do
          end if
          if (j /= n) then
             j1 = j + 1
             do k = j1, n
                sum_val = sum_val + a(ks + i) * b(ks + j)
                ks = ks + k
             end do
          end if
          c(js + i) = sum_val
          is = is + i
       end do
       js = js + j
    end do
  end subroutine sprod

  !=============================================================================
  ! Subroutine: ssinv
  ! Inversion of a symmetric matrix in packed storage mode.
  ! Source code adapted from CERN Computer Library.
  !=============================================================================
  subroutine ssinv(ai, n)
    integer(ip), intent(in) :: n
    real(dp), dimension(*), intent(inout) :: ai

    integer(ip) :: i, j, k, km1, kp1, l
    integer(ip), dimension(100) :: r_stat
    real(dp), dimension(6, 6)   :: a_mat
    real(dp), dimension(100)    :: p, q
    real(dp) :: big, pivot1, test
    real(dp), parameter :: epsiln = 1.0_dp - 20 ! 1.0d-20

    if (n < 1) then
       return
    else if (n == 1) then
       ai(1) = 1.0_dp / ai(1)
       return
    else
       l = 1
       do i = 1, n
          do j = 1, i
             a_mat(j, i) = ai(l)
             l = l + 1
          end do
       end do

       r_stat(1:n) = 1

       do i = 1, n
          k = 0
          big = 0.0_dp
          do j = 1, n
             test = abs(a_mat(j, j))
             if (test > big) then
                if (r_stat(j) < 0) goto 100
                if (r_stat(j) /= 0) then
                   big = test
                   k = j
                end if
             end if
          end do

          if (k <= 0) goto 100

          if (i == 1) pivot1 = a_mat(k, k)
          if (abs(a_mat(k, k) / pivot1) < 1.0d-20) goto 100

          r_stat(k) = 0
          q(k) = 1.0_dp / a_mat(k, k)
          p(k) = 1.0_dp
          a_mat(k, k) = 0.0_dp
          kp1 = k + 1
          km1 = k - 1
          if (km1 < 0) goto 100
          if (km1 /= 0) then
             do j = 1, km1
                p(j) = a_mat(j, k)
                q(j) = a_mat(j, k) * q(k)
                if (r_stat(j) < 0) goto 100
                if (r_stat(j) /= 0) q(j) = -q(j)
                a_mat(j, k) = 0.0_dp
             end do
          end if
          if (k < n) then
             do j = kp1, n
                p(j) = a_mat(k, j)
                q(j) = -a_mat(k, j) * q(k)
                if (r_stat(j) < 0) goto 100
                if (r_stat(j) == 0) p(j) = -p(j)
                a_mat(k, j) = 0.0_dp
             end do
          else if (k /= n) then
             goto 100
          end if

          do j = 1, n
             do k = j, n
                a_mat(j, k) = a_mat(j, k) + p(j) * q(k)
             end do
          end do
       end do

       l = 1
       do j = 1, n
          do k = 1, j
             ai(l) = a_mat(k, j)
             l = l + 1
          end do
       end do
       return
    end if

100 write(*,*) ' *** Failure return from subroutine ssinv ***'
  end subroutine ssinv

end module linalg_mod
