!===============================================================================
! Module: arrays_mod
! Description: Dynamic arrays for correlation functions, potentials,
!              screening factors, Newton-Raphson system, and thermodynamics.
!===============================================================================
module arrays_mod
  use types_mod
  use params_mod
  use grid_mod, only: DR
  use fft_mod,  only: NR
  implicit none
  public

  ! Index and bookkeeping matrices
  integer(ip), dimension(:,:), allocatable :: ic1, ic2  ! Pair interaction index mappings
  integer(ip), dimension(:,:), allocatable :: npos      ! Position mapping npos(i,j)
  integer(ip), dimension(:),   allocatable :: ipvt      ! LAPACK pivoting vector
  integer(ip), dimension(:),   allocatable :: ncore     ! Core radius index
  integer(ip), dimension(:),   allocatable :: imin      ! Minimum index

  ! Correlation and screening functions on radial/reciprocal grid (NR, NIT)
  real(dp), dimension(:,:), allocatable :: dcf   ! Direct correlation function c_ij(r)
  real(dp), dimension(:,:), allocatable :: tdcf  ! Fourier transformed c_ij(k)
  real(dp), dimension(:,:), allocatable :: ts    ! Series function s_ij(k) = h - c
  real(dp), dimension(:,:), allocatable :: h     ! Total correlation function h_ij(r)
  real(dp), dimension(:,:), allocatable :: tsn   ! New series function guess in k-space
  real(dp), dimension(:,:), allocatable :: s     ! Series function in r-space: r * (h - c)
  real(dp), dimension(:,:), allocatable :: snew  ! Updated series function in r-space
  real(dp), dimension(:,:), allocatable :: clr   ! Long-range component of c_ij(r)
  real(dp), dimension(:,:), allocatable :: sclr  ! Long-range component of s_ij(r)
  real(dp), dimension(:,:), allocatable :: vk    ! Chain-sum hypervertex matrix v(k)
  real(dp), dimension(:,:), allocatable :: vki   ! Inverse chain-sum matrix [v(k)]^(-1)
  real(dp), dimension(:,:), allocatable :: th    ! Fourier transform of total correlation h
  real(dp), dimension(:,:), allocatable :: sint  ! Interpolated initial estimate
  real(dp), dimension(:,:), allocatable :: d     ! d(l) functions for Newton-Raphson convolution

  ! Pair potential and virial matrices
  real(dp), dimension(:,:), allocatable :: ulj   ! Truncated & shifted LJ potential
  real(dp), dimension(:,:), allocatable :: wlj   ! Virial -r * du_LJ/dr
  real(dp), dimension(:,:), allocatable :: sk    ! Structure factor matrix S_ij(k)
  real(dp), dimension(:,:), allocatable :: tds   ! Integral of direct correlation function

  ! Newton-Raphson Jacobian system
  real(dp), dimension(:,:), allocatable :: rja   ! Jacobian matrix (ND, ND)
  real(dp), dimension(:,:), allocatable :: rjai  ! Inverse Jacobian matrix (ND, ND)
  real(dp), dimension(:),   allocatable :: dels  ! Residual vector (ND)
  real(dp), dimension(:),   allocatable :: delg  ! Newton correction vector (ND)
  real(dp), dimension(:),   allocatable :: work  ! LAPACK workspace (ND)

  ! Interaction-dependent parameters (NIT)
  real(dp), dimension(:), allocatable :: sigma   ! Contact distances / collision diameters
  real(dp), dimension(:), allocatable :: eps     ! Lennard-Jones well depths
  real(dp), dimension(:), allocatable :: gamm    ! Coulomb coupling parameters
  real(dp), dimension(:), allocatable :: cr      ! Forward FFT density scaling factor
  real(dp), dimension(:), allocatable :: cq      ! Backward FFT density scaling factor
  real(dp), dimension(:), allocatable :: tcn, xg, tc2, tsr
  real(dp), dimension(:), allocatable :: ekt     ! Coulomb energy contributions
  real(dp), dimension(:), allocatable :: ektlj   ! LJ energy contributions
  real(dp), dimension(:), allocatable :: pr      ! Virial pressure contributions
  real(dp), dimension(:), allocatable :: tdcf0, diff, gammi, tmp, chm, a1kt, aklr
  real(dp), dimension(:), allocatable :: ucut    ! Potential cutoff shifts
  real(dp), dimension(:), allocatable :: rcut    ! Interaction cutoff radii

  ! Species-dependent parameters (NSP)
  real(dp), dimension(:), allocatable :: qq      ! Ionic charges in elementary units
  real(dp), dimension(:), allocatable :: chemp   ! Excess chemical potentials
  real(dp), dimension(:), allocatable :: XXc     ! Mole fractions
  real(dp), dimension(:), allocatable :: rho     ! Partial number densities
  real(dp), dimension(:), allocatable :: rhoi    ! Initial reference densities

  ! Input full species-pair matrices (NSP, NSP)
  real(dp), dimension(:,:), allocatable :: rct   ! Input cutoffs matrix
  real(dp), dimension(:,:), allocatable :: ssig  ! Input sigma matrix
  real(dp), dimension(:,:), allocatable :: eeps  ! Input epsilon matrix

  ! Density and temperature series (NTEMP)
  real(dp), dimension(NTEMP) :: densid, tk
  logical,  dimension(:),   allocatable :: ipg

contains

  subroutine allocate_arrays()
    implicit none

    allocate(rcut(NIT), ucut(NIT))
    allocate(rct(NSP, NSP), eeps(NSP, NSP), ssig(NSP, NSP))

    allocate(tds(NSP, NSP))
    allocate(dcf(NR, NIT))
    allocate(ulj(NR, NIT))
    allocate(wlj(NR, NIT))
    allocate(tdcf(NR, NIT))
    allocate(sk(NR, NIT))
    allocate(ts(NR, NIT))
    allocate(h(NR, NIT))
    allocate(tsn(NR, NIT))
    allocate(s(NR, NIT))
    allocate(snew(NR, NIT))
    allocate(clr(NR, NIT))
    allocate(sclr(NR, NIT))
    allocate(vk(NR, NIT))
    allocate(vki(NR, NIT))
    allocate(th(NR, NIT))
    allocate(sint(NR, NIT))
    allocate(d(0:NR, NIT))

    allocate(rja(ND, ND))
    allocate(rjai(ND, ND))
    allocate(dels(ND))
    allocate(delg(ND))
    allocate(work(ND))

    allocate(gamm(NIT))
    allocate(cr(NIT))
    allocate(cq(NIT))
    allocate(tcn(NIT))
    allocate(xg(NIT))
    allocate(tc2(NIT))
    allocate(tsr(NIT))
    allocate(ekt(NIT))
    allocate(ektlj(NIT))
    allocate(sigma(NIT))
    allocate(eps(NIT))

    allocate(qq(NSP))
    allocate(XXc(NSP))
    allocate(rhoi(NSP))
    allocate(rho(NSP))
    allocate(pr(NIT))
    allocate(tdcf0(NIT))
    allocate(diff(NIT))
    allocate(gammi(NIT))
    allocate(tmp(NIT))
    allocate(chm(NIT))
    allocate(a1kt(NIT))
    allocate(aklr(NIT))
    allocate(chemp(NSP))

    allocate(ic1(NIT, NIT))
    allocate(ic2(NIT, NIT))
    allocate(ipvt(ND))
    allocate(ncore(NIT))
    allocate(imin(NIT))
    allocate(ipg(NIT))
    allocate(npos(NSP, NSP))
  end subroutine allocate_arrays

  subroutine deallocate_arrays()
    implicit none

    if (allocated(rcut))   deallocate(rcut)
    if (allocated(ucut))   deallocate(ucut)
    if (allocated(rct))    deallocate(rct)
    if (allocated(eeps))   deallocate(eeps)
    if (allocated(ssig))   deallocate(ssig)
    if (allocated(tds))    deallocate(tds)
    if (allocated(dcf))    deallocate(dcf)
    if (allocated(ulj))    deallocate(ulj)
    if (allocated(wlj))    deallocate(wlj)
    if (allocated(tdcf))   deallocate(tdcf)
    if (allocated(sk))     deallocate(sk)
    if (allocated(ts))     deallocate(ts)
    if (allocated(h))      deallocate(h)
    if (allocated(tsn))    deallocate(tsn)
    if (allocated(s))      deallocate(s)
    if (allocated(snew))   deallocate(snew)
    if (allocated(clr))    deallocate(clr)
    if (allocated(sclr))   deallocate(sclr)
    if (allocated(vk))     deallocate(vk)
    if (allocated(vki))    deallocate(vki)
    if (allocated(th))     deallocate(th)
    if (allocated(sint))   deallocate(sint)
    if (allocated(d))      deallocate(d)
    if (allocated(rja))    deallocate(rja)
    if (allocated(rjai))   deallocate(rjai)
    if (allocated(dels))   deallocate(dels)
    if (allocated(delg))   deallocate(delg)
    if (allocated(work))   deallocate(work)
    if (allocated(gamm))   deallocate(gamm)
    if (allocated(cr))     deallocate(cr)
    if (allocated(cq))     deallocate(cq)
    if (allocated(tcn))    deallocate(tcn)
    if (allocated(xg))     deallocate(xg)
    if (allocated(tc2))    deallocate(tc2)
    if (allocated(tsr))    deallocate(tsr)
    if (allocated(ekt))    deallocate(ekt)
    if (allocated(ektlj))  deallocate(ektlj)
    if (allocated(sigma))  deallocate(sigma)
    if (allocated(eps))    deallocate(eps)
    if (allocated(qq))     deallocate(qq)
    if (allocated(XXc))    deallocate(XXc)
    if (allocated(rhoi))   deallocate(rhoi)
    if (allocated(rho))    deallocate(rho)
    if (allocated(pr))     deallocate(pr)
    if (allocated(tdcf0))  deallocate(tdcf0)
    if (allocated(diff))   deallocate(diff)
    if (allocated(gammi))  deallocate(gammi)
    if (allocated(tmp))    deallocate(tmp)
    if (allocated(chm))    deallocate(chm)
    if (allocated(a1kt))   deallocate(a1kt)
    if (allocated(aklr))   deallocate(aklr)
    if (allocated(chemp))  deallocate(chemp)
    if (allocated(ic1))    deallocate(ic1)
    if (allocated(ic2))    deallocate(ic2)
    if (allocated(ipvt))   deallocate(ipvt)
    if (allocated(ncore))  deallocate(ncore)
    if (allocated(imin))   deallocate(imin)
    if (allocated(ipg))    deallocate(ipg)
    if (allocated(npos))   deallocate(npos)
  end subroutine deallocate_arrays

end module arrays_mod
