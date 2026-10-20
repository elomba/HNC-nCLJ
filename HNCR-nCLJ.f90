!             Program HNCR
!--------------------------------------------------------------------
!
!
!    Program to solve the OZ equation with a given closure
!    (hnc as currently programmed) for mixtures of charged LJ
!    fluids, using Allnat's renomalized OZ equation
!    (A.R. Allnat, Mol.Phys. 8, 533 (1964). the numerical iteration
!    scheme is based on a Newton-Raphson like procedure due
!    S. Labik, A. Malijevsky and P. Vonka, Mol.Phys. 56, 709 (1985)
!    adapted to fluid mixtures and incorporating the
!    iteration strategy for multiple calculations proposed by
!    M. Kinoshita and M. Harada, Mol.Phys.,65,599 (1988) following
!    the procedure suggested by E.Lomba for simple fluids, Mol.Phys.
!    68, 87 (1989).
!
!    Note :
!
!    This is a 4-component version. For multiple components
!    the program input must be modified and a new subroutine
!    for the computation of the initial solution is required.
!
!    Note about program units:
!
!    This program uses reduced units of length, it is adviceable
!    to take the largest ionic diameter as lenght unit. The inverse
!    temperature parameter must be then given in the corresponding
!    reduced units in accordance with the choice of unit length.
!
!
!    Programmed by E. Lomba and J.S. Hoye
!
!
!                  Trondheim 19.5.1991
!
!    F90 version: Madrid, 11.09.2007 (E. Lomba and C. Martin)
!
!    Lennard-Jones version: Ferrol, 05.2020 (E. Lomba)
!
!    Libraries required : this program uses subroutines dgetri and dgetrf from 
!                         LAPACK and BLAS computer libraries (Atlas
!
!                         or Intel's MathKernel).
!    For questions related to the program contact E.Lomba
!    e-mail: enrique.lomba@csic.es
!
!----------------------------------------------------------------
!
! ************************************************************
! *                                                          *
! *   input:                                                 *
! *                                                          *
! *     card #1:                                             *
! *                                                          *
! *       title :  character variable describing the comput- *
! *                ation                                     *
! *                                                          *
! *     card #2:                                             *
! *                                                          *
! *          nr :  point number in fourier transform         *
! *          dr :  grid size                                 *
! *       tdiff :  final error                               *
! *       pmixt :  parameter of mixture in broyles's method  *
! *        tmax :  maximum allowed time for computation (sec)*
! *         amx :  maximum error to start the change in pmixt*
! *                                                          *
! *     card #3:                                             *
! *                                                          *
! *         num :  number of points in n-r iterations        *
! *                (num < or equal numx)		     *
! *      diffnr :  error in n-r iterations                   *
! *                                                          *
! *     card #4:                                             *
! *         nsp : number of species (<=4)                    *
! *                                                          *
! *     card #5:                                             *
! *                                                          *
! *       sigma(j) :   atomic diameter of species  j=1,nsp   *
! *                    It is adviceable to set one of the    *
! *                    diameters to unity (typically the     *
! *                    largest one, sigma)		     *
! *                                                          *
! *     card #6:                                             *
! *                                                          *
! *          qq(j) :   charge in e. units of species         *
! *                    j (j=1,nsp)                           *
! *                                                          *
! *     card #7:                                             *
! *                                                          *
! *     nrho,ntr   :  number of densities and inverse temp.  *
! *                                                          *
! *     card #8:                                             *
! *                                                          *
! *     densid(j) : total density (sum of all ionic species) *
! *                 (corresponds to (N/V)*sigma**3           *
! *                 (j = 1,nrho)                             *
! *                                                          *
! *     card #8:                                             *
! *                                                          *
! *     tk(j)     : inverse temperature (e**2/kt*sigma*eps)  *
! *                 eps is the solvent dielectric constant   *
! *                 (j = 1, ntr)                             *
! *                                                          *
! *     card #10:                                             *
! *                                                          *
! *         compsi :  if .t. an MSA initial solution is sup- *
! *                   plied, otherwise is read from disk     *
! *                                                          *
! *     card #11:                                            *
! *                                                          *
! *         XXc    :  molar fraction (comp)                        *
! *                                                          *
! ************************************************************
! *                                                          *
! *  input files :                                           *
! *               hncrin.hncr : contains input data          *
! *                solin.hncr : contains initial short range *
! *                             component of r*s(r)          *
! *                             (required if compsi=.f.)     *
! *                                                          *
! *                                                          *
! *  output files :                                          *
! *            Salida.hncr : Main output file                *
! *            Thermohnc.dat:Thermodynamics of the systems   *
! *            RDFhncr.dat : Ion-Ion correlation functions   *
! *            lsol.hncr   : short range component of r*s(r) *
! *                                                          *
! ************************************************************
! *                                                          *
! *  Machine dependent routines : second() CPU time in secs. *
! *                           See function second in         *
! *                           this code for various options  *
! *                                                          *
! ************************************************************

!   Module section for intercommunication of routines

Module Params

  !  Specific dimensions are set here

  Integer (kind=4) :: NSP , NUM , NIT , ND 
  Integer (kind=4) , Parameter :: NTEMP = 2000
  Real (kind=8) :: Pi
End Module Params

Module Funciones

  Use Params

  Integer (kind=4) , Dimension (:,:) , Allocatable  :: ic1 , ic2
  Integer (kind=4) , Dimension (:,:) , Allocatable :: npos
  Integer (kind=4) , Dimension (:) , Allocatable:: ipvt
  Integer (kind=4) , Dimension (:) , Allocatable :: ncore , imin 
  Real (kind=8), Dimension (:,:) , Allocatable :: dcf , tdcf , ts , sk, &
       &          h , tsn , s , snew , clr , sclr , &
       &          vk , vki , th, ulj, wlj, tds
  Real (kind=8) , Dimension (:,:) , Allocatable :: sint , d
  Real (kind=8) , Dimension (:,:) , Allocatable :: rja , rjai
  Real (kind=8) , Dimension (:) , Allocatable :: dels , delg , work
  Real (kind=8) , Dimension (:) , Allocatable  :: sigma , gamm ,  &
       &   cr , cq , tcn , xg , tc2 , tsr ,  &
       &   ekt , ektlj, pr , tdcf0 , diff  , gammi , tmp , &
       &   chm  , a1kt , aklr, eps
  Real (kind=8) , Dimension (:) , Allocatable :: qq , chemp , XXc , rho , rhoi
  Real (kind=8) , Dimension (NTEMP) :: densid , tk
  Logical , Dimension (:) , Allocatable:: ipg
  Real (kind=8) :: gama 
End Module Funciones

Module RQst
  Use Params

  !     this module transfers variables to pyhs and msa

  Real (kind=8), Dimension(:) , Allocatable :: R , Q
  Real (kind=8) :: DQ , DR , CQ1 , CR1
End Module RQst

Module FFT

  !     this module transfers variables to fftsin and fftcos and pyhs

  Use Params

  Real (kind=8), Dimension(:) , Allocatable :: SINe , COSe , FR , FI
  Integer (kind=4) :: N , NR , NU
End Module FFT






!      Main program

Program HNCR
  Use Params
  Use RQst
  Use FFT
  Use Funciones
  Implicit None

  Integer (kind=4) :: i , i1 , i2 , ic , ics , id1 , id2 , &
       &        idn , ik , in1 , in2 , inr , ipk , is , it1
  Integer (kind=4) :: it2 , iter , itr , j , js , k , kcut , kwit , l , l1 , &
       &        l2 , mx ,  ncj , ncmin , newton , ni , nps  , &
       &        nrho , nrmax  , ntr , ntt ,  nx
  Real (kind=8) :: a1 , a1sum , a2kt , aa1 , akcor , &
       &                 alfa , alr , alrsum , amx , ankt , asr ,  &
       &                 chss , chsum
  Real (kind=8) :: csv , cv , dcfc , &
       &                 dcflow , delnr , dethk , &
       &                 diffm , diffnr , dnt , dnt2 
  Real (kind=8) :: dri , eik , elr , esum , eljsum, fcnr , fklr , &
       &                 fklrii , hc , hij , &
       &                 hij2 , hlow , pmixt
  Real (kind=8) :: pres , psum , qi , &
       &                 rhot , ri , rin , rj , rkd , rkdx2 , &
       &                 rkx , rms1 , rms2
  Real (kind=8) :: rmusum , rr , sek , sekc, seklj, sigx , &
       &                 sij  , sm1 , sm11 , sm3 , sncj , &
       &                 ssm1 , ssum , suma , tcc
  Real (kind=8) :: tdiff , tfin , thij , &
       &                 thlr , tim , timer , tmax , tr , &
       &                 trachk , tri , tt
  Real (kind=8) :: xc , xcsum , xe , xm , &
       &                 xm2 , xp , xp2 , za

  ! ************************************************************
  ! *							       *
  ! *   Most important program parameters :		       *
  ! *		(to be set in module Params)                   *
  ! *   npmx = maximum no. of integration points	       *
  ! *   numx = maximum no. of points in NR iterations	       *
  ! *   nsp  = number of ionic species			       *
  ! *   nit  = number of different interactions	   	       *
  ! *							       *
  ! *   Note: npmx must be the same in hncr,pyhs,msa,fftcos    *
  ! *	    and fftsin					       *
  ! ************************************************************


  !     Correlation functions (and related quantities) are stored in
  !     symmetric mode, i.e.
  !
  !            A(I,J) = A(J,I) ---> A((J*J-J)/2+I) (I < or = J)


  !   Save some storage (not too much !)
  !
  Real (kind=8) :: ep4, sig,  r6, r12, sksum, scccon, det, scc0, xd,&
       & xion, cionion, Mrr0, Mcc0, Mrc0, lamb1, lamb2, vcc, vrr
  real (kind=8), dimension(:), allocatable :: ucut, rcut
  real (kind=8), dimension(:,:), allocatable :: rct, ssig, eeps
  Real (kind=4) :: SECOND
  integer :: info, lx
  Logical ::  compsi , compj , rgrid
  Character (80) :: title
  Character (1) ::  tit

  !*******************************************************************
  !                                                                  *
  !                    Program Initialization                        *
  !                                                                  *
  !*******************************************************************

  Pi = Acos(-1.00d0)
  !     Initialize time counter

  tt = SECOND()
  kwit       = 0
  ntt        = 0
  compj      = .True.
  rgrid      = .false.

  Open (3,FILE='Salida.hncr',STATUS='unknown')
  Open (8,FILE='Thermohncr.dat',STATUS='unknown')
  Open (9,FILE='datosLJ.hncr',STATUS='old')

  !     read input data

  Read (9,99000) title
  Read (9,*) NR , DR , tdiff , pmixt , tmax , amx
  Read (9,*) num , diffnr
  Read (9,*) NSP

  NIT = (NSP*NSP + NSP)/2
  allocate(rcut(nit), ucut(nit), rct(nsp,nsp),eeps(nsp,nsp),ssig(nsp,nsp))
  ND  = NIT*NUM

  !    distribucion de las matrices (allocate)

  Call Distribucion
  ulj(:,:) = 0.0d0
  wlj(:,:) = 0.0d0

!  Read (9,*) (sigma((j*j+j)/2),j=1,NSP)
!  Read (9,*) (eps((j*j+j)/2),j=1,NSP)
  do i = 1, nsp
     read(9,*) ssig(i,1:nsp)
  enddo
  do i = 1, nsp
     read(9,*) eeps(i,1:nsp)
  enddo
  Read (9,*) qq(1:NSP)
  Read (9,*) nrho , ntr
  Read (9,*) (densid(j),j=1,nrho)
  Read (9,*) (tk(j),j=1,ntr)
  Read (9,*) compsi
  If(NSP <= 3) Then
     Read (9,*) XXc(1)
  Else If(NSP == 4) Then
     Read (9,*) XXc(1:NSP-2)
  Endif
  do i=1, nsp
     read(9,*) rct(i,1:nsp)
  enddo
  Read(9,*) gama
  Close(9)
  lx = 1
  do i=1, nsp
     do j=1, i
        rcut(lx) = rct(i,j)
        lx=lx+1
     enddo
  enddo

  If(NSP > 4) Stop " *** N0 de especies > 4 *** "

  Allocate (R(NR))
  Allocate (Q(NR))
  Allocate (SINe(NR))
  Allocate (COSe(NR))
  Allocate (FR(NR))
  Allocate (FI(NR))

  !   calculo de las fracciones molares restantes
  !       Suma(XXC(1:NSP)) = XXc(1) + XXc(2) + .... = 1.00
  !       Suma(XXc(1:NSP)*QQ(1:BSP)) = QQ(1)*XXc(1) + QQ(2)*XXC(2) + ....  = 0.00
  !       Suma(Rho(1:NSP)*QQ(1:NSP)) = Rho(1)*QQ(1) + Rho(2)*QQ(2) + ....  = 0.00

  If(NSP == 2) Then
     XXc(1) = 1.0d0/(1.00d0 - QQ(1)/QQ(2))
     XXc(2) = 1.0d0 - XXc(1)
  Else If(NSP == 3) Then
     XXc(3) = (QQ(2)*(XXc(1) - 1.00d0) - QQ(1)*XXc(1))/(QQ(3) - QQ(2))
     XXc(2) = 1.00d0 - XXc(1) - XXc(3)
  Else If(NSP == 4) Then
     XXc(3) = (QQ(3)*(XXc(1) + XXc(2) - 1.0d0) - XXc(1)*QQ(1) - XXc(2)*QQ(2))/(QQ(4) - QQ(3))
     XXc(4) = 1.00d0 - XXc(1) - XXc(2) - XXc(3)
  End If

  !   inicializacion de funciones

  dcf(:,:)   = 0.00d0
  s(:,:)     = 0.00d0
  sint(:,:)  = 0.00d0
  ipg(1:NIT) = .False.

  !     initialize output


  !     initialize parameters for ft's

  NU = Int(Log(Dble(NR))/Log(2.0d0) + 0.540)
  NR = 2**NU
  N = NR - 1
  rin = 1.0d0/Dble(4*NR)
  a1 = pi/Dble(NR)
  Do i = 1 , NR
     SINe(i) = Sin(Dble(i)*a1)
     COSe(i) = Cos(Dble(i)*a1)
  Enddo
  DQ = pi/(Dble(NR)*DR)
  CR1 = pi*DR
  CQ1 = DQ/(8.0d0*pi**2)

  !     initialize control matrices for NR iterations

  mx = 1
  Do i = 1 , NSP
     ipg((i*i+i)/2) = .True.
     Do j = 1 , i
        npos(i,j) = mx
        npos(j,i) = mx
!        If ( mx /= (i*i+i)/2 ) then
           !
           ! LB mixing rules for LJ parameters
!           sigma(mx) = (sigma((i*i+i)/2)+sigma((j*j+j)/2))/2.00d0
!           eps(mx) = sqrt(eps((i*i+i)/2)*eps((j*j+j)/2))
           !
 !       Endif
        sigma(mx) = ssig(i,j)
        eps(mx)   = eeps(i,j)
        nx = 1
        Do k = 1 , NSP
           Do l = 1 , k
              If ( k > i ) Then
                 in2 = (k*k-k)/2 + i
              Else
                 in2 = (i*i-i)/2 + k
              Endif
              If ( l > j ) Then
                 in1 = (l*l-l)/2 + j
              Else
                 in1 = (j*j-j)/2 + l
              Endif
              ic1(mx,nx) = in1
              ic2(mx,nx) = in2
              nx = nx + 1
           Enddo
        Enddo
        mx = mx + 1
     Enddo
  Enddo

  Write (3,99003) title , (j,qq(j),j,sigma((j*j+j)/2),j,XXc(j),j=1,NSP)
  Write (8,99003) title , (j,qq(j),j,sigma((j*j+j)/2),j,XXc(j),j=1,NSP)
  Write (*,99003) title , (j,qq(j),j,sigma((j*j+j)/2),j,XXc(j),j=1,NSP)
  Write (3,99004) tdiff , pmixt , NR , DR , num , diffnr
  Write (8,99004) tdiff , pmixt , NR , DR , num , diffnr
  Write (*,99004) tdiff , pmixt , NR , DR , num , diffnr
  If(NSP == 2) Then
     Write (8,99023) (j,j=1,NSP)
     Write (8,99024)
  Else If(NSP == 3) Then
     Write (8,99025) (j,j=1,NSP)
     Write (8,99026)
  Else If(NSP == 4) Then
     Write (8,99013) (j,j=1,NSP)
     Write (8,99014)
  End If


  !   initialize cutoff radii

  Do j = 1 , NIT
     ncore(j) = Nint(0.6*sigma(j)/DR)

  Enddo

  Do i = 1 , NR
     R(i) = Dble(i)*DR
     Q(i) = Dble(i)*DQ
  Enddo

  !********************************************************************
  !                                                                   *
  !    STEP 1: Enter or compute initial estimate of s(r)		    *
  !               						    *
  !********************************************************************

  If ( compsi ) Then

     !     msa initial s(r) (reformulate rho(i) for nsp > 2)

     rhot = densid(1)

     rho(1:NSP) = Rhot*XXc(1:NSP)

     If(NSP == 2) Then
        Call MSA(s,NR,rho,qq,1/tk(1),sigma(1))
     Else
        s(1:NR,1:NIT) = 0.00d0
     End If
  Else

     !   read initial estimate from file

     Open (10,FILE='solin.hncr',STATUS='old')
     Read (10,*) rhot , dri , ni , tri , (gammi(j),j=1,NIT)
     Write (3,99016)
     Write (*,99016)
     If ( Abs(DR-dri) > 1.d-6 ) Then
        Write (3,99017) 
        Write (*,99017)
        nrmax = Min(ni,NR-1)
        Do i = 1 , nrmax
           Read (10,*) rr , (sint(i,j),j=1,NIT)
        Enddo
        nrmax = ni
        If ( R(NR-1) < rr ) nrmax = NR - 1
        Do i = 1 , nrmax
           inr = Int(R(i)/dri)
           Do j = 1 , NIT
              s(i,j) = sint(inr,j) + (sint(inr+1,j)-sint(inr,j))&
                   &      *(Dble(i)*DR/dri-inr)
           Enddo
        Enddo
     Else
        Do i = 1 , Min(ni,NR)
           Read (10,*) rr , (s(i,j),j=1,NIT)
        Enddo
     Endif
  Endif

  !*********************************
  !     loop over densities        *
  !*********************************

  Do idn = 1 , nrho
     rhot = densid(idn)

     !    compute ion densities from total density using charge neutrality
     !    condition (only for nsp = 2)

     Rho(1:NSP) = Rhot*XXc(1:NSP)

     Write (*,99018) Rho(1:NSP) 
     Write (*,99022)  pi*(Sum(Rho(1:NSP)*sigma(1:NSP)**3))

     rkx = Sum(Rho(1:NSP)*QQ(1:NSP)**2)

     !*********************************
     !     loop over temperatures     *
     !*********************************

     Do itr = 1 , ntr
        tr = 1.0d0/tk(itr)
        !
        ! Compute LJ potential and virial
        !
        do j = 1, nit
           ep4 = eps(j)*4.d0/tr
           sig = sigma(j)
           ucut(j) = ep4*(sig/rcut(j))**6*((sig/rcut(j))**6-1.0d0)
           do i = 1,nr
              r6 = (sig/r(i))**6
              r12 = r6*r6
              if (r(i) .le.rcut(j)) then
                 ulj(i,j) = ep4*(r12-r6)-ucut(j)
                 wlj(i,j) = -ep4*(12.d0*r12-6.d0*r6)
              Endif
           Enddo
        Enddo
        do i=1,nr
           write(777,'(8f15.7)')r(i),exp(-ulj(i,1:nit))
           write(778,'(8f15.7)')r(i),wlj(i,1:nit)
        Enddo



        !     compute debye inverse screening length

        rkd = Sqrt(4*pi*rkx/tr)

        !     compute adequate yukawa parameter to supress singularity at r=0
        !     in short range functions

        If ( rkd > 1.0d0 ) Then
           xe = 4*rkd
        Else
           xe = 2.00d0
        Endif
        xp2 = (xe**2 - Sqrt(xe**4 - 4*(xe*rkd)**2))/2.00d0
        xm2 = (xe**2 + Sqrt(xe**4 - 4*(xe*rkd)**2))/2.00d0
        xp = Sqrt(xp2)
        xm = Sqrt(xm2)
        Write (3,99012) 1/tr , rhot
        Write (*,99012) 1/tr , rhot
        Write (3,99001)
        Write (*,99001)
        mx = 1

        !     Set density and temperature dependent parameters

        Do id1 = 1 , NSP
           Do id2 = 1 , id1

              !     cr(i),cq(i) are density scaled factors for forward and backward
              !     FFT respectively

              cr(mx)   = CR1*Sqrt(rho(id1)*rho(id2))
              cq(mx)   = CQ1/Sqrt(rho(id1)*rho(id2))
              dnt      = cr(mx)/CR1
              gamm(mx) = -qq(id1)*qq(id2)*gama/tr
              cv       = 4*pi*dnt*gamm(mx)*xe**2
              csv      = gamm(mx)*xe**2/(xm2-xp2)
              rkdx2    = (rkd*xe)**2
              Do i = 1 , NR
                 qi = Q(i)
                 ri = R(i)

                 !******************************************************************
                 !          STEP 2: Compute and store long range corrections,      *
                 !                  v(k) and [v(k)]**(-1) elements of the OZ eq.   *
                 !******************************************************************

                 !   compute the lr behaviour of c(12) and h(12) (note that
                 !   a yukawa term is added so as to make c(12) - clr(12) finite
                 !   at r = 0).

                 clr(i,mx)  = gamm(mx)*(1-Exp(-xe*ri))
                 sclr(i,mx) = csv*(Exp(-xp*ri)-Exp(-xm*ri))

                 !    compute the chain sum hypervertex v(k) matrix and its inverse

                 vk(i,mx)   = cv/(qi**4+(qi*xe)**2+rkdx2)
                 vki(i,mx)  = -cv/(qi**4+(qi*xe)**2)
                 If ( ipg(mx) ) Then
                    vk(i,mx)  = 1 + vk(i,mx)
                    vki(i,mx) = 1 + vki(i,mx)
                 Endif
              Enddo
              mx = mx + 1
           Enddo
        Enddo
        ncmin = NR
        Do j = 1 , NIT
           ncmin = Min(ncmin,ncore(j))

           !       transform series function

           Call FFTSIN(s(1,j),cr(j),ts(1,j))
        Enddo

        ncmin = ncmin - 10

        !     begin iterative solution

        iter = 1

        !     begin broyles cycle

        !     initiate loop over interactions

20      Do j = 1 , NIT

           !******************************************************************
           !         STEP 3: Closure relation (c(r) and h(r) computed here)  *
           !******************************************************************

           !     calculate  the direct correlation function  (hnc closure)
           !     (lr are subtracted)

           Do i = 1 , ncore(j) - 1
              dcf(i,j) = -R(i) - sclr(i,j) - s(i,j)
           Enddo
           Do i = ncore(j) , N
              dcf(i,j) = R(i)*(Exp(-ulj(i,j)+(s(i,j) +gamm(j) - clr(i,j)+sclr(i,j))&
                   &    /R(i)) - 1.0d0) - sclr(i,j) - s(i,j)
           Enddo

           !     compute h(12) (including the lr behaviour)

           h(1:N,j) = (s(1:N,j)+dcf(1:N,j)+sclr(1:N,j))/R(1:N)

           !*******************************************************************
           !     STEP 4: Transform direct correlation function c(12)          *
           !             (the core discontinuity must be treated adequately)  *
           !*******************************************************************

           dcflow = -R(ncore(j)) - sclr(ncore(j),j) - s(ncore(j),j)
           dcfc = dcf(ncore(j),j)
           dcf(ncore(j),j) = (dcf(ncore(j),j)+dcflow)/2.0d0
           Call FFTSIN(dcf(1,j),cr(j),tdcf(1,j))
           dcf(ncore(j),j) = dcfc

           !*******************************************************************
           !     STEP 5: Initialize functions for NR loop                     *
           !*******************************************************************

           !       compute d(l) functions (hnc closure)

           hc = h(ncore(j),j)
           hlow = -1.0d0 - sclr(ncore(j),j)
           h(ncore(j),j) = (h(ncore(j),j)+hlow)/2.0d0
           Call FFTCOS(0.0d0,h(1,j),d(1,j))

           suma = Sum(h(1:N,j))

           d(0,j) = 4.00d0*suma
           h(ncore(j),j) = hc

           !      Initial NR guess

           tsn(1:Num,j) = ts(1:Num,j)
        Enddo

        !       loop over interactions closed above

        !*******************************************************************
        !        STEP 6: Newton Raphson procedure over j=1,..NUM points    *
        !*******************************************************************

        newton = 0
        Do While ( .True. )

           !     start newton raphson iterations

           !     6.1 *** Check for excess of CPU time

           If ( kwit == 1 ) Then
              Write (3,99019) newton , iter 
              Write (*,99019) newton , iter
              Stop
           Endif

           !    6.2 *** Compute the jacobian matrix

           Do i = 1 , num
              Do j = 1 , NIT
                 suma = 0.0d0
                 Do k = 1 , num
                    ik = Abs(i-k)
                    suma = suma + (d(ik,j)-d(i+k,j))*(tsn(k,j)-ts(k,j))
                 Enddo
                 tcc = tdcf(i,j) + suma*rin
                 xg(j) = Q(i)*vki(i,j) - tcc
                 tsr(j) = tcc + Q(i)*vk(i,j)
              Enddo
              Call SSINV(xg,NSP)
              Do j = 1 , NIT
                 tsr(j) = xg(j)*Q(i)**2 - tsr(j)
                 dels(num*(j-1)+i) = tsr(j) - tsn(i,j)
              Enddo

              !       perform the explicit computation and inversion only if
              !       this is the first system

              If ( compj ) Then
                 Do l1 = 1 , NIT
                    Do l2 = l1 , NIT
                       in1 = ic1(l1,l2)
                       in2 = ic2(l1,l2)
                       rj  = -Q(i)**2*xg(in1)*xg(in2)
                       If ( l1 == l2 ) rj = rj + 1.00d0
                       it1 = num*(l1-1)
                       i1  = it1 + i
                       it2 = num*(l2-1)
                       Do k = 1 , num
                          i2 = it2 + k
                          ik = Abs(i-k)
                          rja(i1,i2) = rj*(d(ik,l2)-d(i+k,l2))*rin
                          If ( l1 /= l2 ) rja(i+it2,k+it1)   &
                               &      = rj*(d(ik,l1)-d(i+k,l1))*rin
                          If ( i1 == i2 ) rja(i1,i1) = rja(i1,i1) + 1.00d0
                       Enddo
                    Enddo
                 Enddo
              Endif
           Enddo

           !     invert jacobian (only if this is the first system)

           If ( compj ) Call INVMAT(rja,ND,NIT*num,rjai,work,ipvt)

           !     6.2 *** Compute new s(k) using n-r strategy

           ssum = 0.0d0
           Do i = 1 , num*NIT
              suma = Sum(rjai(i,1:Num*NIT)*dels(1:Num*NIT))
              delg(i) = suma
              ssum = ssum + suma*suma
           Enddo
           delnr = Sqrt(ssum)
           newton = newton + 1
           ntt = ntt + 1

           !     6.3 *** Check convergence and time

           timer = SECOND() - tt
           If ( timer > tmax ) kwit = 1
           If ( delnr > 400 .And. ntt > 1000 ) kwit = 1

           If ( delnr > diffnr ) Then
              Do j = 1 , NIT
                 Do i = 1 , num
                    tsn(i,j) = tsn(i,j) + delg(i+num*(j-1))
                 Enddo
              Enddo

              !     if not converged go back to step 6.1

              Goto 40
           Endif

           !*******************************************************************
           !         STEP 7: Direct OZ equation solution for i=num+1,N        *
           !*******************************************************************

           aa1 = 0.0d0

           !     check for maximum significant k value

           Do kcut = N , num + 1 , -1
              Do j = 1 , NIT
                 aa1 = Max(aa1,Abs(tdcf(kcut,j)))
              Enddo
              If ( aa1 > 5.0d-4 ) Goto 60
           Enddo
           Goto 60
40      Enddo
60      Do i = num + 1 , N

           !       Here the matrix OZ is solved in k-space

           Do j = 1 , NIT
              tcc = tdcf(i,j)
              xg(j) = Q(i)*vki(i,j) - tcc
              tsr(j) = tcc + Q(i)*vk(i,j)
           Enddo
           Call SSINV(xg,NSP)
           tsn(i,1:NIT) = xg(1:NIT)*Q(i)**2 - tsr(1:NIT)
        Enddo

        !*******************************************************************
        ! STEP 8: Transform back to real space s(r) including all k values *
        !*******************************************************************

        Do j = 1 , NIT
           Call FFTSIN(tsn(1,j),cq(j),snew(1,j))
        Enddo

        !*******************************************************************
        !          STEP 9: Check for final convergence			   *
        !*******************************************************************

        diffm = 0.0d0
        Do j = 1 , NIT
           diff(j) = 0.0d0
           Do i = 1 , N
              diff(j) = diff(j) + (snew(i,j) - s(i,j))**2
           Enddo
           diff(j) = Sqrt(DR*diff(j))
           diffm = Max(diffm,diff(j))
        Enddo

        !*******************************************************************
        !     STEP 10: Broyles mixing iterates step  in r- and k-space     *
        !              with self-adjusting mixture parameter		   *
        !*******************************************************************

        Do j = 1 , NIT
           If ( diff(j) > amx ) Then
              alfa = pmixt
           Else
              alfa = 1.0d0 - (1.00d0-pmixt)*(diff(j)/amx)**2
           Endif
           Do i = 1 , N
              s(i,j)  = s(i,j)  + alfa*(snew(i,j) - s(i,j))
              ts(i,j) = ts(i,j) + alfa*(tsn(i,j)  - ts(i,j))
           Enddo
        Enddo
        tim = SECOND() - tt
        Write (3,99002) iter , newton , diffm , tim
        Write (*,99002) iter , newton , diffm , tim

        !    Branch according to convergence

        If ( diffm < tdiff ) Then

           !     end broyles cycle

           !     set compj so to avoid recomputation of the jacobian inside
           !     the NR loop

           compj = .false.

           !*******************************************************************
           !             Compute the jacobian matrix for the last iteration   *
           !             for each system to guarantee stability		   *
           !*******************************************************************

           Do i = 1 , num
              Do j = 1 , NIT
                 suma = 0.0d0
                 Do k = 1 , num
                    ik = Abs(i-k)
                    suma = suma + (d(ik,j)-d(i+k,j))*(tsn(k,j)-ts(k,j))
                 Enddo
                 tcc = tdcf(i,j) + suma*rin
                 xg(j) = Q(i)*vki(i,j) - tcc
                 tsr(j) = tcc + Q(i)*vk(i,j)
              Enddo
              Call SSINV(xg,NSP)
              Do j = 1 , NIT
                 tsr(j) = xg(j)*Q(i)**2 - tsr(j)
                 dels(num*(j-1)+i) = tsr(j) - tsn(i,j)
              Enddo
              Do l1 = 1 , NIT
                 Do l2 = l1 , NIT
                    in1 = ic1(l1,l2)
                    in2 = ic2(l1,l2)
                    rj = -Q(i)**2*xg(in1)*xg(in2)
                    If ( l1 == l2 ) rj = rj + 1.00d0
                    it1 = num*(l1-1)
                    i1 = it1 + i
                    it2 = num*(l2-1)
                    Do k = 1 , num
                       i2 = it2 + k
                       ik = Abs(i-k)
                       rja(i1,i2) = rj*(d(ik,l2)-d(i+k,l2))*rin
                       If ( l1 /= l2 ) rja(i+it2,k+it1) = rj*(d(ik,l1)-d(i+k,l1))*rin
                       If ( i1 == i2 ) rja(i1,i1) = rja(i1,i1) + 1.00d0
                    Enddo
                 Enddo
              Enddo
           Enddo

           !     Invert jacobian and sore it in system matrix (rjai)

           Call INVMAT(rja,ND,NIT*num,rjai,work,ipvt)

           !*******************************************************************
           !                  STEP 11: compute thermodynamics                 *
           !*******************************************************************

           xcsum = 0.00d0
           Do j = 1 , NIT
              ncj = ncore(j)
              dnt2 = (cr(j)/CR1)**2
              ic = 0
              If ( ipg(j) ) ic = 1

              !     restore direct correlation function (subtracting the
              !     coulomb tail for computation of isothermal compressibility)
              !     and fourier transform of h(r)

              Do i = 1 , N
                 thlr     = (vk(i,j) - Dble(ic))*Q(i)
                 th(i,j)  = ts(i,j) + tdcf(i,j) + thlr
                 dcf(i,j) = dcf(i,j) + clr(i,j) - gamm(j)
              Enddo

              !       contributions to energy, a/nkt and chemical potentials

              sncj = s(ncj,j) + sclr(ncj,j) - clr(ncj,j) + gamm(j)
              rms1 = -R(ncj) - sncj

              !       Backward Simpson's rule for r < or = sigma_ij

              Do i = ncj - 1 , 1 , -1
                 ics  = 2*(Mod(ncj-i,2)+1)
                 rms1 = rms1 + ics*dcf(i,j)*R(i)
              Enddo
              esum   = -R(ncj)*gamm(j)*h(ncj,j)
              hij2   = (h(ncj,j)*R(ncj))**2
              a1sum  = -(h(ncj,j)+1)*dcf(ncj,j)*R(ncj)
              rmusum = R(ncj)*h(ncj,j)*sncj - 2*R(ncj)*dcf(ncj,j)
              rms2   = dcf(ncj,j)*R(ncj)
              alrsum = hij2

              !       Forward Simpson's rule for r > sigma_ij

              Do i = ncj + 1 , N
                 ri     = R(i)
                 hij    = h(i,j)*ri
                 sij    = hij - dcf(i,j)
                 ics    = 2*(Mod(i-ncj,2)+1)
                 esum   = esum  -ics*gamm(j)*hij
                 hij2   = hij*hij
                 a1sum  = a1sum - ics*(hij+ri)*dcf(i,j)
                 alrsum = alrsum + ics*hij2
                 rmusum = rmusum + ics*(hij*sij-2*ri*dcf(i,j))
                 rms2   = rms2 + ics*dcf(i,j)*ri
              Enddo
              psum = 0
              eljsum = 0
              do  i=ncore(j),n-2,2
                 l = i+1
                 k = i+2
                 psum = r(i)**2*wlj(i,j)*(h(i,j)+1.d0)+4*r(l)**2&
                      &*wlj(l,j)*(h(l,j)+1.d0)+r(k)**2*wlj(k,j)&
                      &*(h(k,j)+1.d0)+psum 
                 eljsum = r(i)**2*ulj(i,j)*(h(i,j)+1.d0)+4*r(l)**2&
                      &*ulj(l,j)*(h(l,j)+1.d0)+r(k)**2*ulj(k,j)&
                      &*(h(k,j)+1.d0)+eljsum 
              Enddo
              pr(j) = -2.d0*pi*psum*dr/9.d0
              ekt(j)  = 2*pi*esum*DR/3 + pi*gamm(j)*sigma(j)**2
              ektlj(j) = 2*pi*eljsum*DR/3
              a1kt(j) = 2*pi*a1sum*DR/3
              chm(j)  = 2*pi*(rmusum-rms1)*DR/3
              xcsum   = (2.0d0 - ic)*(rms1 + rms2)*dnt2 + xcsum
              aklr(j) = 2*pi*alrsum*DR/3
           Enddo

           a2kt = 0.00d0

           !     compute integral in k-space for a/nkt (Simpson's rule)

           Do i = 1 , N
              qi = Q(i)
              ics = 2*(Mod(i,2)+1)
              trachk = 0.00d0
              fklr = 0.00d0
              fklrii = 0.00d0
              Do j = 1 , NIT
                 thij = th(i,j)/qi
                 ic = 0
                 If ( ipg(j) ) Then
                    ic = 1
                    trachk = trachk + thij
                    fklrii = fklrii + thij*thij
                 Endif
                 xg(j) = ic + thij
                 fklr  = fklr + thij*thij
              Enddo
              fklr = -(fklr-fklrii) - 0.5*fklrii
              Call SDET(xg,NSP,dethk)
              a2kt = a2kt + ics*(Log(dethk)-trachk-fklr)*qi*qi
           Enddo
           a2kt = -a2kt*DQ/(3.0d0*4*rhot*pi**2)
           sek   = 0.00d0
           seklj = 0.00d0
           psum  = 0.00d0
           ankt  = 0.00d0
           akcor = 0.00d0

           !     add up contributions from each interaction

           l = 1
           Do i = 1 , NSP
              Do j = 1 , i
                 ipk     = j/i
                 ic      = 2 - ipk
                 fcnr    = Dble(ic)
                 a1kt(l) = a1kt(l)*rho(i)*rho(j)/rhot
                 aklr(l) = fcnr*aklr(l)*rho(i)*rho(j)/rhot
                 pr(l)   = rho(i)*rho(j)*pr(l)/rhot 
                 sek     = sek + ic*ekt(l)*rho(i)*rho(j)/rhot
                 seklj     = seklj + ic*ektlj(l)*rho(i)*rho(j)/rhot
                 ankt    = ankt + a1kt(l)*ic
                 psum    = psum + pr(l)*ic
                 akcor   = akcor + aklr(l)
                 tds(i,j) = 0
                 do k=1, nr
                    tds(i,j) = tds(i,j)+ dcf(k,l)*r(k)
                 Enddo
                 tds(i,j)  = 4*pi*tds(i,j)*dr
                 if (i .ne. j) tds(j,i) = tds(i,j)
                 l       = l + 1
              Enddo
           Enddo
           xd = xxc(1)
           xion = xxc(2)+xxc(3)
           cionion = 0.5*(tds(2,2)+tds(2,3))
           Mrr0 = 1-rhot*(xion**2*cionion+xd**2*tds(1,1)+2*xion*xd*tds(1,2))
           Mcc0 = 1 -rhot*xion*xd*(cionion+tds(1,1)-2*tds(1,2))
           Mrc0 = rhot*sqrt(xion*xd)*(xd*tds(1,2)-xion*cionion-(xd&
                &-xion)*tds(1,2))
           lamb1 = (Mrr0+Mcc0-sqrt((Mrr0-Mcc0)**2+4*Mrc0**2))/2
           lamb2 = (Mrr0+Mcc0+sqrt((Mrr0-Mcc0)**2+4*Mrc0**2))/2
           vrr = (lamb1-Mcc0)/sqrt((lamb1-Mcc0)**2+Mrc0**2)
           vcc = Mrc0/sqrt((lamb1-Mcc0)**2+Mrc0**2)
           print *, (Mrr0-lamb1)*vrr+Mrc0*vcc
           print *, (Mcc0-lamb1)*vcc+Mrc0*vrr
           sm1  = 0.00d0
           sm11 = 0.00d0
           sm3  = 0.00d0
           chss = 0.00d0

           !*******************************************************************
           !              Evaluation of final thermodynamic properties        *
           !*******************************************************************

           !     add hard core and lr contributions

           Do is = 1 , NSP
              ssm1 = 0.0d0
              chsum = 0.0d0
              Do js = 1 , NSP
                 nps   = npos(is,js)
                 sigx  = sigma(nps)
                 sm1   = sm1 + rho(is)*rho(js)*sigx**3
                 ssm1  = ssm1 + rho(js)*sigx**3
                 chsum = rho(js)*(chm(nps)+ekt(nps)) + chsum
              Enddo
              sm3  = sm3 + (rho(is)*qq(is)*qq(is))**2
              sm11 = sm11 + rho(is)*rho(is)*(sigma(npos(is,is)))**3
              eik  = Exp(-rkd*R(N))

              !     1. Excess chemical potentials (reference ideal gass)

              chemp(is) = 0.5*eik*(0.5*eik-1)*rkd*qq(is)**2  &
                   &       /tr + chsum + 2*pi*ssm1/3
              chss      = chss + chemp(is)*rho(is)
           Enddo
           alr = (2*rkd**3/(8*pi*rhot)+pi*sm3/(rkd*rhot*tr**2))  &
                &               *Exp(-2*rkd*R(N))
           asr = pi*(2*(sm1-sm11)+2*sm11)/(3*rhot)
           elr = -rkd**3*Exp(-rkd*R(N))/(8*pi*rhot)

           !     2. Reduced inverse isothermal compressibility

           xc = 1.00d0 - 4*pi*xcsum*DR/(3*rhot)
           write(*,'(3f15.7)')Mrr0, Mcc0, Mrc0
           write(*,'(3f15.7)')lamb1, lamb2
           det=1.0
           do i=1,nsp
              det=det*tds(i,i)
           Enddo


           !     3. Excess internal energy (Coulomb+LJ)

           sekc = sek + elr 
           sek = sekc+seklj

           !     4. Excess Helmholtz free energy (reference ideal gas)

           ankt = ankt + a2kt + asr + sekc + akcor + alr

           !     5. Compressibility factor from virial equation

           pres = 1.0d0 + psum + sekc/3.0d0

           !     6. Compressibility factor from A/NkT and chemical potential

           za = 1 + chss/rhot - ankt

           tfin = Second() - tt

           !*******************************************************************
           !                 STEP 12: Final print out			   *
           !*******************************************************************

           !     output thermodynamics

           Write (3,99009) ntt , iter , tfin
           Write (8,99015) XXc(1) , rhot , 1/tr , pres , sek , sekc, ankt , &
                &                         za , xc , (chemp(j),j=1&
                &,NSP), lamb1, lamb2, vrr, vcc

           !     print out correlation functions

           !     store short ranged s(r) and normal s(r)

           Open (2,FILE='lsol.hncr',STATUS='unknown')

           Write (2,99031) rhot , DR , N , tr , (gamm(j),j=1,NIT)
           Do i = 1 , N
              Write (2,99021) R(i) , (s(i,j),j=1,NIT)
           Enddo

           Close (2)

           !     print out gab(r)

           Open (4,FILE='RDFhncr.dat',STATUS='unknown')

           tit = 'g'
!!$           Write (4,99032) title , (j,qq(j),j,sigma((j*j+j)/2), &
!!$                &      j,XXc(j),j=1,NSP)
!!$           Write (4,99033) tdiff , pmixt , NR , DR , num , diffnr
!!$           Write (4,99034)
!!$           print *, ' Pero bueno ...'
!!$           If(NSP == 2) Then
!!$              print *, ' bien'
!!$              Write (4,99005) ((tit,j,i,j=1,i),i=1,NSP)
!!$              Write (4,99006)
!!$              print *, ' o??'
!!$           Else If(NSP == 3) Then
!!$              Write (4,99027) ((tit,j,i,j=1,i),i=1,NSP)
!!$              Write (4,99028)
!!$           Else If(NSP == 4) Then
!!$              Write (4,99029) ((tit,j,i,j=1,i),i=1,NSP)
!!$              Write (4,99030)
!!$           End If

           Do i = ncmin , N
              Write (4,99020)  R(i) , (h(i,j)+1,j=1,NIT)
           Enddo
           Close (4)
        Else

           iter = iter + 1
           Goto 20
        Endif
     Enddo
  Enddo
  open(99,file='skhnc.dat')
!!$  tit = 's'
!!$  If(NSP == 2) Then
!!$     Write (99,99005) ((tit,j,i,j=1,i),i=1,NSP),'S_NN','S_cc'
!!$     Write (99,99006)
!!$  Else If(NSP == 3) Then
!!$     Write (99,99027) ((tit,j,i,j=1,i),i=1,NSP),'S_NN','S_cc'
!!$     Write (99,99028)
!!$  Else If(NSP == 4) Then
!!$     Write (99,99029) ((tit,j,i,j=1,i),i=1,NSP), 'S_NN','S_cc'
!!$     Write (99,99030)
!!$  End If

  do i=1, n
     mx = 1
     sksum = 0
     do k=1, nsp
        do l=1, k
           if (k==l) Then
              sk(i,mx) = rho(l)*(1+(th(i,mx))/q(i))/rhot
              sksum = sksum + sk(i,mx)
           else
              sk(i,mx) = sqrt(rho(k)*rho(l))*((th(i,mx))/q(i))/rhot
              sksum = sksum + 2*sk(i,mx)
           endif

           mx=mx+1
        Enddo
     Enddo
     call skcc(sk(i,1:nit),rho(1:nsp),rhot,3,nsp,scccon,nit)
     write(99,99020)q(i),sk(i,1:nit), sksum, scccon
  Enddo



  Stop " *** This is then end, my only friend, the end ... *** "

  !*******************************************************************
  !                 Formats and error messages section               *
  !*******************************************************************

99000 Format(a80)
99001 Format (/2x,'broyles it.   n-r it      error      time'/1x,50('='))
99002 Format (3x,i3,10x,i3,6x,f9.6,2x,f8.2)
99003 Format (/,'#      program hncr: '/"#   ",a80,/"#  ",60('*')&
       &        /6('#   qq(',i1,') = ',f9.4,' sigma(',i1,') = ',f7.4 , &
       &           ' Xc(',i1,') = ',f8.5/:))
99004 Format ('# tdiff = ',f10.7,' pmixt = ',f6.3,' nr = ',i5,' dr = ',f7.4,&
       &        ' num = ',i3,' diffnr = ',f10.7//)
99005 Format (/'#     r',12x,3(a1,'(',2I1,')',7x:,2a6))
99006 Format ('# ',49('-'))
99008 Format (/' *** error: grid size must meet the core conditions '/)
99009 Format (/1x,' * end of computation after ',i4,' n-r iter.',&
       &        ' and ',i3,' broyles iter.'/4x,'cpu time :',f7.2,' sec.'/)
99012 Format (/2x,'computing for beta* = ',f9.4,5x,'rho (tot) = ',f8.5/1x,&
       &        60('-'))
99013 Format ("#  ",'Thermodynamics:'/"# ",17('-')/"#   ",'time',3x,'rho',2x,&
       &        '  beta* ',4x,'z',5x,'u/nkt',3x,'u_c/nkt',3x,'a/nkt',5x,'z(a)',6x,'xc',&
       &        4x,4('mu/kt(',i1,')',2x))
99014 Format ("# ",108('*'))
99015 Format (1x,f6.2,1x,3(f7.4,1x),7(f8.4,1x),1x,6(f8.4,2x))
99016 Format(/' *** Reading initial solution from disk ***.',/)
99017 Format(/'  ** Interpolating S(I,J) to match input dr **'/)
99018 Format(/ " Calculating for : ",// " rho(1) = ",f10.6 , &
       &      /," rho(2) = ",f10.6,/," rho(3) = ",f10.6, & 
       &      /," rho(4) = ",f10.6)
99019 Format(/"  *** convergence not achieved after ", i3, &
       & " n-r iter and ",i3," broyles iter "/)
99020 Format(1x,f9.5,3x,10(1x,f11.5))
99021 Format(f10.5,4x,10(f12.5,1x))
99022 Format(" eta    = ",f10.6/)
99023 Format ("#  ",'Thermodynamics:'/2x,17('-')/"#   ",'x(1)',3x,'rho',2x,&
       &        '  beta* ',5x,'z',5x,'u/nkt',3x,'u_c/nkt',3x,'a/nkt',5x,'z(a)',4x,'xc',&
       &        5x,2('mu/kt(',i1,')',2x))
99024 Format ("#",88('*'))
99025 Format ("#  ",'Thermodynamics:'/"# ",17('-')/"#   ",'x(1)',3x,'rho',2x,&
       &        '  beta* ',5x,'z',5x,'u/nkt',3x,'u_c/nkt',3x,'a/nkt',5x,'z(a)',4x,'xc',&
       &        5x,3('mu/kt(',i1,')',2x))
99026 Format ("#",98('*'))
99027 Format (/'#     r',12x,6(a1,'(',2I1,')',7x:),2A7)
99028 Format ('# ',134('-'))
99029 Format (/'#     r',12x,10(a1,'(',2I1,')',7x:),2A7)
99030 Format ('# ',134('-'))
99031 Format(2(f8.5,1x),i5,1x,11(f14.9))
99032 Format (/"#   ",a80,/"# ",60('*') ,//, &
       &  6("#     qq(",i1,') = ',f9.4,'   sigma(',i1,') = ', &
       &     f7.4,'   Xc(',i1,') = ',f8.5/:))
99033 Format ('#   tdiff = ',f10.7,' pmixt = ',f6.3,' nr = ',i5, &
       &     ' dr = ',f7.4,' num = ',i3,' diffnr = ',f10.7//)
99034 Format("#       FUNCIONES DE DISTRIBUCION      ",/,  &
       & "# ====================================="//)

End Program HNCR





Subroutine MSA(Sout,Nr,Rhoe,Qq,Tr,Sigma)
  Use RQst
  Use Params

  !     programmed by e. lomba, trondheim 23.02.90
  !     Computes the MSA solution for equal sized charged hard spheres
  !     Return value sout(i,j) contains the function h-c.

  Implicit None

  Integer (kind=4) :: i , i1 , i2 , mx , n  , Nr 
  Real (kind=8) :: b , ck , CQ , CR , cst , ct , e , &
       &           rho , rk , rlk , Sigma , sk
  Real (kind=8) :: Tr , v , vk , x , xe
  Real (kind=8) :: k
  Real (kind=8) , Dimension (NR) :: ts , s , ss 
  Real (kind=8) , Dimension (NR,NIT) :: Sout
  Real (kind=8) , Dimension (2) :: Qq , Rhoe

  cr = cr1
  cq = cq1

  rk = -1.00/Tr
  rho = 0.00d0

  rho = Sum(Rhoe(1:2)*Qq(1:2)**2)

  e = rho*pi/6.00d0
  x = Sqrt(24*e/Tr)
  If ( x > 1 ) Then
     xe = 4*x
  Else
     xe = 2.00d0
  Endif
  b = (1.0d0 + x - Sqrt(1.0d0 + 2*x))/x
  v = 2*b/Tr
  n = Nr - 1

  !     compute c(k) and solve oz eq. for ts(k)

  Do i = 1 , Nr
     k  = Q(i)
     sk = Sin(k)
     ck = Cos(k)
     ct = (v*(-k*ck+sk) + v**2*(-k**2*ck+2*k*sk+2*(ck-1))/(4*rk*k))&
          &        /k**2 - rk*ck/k
     ct = -4*pi*ct

     !       subtract fourier transform of the long range component of c(k)

     rlk = 4*pi*rk*(1/k-k/(k**2+xe**2))
     cst = ct - rlk

     !       solve renormalized oz eq.

     vk    = k*rho/(k-rho*rlk)
     ts(i) = k*cst*(vk/rho)**2/(k-cst*vk) - cst
  Enddo

  !     transform to real space

  Call FFTSIN(ts,CQ,s)
  Call PYHS(ss,Nr,rho,1.00d0)

  mx = 1
  Do i1 = 1 , 2
     Do i2 = 1 , i1
        Do i = 1 , n
           Sout(i,mx) = ss(i) + Qq(i1)*Qq(i2)*s(i)
        Enddo
        mx = mx + 1
     Enddo
  Enddo

End Subroutine MSA





Subroutine PYHS(S,Nr,Rho,Dhs)
  Use RQst
  Use Params, Only : Pi

  !     this soubroutine computes the hard spheres function r*s(r),
  !     dhs is the hard sphere diameter. rho is the hard sphere density
  !     the approximation used is the percus-yevick.
  !     Programmed by E.Lomba Trondheim 20.09.90

  Implicit None

  Integer (kind=4) :: k , n , Nr
  Real (kind=8) :: ck , cpi , CQ , CR , Dhs , et , f1 , &
       &                 f2 , f3 , qi , qi2 , qk , qr , qr2 , &
       &                 qr3 , qr4
  Real (kind=8) :: Rho , rl1 , rl2 , rtc , S(*) , sk , tdcf , ts(NR)

  cr  = cr1
  cq  = cq1
  et  = pi*Rho*Dhs**3/6.00d0
  cpi = -4.00d0*pi*Dhs*Dhs
  rl1 =  (1.00d0 + 2.00d0*et)**2/(1.00d0-et)**4
  rl2 = -(1.00d0 + 0.5D0*et)**2/(1.00d0-et)**4
  n = Nr - 1
  Do k = 1 , n
     qk    = Q(k)
     qi    = qk*Dhs
     qi2   = qi*qi
     sk    = Sin(qi)
     ck    = Cos(qi)
     qr    = 1.00d0/qi
     qr2   = qr*qr
     qr3   = qr2*qr
     qr4   = qr2*qr2
     f1    = (sk-qi*ck)*qr2
     f2    = (-2.00d0 + 2.00d0*qi*sk-qi2*ck+2.00d0*ck)*qr3
     f3    = (4.0d0*(6.0d0+qi*(qi2 - 6.0d0)*sk) + &
          &      (12.0d0*qi2 - qi2*qi2 - 24.0d0)*ck)*qr4*qr
     tdcf  = cpi*(rl1*f1 + 6.0d0*et*rl2*f2 + 0.50d0*et*rl1*f3)
     rtc   = Rho*tdcf
     ts(k) = tdcf*rtc/(qk-rtc)
  Enddo

  Call FFTSIN(ts,CQ,S)

End Subroutine PYHS







Subroutine FFTSIN(F,C,Tf)
  Use FFT

  Implicit None

  Integer (kind=4) :: i , IBITR , ii , k , k0 , k1 , k2 , l , nr1 
  Real (kind=8) :: C ,  F(NR)  , Tf(NR) , ti , tr

  !     this subroutine returns
  !                      nr-1
  !       tf(i) = c(i)*4*sumk(f(k)*sin(pi*i*k/nr))

  !     f(0) and f(nr) are both set to 0.

  !     see e.o.brigham, "the fast fourier transform", chapter 10.
  !     Programmed by F.Lado

  k0 = NR/2
  nr1 = NR + 1
  FR(1) = 0.0d0
  F(NR) = 0.0d0
  Do i = 1 , k0
     FR(i+1) = F(i+i)
     FI(i) = F(i+i-1)
     FR(nr1-i) = -FR(i+1)
     FI(nr1-i) = -FI(i)
  Enddo

  k = 0
  Do l = 1 , NU
     Do i = 1 , k0
        k1 = k + 1
        k2 = k1 + k0
        tr = FR(k2)
        ti = FI(k2)
        FR(k2) = FR(k1) - tr
        FI(k2) = FI(k1) - ti
        FR(k1) = FR(k1) + tr
        FI(k1) = FI(k1) + ti
        k = k + 1
     Enddo
     Do While ( .True. )
        k = k + k0
        If ( k < NR ) Then
           Do i = 1 , k0
              ii = IBITR(k/k0,NU)
              ii = ii + ii
              k1 = k + 1
              k2 = k1 + k0
              tr = FR(k2)*COSe(ii) + FI(k2)*SINe(ii)
              ti = FI(k2)*COSe(ii) - FR(k2)*SINe(ii)
              FR(k2) = FR(k1) - tr
              FI(k2) = FI(k1) - ti
              FR(k1) = FR(k1) + tr
              FI(k1) = FI(k1) + ti
              k = k + 1
           Enddo
        Else
           k = 0
           k0 = k0/2
           Goto 100
        Endif
     Enddo
100 Enddo

  Do k = 1 , NR
     i = IBITR(k-1,NU) + 1
     If ( i > k ) Then
        tr = FR(k)
        ti = FI(k)
        FR(k) = FR(i)
        FI(k) = FI(i)
        FR(i) = tr
        FI(i) = ti
     Endif
  Enddo

  Do i = 1 , N
     Tf(i) = C*(-(FI(i+1)-FI(nr1-i))+(FR(i+1)-FR(nr1-i))*COSe(i)&
          &           +(FI(i+1)+FI(nr1-i))*SINe(i))
  Enddo
  Return
End Subroutine FFTSIN







Subroutine FFTCOS(F0,F,Tf)
  Use FFT

  Implicit None

  Integer (kind=4) :: i , IBITR , ii , k , k0 , k1 , k2 , l , nr1
  Real (kind=8) ::  F(NR) , F0  , Tf(NR) , ti , tr

  !     this subroutine returns
  !                         nr-1
  !       tf(i) = 4*(f(0)/2+sumk(f(k)*cos(pi*i*k/nr)))

  !     f0 is f(0) and f(nr) is set to 0.

  !     see e.o.brigham, "the fast fourier transform", chapter 10.
  !     Programmed by F.Lado

  k0 = NR/2
  nr1 = NR + 1
  FR(1) = F0
  F(NR) = 0.0d0
  Do i = 1 , k0
     FR(i+1) = F(i+i)
     FI(i) = F(i+i-1)
     FR(nr1-i) = FR(i+1)
     FI(nr1-i) = FI(i)
  Enddo

  k = 0
  Do l = 1 , NU
     Do i = 1 , k0
        k1 = k + 1
        k2 = k1 + k0
        tr = FR(k2)
        ti = FI(k2)
        FR(k2) = FR(k1) - tr
        FI(k2) = FI(k1) - ti
        FR(k1) = FR(k1) + tr
        FI(k1) = FI(k1) + ti
        k = k + 1
     Enddo
     Do While ( .True. )
        k = k + k0
        If ( k < NR ) Then
           Do i = 1 , k0
              ii = IBITR(k/k0,NU)
              ii = ii + ii
              k1 = k + 1
              k2 = k1 + k0
              tr = FR(k2)*COSe(ii) + FI(k2)*SINe(ii)
              ti = FI(k2)*COSe(ii) - FR(k2)*SINe(ii)
              FR(k2) = FR(k1) - tr
              FI(k2) = FI(k1) - ti
              FR(k1) = FR(k1) + tr
              FI(k1) = FI(k1) + ti
              k = k + 1
           Enddo
        Else
           k = 0
           k0 = k0/2
           Goto 100
        Endif
     Enddo
100 Enddo

  Do k = 1 , NR
     i = IBITR(k-1,NU) + 1
     If ( i > k ) Then
        tr = FR(k)
        ti = FI(k)
        FR(k) = FR(i)
        FI(k) = FI(i)
        FR(i) = tr
        FI(i) = ti
     Endif
  Enddo

  Do i = 1 , N
     Tf(i) = FR(i+1) + FR(nr1-i) + (FI(i+1)+FI(nr1-i))*COSe(i)&
          &           - (FR(i+1)-FR(nr1-i))*SINe(i)
  Enddo
  Return
End Subroutine FFTCOS






Function IBITR(J,Nu)

  Implicit None

  Integer (kind=4) :: i , IBITR , J , j1 , j2 , Nu

  !     Bit reversal function. Programmed by F.Lado

  j1 = J
  IBITR = 0
  Do i = 1 , Nu
     j2 = j1/2
     IBITR = IBITR + IBITR + (j1 - j2 - j2)
     j1 = j2
  Enddo
  Return
End Function IBITR





Subroutine SDET(A,N,Det)

  Implicit None

  Integer (kind=4) :: i , j , j1 , j2 , jm , jp , k , k1 , N
  Real (kind=8) :: A(N) , a1 , Det , sum

  !     this subroutine returns the determinant of the symmetric nxn
  !     matrix a. symmetric storage mode is used.
  !     programmed by f.lado

  If ( A(1) >= 0.0d0 ) Then
     Det = A(1)
     If ( N == 1 ) Return
     A(1) = Sqrt(A(1))
     k1 = 1
     Do k = 2 , N
        A(k1+1) = A(k1+1)/A(1)
        k1 = k1 + k
     Enddo
     j1 = 1
     Do j = 2 , N
        sum = 0.0d0
        jm = j - 1
        Do i = 1 , jm
           sum = sum + A(j1+i)**2
        Enddo
        j2 = j1 + j
        a1 = A(j2) - sum
        If ( a1 < 0.0d0 ) Goto 100
        Det = Det*a1
        A(j2) = Sqrt(a1)
        If ( j /= N ) Then
           jp = j + 1
           k1 = j2
           Do k = jp , N
              sum = 0.00d0
              Do i = 1 , jm
                 sum = sum + A(j1+i)*A(k1+i)
              Enddo
              A(k1+j) = (A(k1+j) - sum)/A(j2)
              k1 = k1 + k
           Enddo
        Endif
        j1 = j2
     Enddo
     Return
  Endif

100 Print * , ' *** error in routine sdet ***'

End Subroutine SDET







Subroutine SPROD(A,B,C,N)

  Implicit None

  Integer i , i1 , is , j , j1 , js , k , ks , N
  Real (kind=8) :: A(10)  , B(10) , C(10) , sum

  !     this subroutine returns the matrix product c = ab for symmetric
  !     nxn matrices a, b, and c. symmetric storage mode is used.

  !     adapted from the ibm "scientific subroutine package" by F.Lado.

  js = 0
  Do j = 1 , N
     is = 0
     Do i = 1 , j
        sum = 0.0d0
        Do k = 1 , i
           sum = sum + A(is+k)*B(js+k)
        Enddo
        ks = is + i
        If ( i /= j ) Then
           i1 = i + 1
           Do k = i1 , j
              sum = sum + A(ks+i)*B(js+k)
              ks = ks + k
           Enddo
        Endif
        If ( j /= N ) Then
           j1 = j + 1
           Do k = j1 , N
              sum = sum + A(ks+i)*B(ks+j)
              ks = ks + k
           Enddo
        Endif
        C(js+i) = sum
        is = is + i
     Enddo
     js = js + j
  Enddo
  Return
End Subroutine SPROD







Subroutine SSINV(Ai,N)

  Implicit None

  Integer (kind=4) :: i , j , k , km1 , kp1 , l , N , r(100)
  Real (kind=8) :: a(6,6) , Ai(10) , big , epsiln , p(100) , pivot1 , q(100) , test

  Data epsiln/1.00d-20/

  !     Inversion of a symmetric matrix. Symmetric storage mode is used.
  !     Orignal source code from CERN computer library.

  If ( N < 1 ) Then
  Elseif ( N == 1 ) Then
     Ai(1) = 1.00d0/Ai(1)
     Return
  Else
     l = 1
     Do i = 1 , N
        Do j = 1 , i
           a(j,i) = Ai(l)
           l = l + 1
        Enddo
     Enddo

     !                   construct truth table

     Do i = 1 , N
        r(i) = 1
     Enddo

     !                   begin programme

     Do i = 1 , N
        k = 0

        !                   search for pivot

        big = 0.00d0
        Do j = 1 , N
           test = Abs(a(j,j))
           If ( test > big ) Then
              If ( r(j) < 0 ) Goto 100
              If ( r(j) /= 0 ) Then
                 big = test
                 k = j
              Endif
           Endif
        Enddo

        !                   test for zero matrix

        If ( k <= 0 ) Goto 100

        !                   test for linearity

        If ( i == 1 ) pivot1 = a(k,k)
        If ( Abs(a(k,k)/pivot1) < epsiln ) Goto 100

        !                   preparation for elimination step1

        r(k) = 0
        q(k) = 1.0d0/a(k,k)
        p(k) = 1.0d0
        a(k,k) = 0.00d0
        kp1 = k + 1
        km1 = k - 1
        If ( km1 < 0 ) Goto 100
        If ( km1 /= 0 ) Then
           Do j = 1 , km1
              p(j) = a(j,k)
              q(j) = a(j,k)*q(k)
              If ( r(j) < 0 ) Goto 100
              If ( r(j) /= 0 ) q(j) = -q(j)
              a(j,k) = 0.0d0
           Enddo
        Endif
        If ( k < N ) Then
           Do j = kp1 , N
              p(j) = a(k,j)
              q(j) = -a(k,j)*q(k)
              If ( r(j) < 0 ) Goto 100
              If ( r(j) == 0 ) p(j) = -p(j)
              a(k,j) = 0.00d0
           Enddo
        Elseif ( k /= N ) Then
           Goto 100
        Endif

        !                   elimination proper

        Do j = 1 , N
           Do k = j , N
              a(j,k) = a(j,k) + p(j)*q(k)
           Enddo
        Enddo
     Enddo

     !                   elements of left diagonal

     l = 1
     Do j = 1 , N
        Do k = 1 , j
           Ai(l) = a(k,j)
           l = l + 1
        Enddo
     Enddo
     Return
  Endif

  !                   failure return

100 Print * , '   *****  failure return from subroutine ssinv  *****'
End Subroutine ssinv








Subroutine INVMAT(A,Lda,N,B,Work,Ipvt)

  !     Computes the  inverse of a matrix
  !     using routines from LAPACK and BLAS
  !     Programmed by E.Lomba as inteface to dgetrf and dgtri routines
  !           ***description

  !     on entry

  !        a       real*8(lda, n)
  !                the output from sgeco or sgefa.

  !        lda     integer
  !                the leading dimension of the array  a .

  !        n       integer
  !                the order of the matrix  a .

  !        ipvt    integer(n)
  !                the pivot vector from sgeco or sgefa.

  !        work    real(n)
  !                work vector.  contents destroyed.

  !     on return

  !        b       inverse of original a matrix

  Implicit None

  Integer i , info , Ipvt(N) , j , Lda , N
  Real (kind=8) :: A(Lda,1) , B(Lda,1)  , det(2) , Work(N)

  Do i = 1 , N
     Do j = 1 , N
        B(i,j) = A(i,j)
     Enddo
  Enddo
  Call DGETRF(N,N,B,Lda,Ipvt,info)
  If ( info == 0 ) Then
     Call DGETRI(N,B,Lda,Ipvt,Work,N,info)
     Return
  Endif

  Stop ' *** error from invmat: a is singular ***'
End Subroutine INVMAT






Subroutine Distribucion
  Use Funciones
  Use FFT, Only : NR
  Use Params, Only : Nit , ND , NSP

  Implicit None
  Allocate (tds(nsp,nsp))
  Allocate (dcf(NR,NIT))
  Allocate (ulj(NR,NIT))
  Allocate (wlj(NR,NIT))
  Allocate (tdcf(NR,NIT))
  Allocate (sk(NR,NIT))
  Allocate (ts(NR,NIT))
  Allocate (h(NR,NIT))
  Allocate (tsn(NR,NIT))
  Allocate (s(NR,NIT))
  Allocate (snew(NR,NIT))
  Allocate (clr(NR,NIT))
  Allocate (sclr(NR,NIT))
  Allocate (vk(NR,NIT))
  Allocate (vki(NR,NIT))
  Allocate (th(NR,NIT))
  Allocate (sint(0:NR,NIT))
  Allocate (d(0:NR,NIT))
  Allocate (rja(ND,ND))
  Allocate (rjai(ND,ND))
  Allocate (dels(ND))
  Allocate (delg(ND))
  Allocate (work(ND))
  Allocate (gamm(NIT))
  Allocate (cr(NIT))
  Allocate (cq(NIT))
  Allocate (tcn(NIT))
  Allocate (xg(NIT))
  Allocate (tc2(NIT))
  Allocate (tsr(NIT))
  Allocate (ekt(NIT))
  Allocate (ektlj(NIT))
  Allocate (sigma(NIT))
  Allocate (eps(NIT))
  Allocate (qq(NSP))
  Allocate (XXc(NSP))
  Allocate (rhoi(NSP))
  Allocate (rho(NSP))
  Allocate (pr(NIT))
  Allocate (tdcf0(NIT))
  Allocate (diff(NIT))
  Allocate (gammi(NIT))
  Allocate (tmp(NIT))
  Allocate (chm(NIT))
  Allocate (a1kt(NIT))
  Allocate (aklr(NIT))
  Allocate (chemp(NSP))
  Allocate (ic1(NIT,NIT))
  Allocate (ic2(NIT,NIT))
  Allocate (ipvt(ND))
  Allocate (ncore(NIT))
  Allocate (imin(NIT))
  Allocate (ipg(NIT))
  Allocate (npos(NSP,NSP))

End Subroutine Distribucion

Subroutine skcc(sk,rho,rhot,itcc,nsp,scc,nit)
  Implicit None
  integer :: itcc, nsp, nit
  real(kind=8) :: sk(nit),rho(nsp),rhot,x(nsp), s11, s12, s22,x1,x2
  integer :: i
  real(kind=8) :: scc
  do i=1, nsp
     x(i)=rho(i)/rhot
  Enddo
  if (itcc==2) Then
     x1=x(1)
     x2=x(2)
     s11=sk(1)
     s22=sk(3)
     s12=sk(2)
  else
     x1=x(1)
     x2=x(2)+x(3)
     s11=sk(1)
     s12=sk(2)+sk(4)
     s22=sk(3)+2*sk(5)+sk(6)
  endif
  scc = x2**2*s11+x1**2*s22- 2*x1*x2*s12
end Subroutine skcc

