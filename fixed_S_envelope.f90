module fixedSprofile
 implicit none

 logical :: do_ledoux = .false.

contains


!-----------------------------------------------------------------------
!+
!  Calculate constant-entropy density and pressure
!
! Indexing convention:
! * N = size(r_prof) - core_idx, such that r(N) is the surface cell face and r(0) = r_prof(core_idx) is the core cell face
! * r(i), m(i) are cell face values at position i+1/2, ranging from i = 0 (core) to N (surface)
! * pres(i), rho(i) are cell-centred values, ranging from i = 0 to N
! * dm(i) := m(i) - m(i-1) is the mass of cell i, ranging from i = 1 to N
! * dm_c(i) := 0.5*( m(i+1) - m(i-1) ) is the mass between cells i and i+1, ranging from i = 0 to N-1
! * Core boundary, rcore, is at r(0) = r_prof(core_idx) while the surface is at r(N) = r(size(r_prof)-core_idx)
!+
!-----------------------------------------------------------------------
subroutine calc_rho_and_pres(m_prof,r_prof,rho_prof,p_prof,X_prof,Y_prof,mu_prof)
 use physcon, only:solarm,solarr,kb_on_mh
 use eos, only:init_eos,entropy,gmw,ieos,gamma,X_in,Z_in
 real, allocatable, dimension(:), intent(inout) :: m_prof,r_prof,rho_prof,p_prof,X_prof,Y_prof,mu_prof
 real, allocatable, dimension(:) :: m,r,rho,pres,dm,dm_c
 real    :: Rstar,rcore,Sc,Sfac,tol,delta,gamma1_0,Rout,r_midenv,s_midenv
 real    :: surfpres,S_lo,S_hi,dS
 integer :: i,core_idx,N,ientropy,counter1,ierr,idx_midenv,maxiter
 logical :: success,have_lo,have_hi


 !-----------------------------------------------------------------------------------------
 ! USER SETTINGS
 !-----------------------------------------------------------------------------------------
 !
 ! Choose core radius
 rcore = 18.5 * solarr
 !
 ! Target surface pressure to shoot for (dyn cm^-2)
 surfpres = 600.
 !
 ! Central Gamma_1 (used when do_ledoux = .true.)
 gamma1_0 = 1.55916
 !
 ! Adjustment factor for shooting
 Sfac = 0.05      ! Step in entropy (units of kB per particle) used while bracketing the target
 tol = 1d-3       ! Relative tolerance for matching the surface pressure
 maxiter = 200    ! Maximum number of shots
 !
 ! EoS options
 ieos = 12
 ientropy = 2  ! Include both gas and radiation entropy
 gamma = 5./3. ! Polytropic index
 !-----------------------------------------------------------------------------------------

 ! Properties in middle of envelope
 r_midenv = 0.5*(r_prof(size(r_prof))+rcore)
 idx_midenv = minloc(abs(r_prof-r_midenv),dim=1)
 X_in = X_prof(idx_midenv)
 Z_in = 1.-X_in-Y_prof(idx_midenv)
 gmw = mu_prof(idx_midenv)
 s_midenv = entropy(rho_prof(idx_midenv),p_prof(idx_midenv),gmw,ientropy)

 print*,'CONSTRUCTING FIXED ENTROPY STAR'
 write(*,'(a30,e13.5,a10)')'rcore =',rcore / solarr,'Rsun'
 write(*,'(a30,e13.5,a10)')'Target surface pressure =',surfpres,'dyn cm^-2'
 write(*,'(a30,e13.5,a10,e13.5,a10,e13.5)')'Composition: X = ',X_in,', Z = ',Z_in,', mu = ',gmw
 write(*,'(a30,e13.5)')'Initial guess for entropy: s = ',s_midenv
 read*
 call init_eos(ieos,ierr)

 core_idx = minloc(abs(r_prof-rcore),dim=1)
 N = size(m_prof) - core_idx

 ! Allocate and set up arrays
 allocate( m(0:N), rho(0:N), pres(0:N), r(0:N) )
 allocate( dm(1:N), dm_c(0:N-1))
 m(0:N) = m_prof(core_idx:)
 rho(0) = rho_prof(core_idx)
 pres(0) = p_prof(core_idx)
 r(0) = r_prof(core_idx)
 do i = 1,N
    dm(i) = m(i) - m(i-1)
 enddo
 do i = 0,N-1
    dm_c(i) = 0.5*(m_prof(core_idx+i+1) - m_prof(core_idx+i-1))
 enddo

 !-----------------------------------------------------------------------------------------
 ! INITIAL GUESSES
 !-----------------------------------------------------------------------------------------
 ! Constant value of entropy
 Sc = s_midenv
 !
 ! Guess for delta factor
 delta = 0.005
 !-----------------------------------------------------------------------------------------
 
 !---------------------------SHOOT-ON-ENTROPY---------------------------------------------
 ! Vary the constant entropy Sc until the surface pressure pres(N) matches surfpres.
 ! pres(N) increases with Sc.
 ! Phase 1: step Sc by dS until the target is bracketed by S_lo (too low) and S_hi (too high).
 ! Phase 2: bisect between S_lo and S_hi.
 dS = kb_on_mh * Sfac / gmw
 have_lo = .false.
 have_hi = .false.
 write(*,'(a4,2x,a8,2x,3(a14,2x))') 'c','status','S','Psurf','R [Rsun]'
 do counter1 = 1,maxiter
    call one_shot(Sc,delta,rcore,m,dm,dm_c,r,rho,pres,ientropy,Rout,success)
   !  if (success) then
   !  else
    endif

    if (success .and. abs(pres(N)-surfpres) < tol*surfpres) exit

    if (.not.success .or. pres(N) < surfpres) then
       S_lo = Sc       ! surface pressure too low: entropy too low
       have_lo = .true.
    else
       S_hi = Sc       ! surface pressure too high: entropy too high
       have_hi = .true.
    endif

    if (have_lo .and. have_hi) then
       if (S_hi - S_lo < 1.e-14*abs(Sc)) then
          print*,'ERROR: entropy bracket collapsed without matching surface pressure'
          stop
       endif
       Sc = 0.5*(S_lo + S_hi)
    elseif (have_lo) then
       Sc = Sc + dS
    else
       Sc = Sc - dS
    endif
 enddo
 if (counter1 > maxiter) then
    print*,'ERROR: surface pressure not matched after maxiter = ',maxiter,' shots'
    stop
 endif
 !-----------------------------------------------------------------------------------------

!  print*,'CONVERGED ONTO DESIRED PROFILE'
 write(*,'(a30,e13.5,a10)') 'Obtained surface pressure =',pres(N),'dyn cm^-2'
 write(*,'(a30,e13.5,a10)') 'Obtained radius =',r(N) / solarr,'Rsun'
 r_prof(core_idx:) = r(0:N)
 rho_prof(core_idx:) = rho(0:N)
 p_prof(core_idx:) = pres(0:N)

end subroutine calc_rho_and_pres


!-----------------------------------------------------------------------
!+
!  One shot
!+
!-----------------------------------------------------------------------
subroutine one_shot(Sc,delta,rcore,m,dm,dm_c,r,rho,pres,ientropy,Rstar,success)
 use physcon, only:gg,pi,solarm,solarr
 use eos, only:get_rho_from_p_s,gmw
 real, intent(in)                               :: Sc,delta,rcore
 real, allocatable, dimension(:), intent(in)    :: m,dm,dm_c
 real, allocatable, dimension(:), intent(inout) :: r,rho,pres
 real, intent(out)                              :: Rstar   ! outermost radius reached
 logical, intent(out)                           :: success ! .false. if the integration broke down before the surface
 integer, intent(in)                            :: ientropy
 integer                                        :: i,N
 real                                           :: gamma1_i,Gamma1_dlnrho_dlnpi

 N = size(m) - 1
 success = .false.
 Rstar = r(0)
 do i = 0,N-1
    pres(i+1) = pres(i) - dm_c(i) * gg*m(i)/(4.*pi*r(i)**4)

    if (pres(i+1)<0.) return
    if (do_ledoux) then
       call get_gamma1_from_rho_p(rho(i),pres(i),gamma1_i)
       Gamma1_dlnrho_dlnpi = Gamma1_dlnrho_dlnp(delta,r(i))
       rho(i+1) = rho(i) + (pres(i+1)-pres(i))*(Gamma1_dlnrho_dlnpi/gamma1_i*rho(i)/pres(i))
    else
       call get_rho_from_p_s(pres(i+1),Sc,rho(i+1),gmw,rho(i),ientropy)
    endif
    if (rho(i+1)<0.) return
    r(i+1) = (3.*dm(i+1)/(4.*pi*rho(i+1)) + r(i)**3)**(1./3.)
    Rstar = r(i+1)
 enddo
 success = .true.

end subroutine one_shot


!-----------------------------------------------------------------------
!+
!  Get Gamma_1 from density and pressure when do_ledoux = .true.
!+
!-----------------------------------------------------------------------
subroutine get_gamma1_from_rho_p(rho,pres,gamma1)
 use eos_mesa, only:get_eos_eT_from_rhop_mesa
 use mesa_microphysics, only:getvalue_mesa
 real, intent(in)  :: rho,pres
 real, intent(out) :: gamma1
 real              :: ui,tempi
 integer           :: ierr

 call get_eos_eT_from_rhop_mesa(rho,pres,ui,tempi)
 call getvalue_mesa(rho,ui,11,gamma1,ierr)

 end subroutine get_gamma1_from_rho_p


!-----------------------------------------------------------------------
!+
!  Calculate pressure given density and Gamma_1, used when setting
!  central boundary conditions
!+
!-----------------------------------------------------------------------
subroutine get_pres_from_rho_gamma1(rho,gamma1,pres,presguess)
 real, intent(in)  :: rho,gamma1,presguess
 real, intent(out) :: pres
 real                :: spres,spres_plus_dspres,gamma1_plus_dgamma1,dGamma1_dspres,gamma1i
 real(kind=8)        :: corr
 real, parameter     :: eoserr=1d-9,dfac=1d-12

 ! We apply the Newton-Raphson method directly to pres^1/2 ("spres") instead
 ! of pres since S(pres) cannot take a negative argument.
 spres = sqrt(presguess) ! Initial guess
 corr = huge(corr);
 do while (abs(corr) > eoserr*abs(spres))
    ! First calculate dGamma_1/dspres
    spres_plus_dspres = spres * (1. + dfac)
    call get_gamma1_from_rho_p(rho,spres_plus_dspres**2,gamma1_plus_dgamma1)
    call get_gamma1_from_rho_p(rho,spres**2,gamma1i)
    dGamma1_dspres = (Gamma1_plus_dGamma1 - gamma1i) / (spres_plus_dspres - spres)
    corr = (gamma1i - gamma1) / dGamma1_dspres
    spres = spres - corr
 enddo
 pres = spres**2

end subroutine get_pres_from_rho_gamma1


!-----------------------------------------------------------------------
!+
!  Calculate pressure given density and entropy using Newton-Raphson
!  method
!+
!-----------------------------------------------------------------------
subroutine get_pres_from_rho_s(rho,S,pres,mu,presguess,ientropy)
 use physcon, only:kb_on_mh
 use eos, only:entropy
 real, intent(in)    :: rho,S,mu,presguess
 real, intent(inout) :: pres
 real                :: spres,spres_plus_dspres,S_plus_dS,dSdspres
 real(kind=8)        :: corr
 real, parameter     :: eoserr=1d-9,dfac=1d-12
 integer, intent(in) :: ientropy
 
 ! We apply the Newton-Raphson method directly to pres^1/2 ("spres") instead
 ! of pres since S(pres) cannot take a negative argument.
 spres = sqrt(presguess) ! Initial guess
 corr = huge(corr);
 do while (abs(corr) > eoserr*abs(spres))
    ! First calculate dS/dspres
    spres_plus_dspres = spres * (1. + dfac)
    S_plus_dS = entropy(rho,spres_plus_dspres**2,mu,ientropy)
    dSdspres = (S_plus_dS - entropy(rho,spres**2,mu,ientropy)) / (spres_plus_dspres - spres)
    corr = ( entropy(rho,spres**2,mu,ientropy) - S ) / dSdspres
    spres = spres - corr
 enddo
 pres = spres**2

end subroutine get_pres_from_rho_s


!-----------------------------------------------------------------------
!+
!  Calculate intended Gamma_1*dlnrho/dlnp
!+
!-----------------------------------------------------------------------
real function Gamma1_dlnrho_dlnp(delta,r)
 real, intent(in) :: r,delta
 real :: p5,p4,p3,p2,p1,p0,mu,sigma,rcut,rnorm,rrsun

 rcut = 308.9
 rrsun = r/6.969e10

 if (rrsun < rcut) then
    Gamma1_dlnrho_dlnp = 1. + delta
 else
    mu = 326
    sigma = 6.753
    p5 = 0.002137
    p4 = -0.008886
    p3 = -0.0426
    p2 = -0.05158
    p1 = -0.06998
    p0 = 1.0700

    rnorm = (rrsun-mu)/sigma

    Gamma1_dlnrho_dlnp = p5*rnorm**5 + p4*rnorm**4 + p3*rnorm**3 + p2*rnorm**2 + p1*rnorm + p0
    Gamma1_dlnrho_dlnp = max(1.+delta, Gamma1_dlnrho_dlnp)
 endif

 end function Gamma1_dlnrho_dlnp


end module fixedSprofile