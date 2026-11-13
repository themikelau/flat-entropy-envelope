program write_fixed_S_profile
 use fixedSprofile, only:calc_rho_and_pres
 use readwrite_mesa,only:write_mesa,read_mesa
 use eos,           only:calc_temp_and_ene,ieos,X_in,Z_in
 use physcon,       only:solarr
 real, allocatable, dimension(:) :: m,r,rho,pres,ene,temp
 character(len=120)              :: filepath,outputpath
 integer                         :: ierr,i,core_idx
 real, allocatable, dimension(:) :: rhostar,rstar,presstar,mstar,enestar,tempstar,Xfrac,Yfrac,mu
 real                            :: Mtot,rcore

 filepath = './profile40.data'
 call read_mesa(filepath,rho,r,pres,m,ene,temp,X_in,Z_in,Xfrac,Yfrac,mu,Mtot,ierr,cgsunits=.true.)
 call calc_rho_and_pres(m,r,rho,pres,Xfrac,Yfrac,mu)
 print*,'Calculating thermal profile'
 do i = 1,size(m)  ! read_mesa allocates arrays as 1:n
     call calc_temp_and_ene(ieos,rho(i),pres(i),ene(i),temp(i),ierr)
 enddo

 print*,'Flipping arrays to sort from surface to centre'
 call flip_star_arrays(m,pres,temp,r,rho,ene)

 outputpath = 'fixedSprofile.dat'
 call write_mesa(outputpath,m,pres,temp,r,rho,ene)
 print*,'Profile written to ',outputpath

 print*,'May this envelope stay put.'

contains


!-----------------------------------------------------------------------
!+
!  Flip arrays
!+
!-----------------------------------------------------------------------
subroutine flip_star_arrays(m,pres,temp,r,rho,ene)
 use table_utils,   only:flip_array
 real, intent(inout), dimension(:) :: m,pres,temp,r,rho,ene

 print*,'Flipping arrays to sort from surface to centre'
 call flip_array(m)
 call flip_array(pres)
 call flip_array(temp)
 call flip_array(r)
 call flip_array(rho)
 call flip_array(ene)

end subroutine flip_star_arrays

end program write_fixed_S_profile