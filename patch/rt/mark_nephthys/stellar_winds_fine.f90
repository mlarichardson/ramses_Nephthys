!################################################################
!################################################################
!################################################################
!################################################################
module stellar_commons
   use amr_commons, ONLY:dp
   integer(kind=4):: nt_SW, nz_SW  ! number of grids for time and metallicities
   real(dp),allocatable,dimension(:):: log_tSW ! log10 yr
   real(dp),allocatable,dimension(:):: log_zSW ! log10 z
   ! Notice that the values below should be normalised to 1Msun
   real(dp),allocatable,dimension(:,:):: log_cmSW  ! cumulative mass fraction 
   real(dp),allocatable,dimension(:,:):: log_ceSW  ! cumulative energy per 1Msun SSP     
   real(dp),allocatable,dimension(:,:):: log_cmzSW  ! cumulative metal mass fraction
   real(dp),allocatable:: log_cmSW_spec(:,:,:)      ! cumulative mass fraction for several species
end module
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_stellar_winds
   use stellar_commons
   use amr_commons
   implicit none
   integer :: iz, ich, i
   real(dp),allocatable,dimension(:,:):: log_cmHSW  ! cumulative H mass fraction 
   real(dp),allocatable,dimension(:,:):: log_cmCSW  ! cumulative C mass fraction 
   real(dp),allocatable,dimension(:,:):: log_cmNSW  ! cumulative N mass fraction 
   real(dp),allocatable,dimension(:,:):: log_cmOSW  ! cumulative O mass fraction 
   real(dp),allocatable,dimension(:,:):: log_cmMgSW ! cumulative Mg mass fraction 
   real(dp),allocatable,dimension(:,:):: log_cmSiSW ! cumulative Si mass fraction 
   real(dp),allocatable,dimension(:,:):: log_cmSSW  ! cumulative S mass fraction 
   real(dp),allocatable,dimension(:,:):: log_cmFeSW ! cumulative Fe mass fraction
   real(kind=8),allocatable:: dum1d(:)
   logical::ok
   if(.not.use_initial_mass)then
      write(*,*) ' ERROR occured in stellar_winds: use_initial_mass=.false.'
      write(*,*) ' set use_initial_mass=.true. or sn2_real_delay=.true.'
      call clean_stop
   endif

   ! Read stellar winds table
   inquire(FILE=TRIM(stellar_winds_file),exist=ok)
   if(.not.ok)then
      if(myid.eq.1)then
         write(*,*)'Cannot access the file ', TRIM(stellar_winds_file)
      endif
      call clean_stop
   endif

   open(unit=10,file=TRIM(stellar_winds_file),status='old',form='unformatted')
   read(10) nt_SW, nz_SW

   allocate(log_tSW  (1:nt_SW))          ! log Gyr
   allocate(log_zSW  (1:nz_SW))          ! log absolute Z
   allocate(log_cmSW (1:nt_SW,1:nz_SW))  ! log cumulative mass fraction per Msun
   allocate(log_ceSW (1:nt_SW,1:nz_SW))  ! log cumulative energy in erg per Msun
   allocate(log_cmzSW(1:nt_SW,1:nz_SW))  ! log cumulative mass fraction per Msun

   if(nchem>0)then
      allocate(log_cmHSW(1:nt_SW,1:nz_SW))
      allocate(log_cmCSW(1:nt_SW,1:nz_SW))
      allocate(log_cmNSW(1:nt_SW,1:nz_SW))
      allocate(log_cmOSW(1:nt_SW,1:nz_SW))
      allocate(log_cmMgSW(1:nt_SW,1:nz_SW))
      allocate(log_cmSiSW(1:nt_SW,1:nz_SW))
      allocate(log_cmSSW(1:nt_SW,1:nz_SW))
      allocate(log_cmFeSW(1:nt_SW,1:nz_SW))
     
      allocate(log_cmSW_spec(1:nchem,1:nt_SW,1:nz_SW))
   endif
   
   allocate(dum1d (1:nt_SW))
   read(10) dum1d
   log_tSW(:) = dum1d(:)
   deallocate(dum1d)
   allocate(dum1d (1:nz_SW))
   read(10) dum1d
   log_zSW(:) = dum1d(:)
   deallocate(dum1d)

   allocate(dum1d (1:nt_SW))
   !  cumulative stellar mass loss
   do iz=1,nz_SW
      read(10) dum1d
      log_cmSW(:,iz) = dum1d(:)
   enddo
   ! cumulative mechanical energy from winds
   do iz=1,nz_SW
      read(10) dum1d
      log_ceSW(:,iz) = dum1d(:)
   enddo
   ! cumulative metal mass from winds
   do iz=1,nz_SW
      read(10) dum1d
      log_cmzSW(:,iz) = dum1d(:)
   enddo

   if(nchem>0)then
      ! cumulative H mass from winds
      do iz=1,nz_SW
         read(10) dum1d
         log_cmHSW(:,iz) = dum1d(:)
      enddo
      ! cumulative He mass from winds
      do iz=1,nz_SW
         read(10) dum1d
         !log_cmHeSW(:,iz) = dum1d(:)
      enddo
      ! cumulative C mass from winds
      do iz=1,nz_SW
         read(10) dum1d
         log_cmCSW(:,iz) = dum1d(:)
      enddo
      ! cumulative N mass from winds
      do iz=1,nz_SW
         read(10) dum1d
         log_cmNSW(:,iz) = dum1d(:)
      enddo
      ! cumulative O mass from winds
      do iz=1,nz_SW
         read(10) dum1d
         log_cmOSW(:,iz) = dum1d(:)
      enddo
      ! cumulative Mg mass from winds
      do iz=1,nz_SW
         read(10) dum1d
         log_cmMgSW(:,iz) = dum1d(:)
      enddo
      ! cumulative Si mass from winds
      do iz=1,nz_SW
         read(10) dum1d
         log_cmSiSW(:,iz) = dum1d(:)
      enddo
      ! cumulative S mass from winds
      do iz=1,nz_SW
         read(10) dum1d
         log_cmSSW(:,iz) = dum1d(:)
      enddo
      ! cumulative Fe mass from winds
      do iz=1,nz_SW
         read(10) dum1d
         log_cmFeSW(:,iz) = dum1d(:)
      enddo

      ich=0
      do i=1,nchem
          ich=ich+1
          if(myid==1) print *, 'SW CHEM(',int(i,kind=2),') = '//TRIM(chem_list(i))
          if(TRIM(chem_list(i))=='H')  log_cmSW_spec(ich,:,:)=log_cmHSW 
          if(TRIM(chem_list(i))=='C')  log_cmSW_spec(ich,:,:)=log_cmCSW 
          if(TRIM(chem_list(i))=='N')  log_cmSW_spec(ich,:,:)=log_cmNSW 
          if(TRIM(chem_list(i))=='O')  log_cmSW_spec(ich,:,:)=log_cmOSW 
          if(TRIM(chem_list(i))=='Mg') log_cmSW_spec(ich,:,:)=log_cmMgSW 
          if(TRIM(chem_list(i))=='Si') log_cmSW_spec(ich,:,:)=log_cmSiSW 
          if(TRIM(chem_list(i))=='S')  log_cmSW_spec(ich,:,:)=log_cmSSW 
          if(TRIM(chem_list(i))=='Fe') log_cmSW_spec(ich,:,:)=log_cmFeSW 
      end do
      
      deallocate(log_cmHSW,log_cmCSW,log_cmNSW,log_cmOSW)
      deallocate(log_cmMgSW,log_cmSiSW,log_cmSSW,log_cmFeSW)
   endif


   deallocate(dum1d)

end subroutine init_stellar_winds
!################################################################
!################################################################
!################################################################
!################################################################
subroutine stellar_winds_fine(ilevel)
  use amr_commons
  use pm_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine computes the thermal energy, the kinetic energy and 
  ! the metal mass dumped in the gas by stars (winds).
  ! This routine is called every fine time step.
  ! NOTICE:
  ! 1) double-check that only stars have tp.ne.0
  ! 2) set use_initial_mass=.true.
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,npart2,icpu
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  ! Gather star particles only.

#if NDIM==3
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        
        ! Count star particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(tp(ipart).ne.0)then ! <-- will remove dark matter and sink particles 
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        ! Gather star particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only star particles
              if(tp(ipart).ne.0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig   
              endif
              if(ip==nvector)then
                 call stellar_winds_dump(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip>0)call stellar_winds_dump(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do 
  ! End loop over cpus

#endif

111 format('   Entering stlelar winds for level ',I2)

end subroutine stellar_winds_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine stellar_winds_dump(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine stellar_winds_fine. Each stellar particle
  ! dumps mass, momentum and energy in the nearest grid cell using array
  ! unew.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc, ich
  real(dp)::dx_min,vol_min
  real(dp)::dx,dx_loc,scale,birth_time,d
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mloss,mzloss,ethermal,ekinetic,dteff
  real(dp),dimension(1:nvector),save::vol_loc
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp)::msun2g=2d33
  ! stellar library
  real(dp)::mejecta,dfmloss,dfmzloss,log_deloss_erg
  real(dp)::mstar_ini,mstar_ini_msun,zstar
  real(dp)::unit_e_code
  real(dp)::dfmloss_spec(1:nchem)
  real(dp),dimension(1:nchem,1:nvector)::mloss_spec
  ! fractional abundances ; for ionisation fraction and ref, etc 
  real(dp),dimension(1:nvector,1:NVAR),save::fractions ! not compatible with delayed cooling
  integer::ivar,i_fractions

  ! starting index for passive variables except for imetal and chem
  i_fractions = imetal+nchem+1  

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc(1:nvector)=dx_loc**ndim
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

#if NDIM==3
  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/dx
     end do
  end do

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=int(x(j,idim))
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igd(j,idim)=id(j,idim)/2
     end do
  end do
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
  do j=1,np
     igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        else
           icd(j,idim)=id(j,idim)/3 ! 0 or 1
           igrid(j) = ind_grid(ind_grid_part(j))
           ok(j) = .true. ! moving to it's last known cell ... better than nothing. 
           ! So never injecting into coarser level
        end if
     end do
  end do
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do

  ! Compute parent cell adresses
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     else
        indp(j) = nbors_father_cells(ind_grid_part(j),kg(j))
        vol_loc(j)=vol_loc(j)*2**ndim ! ilevel-1 cell volume
     end if
  end do

  ! Compute individual time steps
  do j=1,np
     dteff(j)=dtnew(levelp(ind_part(j)))
  end do

  if(use_proper_time)then
     do j=1,np
        dteff(j)=dteff(j)*aexp**2
     end do
  endif

  ! Reset ejected mass, metallicity, thermal energy
  do j=1,np
     mloss(j)=0d0
     mzloss(j)=0d0
     ethermal(j)=0d0
     mloss_spec(:,j)=0d0
  end do

  ! Compute stellar mass loss and thermal feedback due to stellar winds
  do j=1,np

     ! properties of a star particle
     birth_time=tp(ind_part(j))
     zstar = zp(ind_part(j))
     mstar_ini = mp0(ind_part(j)) ! use_initial_mass=.true. 
     mstar_ini_msun = mstar_ini*(scale_l**3*scale_d/msun2g)

     call cmp_stellar_wind_props (birth_time,dteff(j),zstar,dfmloss,log_deloss_erg,dfmzloss,dfmloss_spec)

     ! Stellar mass loss
     mejecta= mstar_ini*dfmloss  ! dfmloss: mass loss fraction during dteff(j)
     mloss(j)=mloss(j)+mejecta/vol_loc(j)
     ! Thermal energy
     unit_e_code = mstar_ini_msun*(10d0**dble(log_deloss_erg)/dble(msun2g)/scale_v**2)
     ethermal(j)=ethermal(j)+unit_e_code*(mejecta/vol_loc(j))
     ! Metallicity
     if(metal)then
        mzloss(j)=mzloss(j)+mstar_ini*dfmzloss/vol_loc(j)
     endif
     ! Chemical species
     if(nchem>0)then
        do ich=1,nchem
           mloss_spec(ich,j)=mloss_spec(ich,j)+mstar_ini*dfmloss_spec(ich)/vol_loc(j)
        end do 
     endif
     ! Reduce star particle mass
     mp(ind_part(j))=mp(ind_part(j))-mejecta
  end do

  ! NOTE: This assumes that anything with ivar>imetal+nchem is a fractional quantity
  ! store ionisation fractions, ref, etc.
  if(i_fractions.le.nvar)then
     do j=1,np
        d = unew(indp(j),1)
        if (d .gt. 0.d0) then
          do ivar=i_fractions,nvar
             fractions(j,ivar) = unew(indp(j),ivar)/d
          end do
        else
          write(6,*) "Error: getting d is <= 0"
          write(6,*) "myid, igrid", myid, ind_grid(ind_grid_part(j))
          write(6,*) " j_part = ", j
          write(6,*) " id_part = ", ind_part(j), idp(ind_part(j))
          write(6,*) " x_part = ", xp(ind_part(j),1), xp(ind_part(j),2), xp(ind_part(j),3)
          write(6,*) " indp = ", indp(j)
          write(6,*) " d(indp) = ", d
          write(6,*) " igrid, icell, ind_grid_part = ", igrid(j), icell(j),ind_grid_part(j)
          write(6,*) " ingb, fath_grid, fath_cpu = ", (i, nbors_father_grids(ind_grid_part(j),i), &
                                cpu_map(nbors_father_grids(ind_grid_part(j),i)),i=1,8)
          write(6,*) " id,igd,idc1 = ", id(j,1),igd(j,1),icd(j,1)
          write(6,*) " id,igd,idc2 = ", id(j,2),igd(j,2),icd(j,2)
          write(6,*) " id,igd,idc3 = ", id(j,3),igd(j,3),icd(j,3)
          write(6,*) " kg = ", kg(j)
          write(6,*) " mass_loss = ", mloss(j)
          stop 
        endif
     end do
  endif
 
  ! Update hydro variables due to feedback
  do j=1,np
     ! Specific kinetic energy of the star
     ekinetic(j)=0.5*(vp(ind_part(j),1)**2 &
          &          +vp(ind_part(j),2)**2 &
          &          +vp(ind_part(j),3)**2)

     ! Update hydro variable in NGP cell
     unew(indp(j),1)=unew(indp(j),1)+mloss(j)
     unew(indp(j),2)=unew(indp(j),2)+mloss(j)*vp(ind_part(j),1)
     unew(indp(j),3)=unew(indp(j),3)+mloss(j)*vp(ind_part(j),2)
     unew(indp(j),4)=unew(indp(j),4)+mloss(j)*vp(ind_part(j),3)
     unew(indp(j),5)=unew(indp(j),5)+mloss(j)*ekinetic(j)+ ethermal(j)
      
  end do

  ! Add metals
  if(metal)then
     do j=1,np
        unew(indp(j),imetal)=unew(indp(j),imetal)+mzloss(j)
     end do
  endif

  ! Add individual species
  if(nchem>0)then
     do ich=1,nchem
        do j=1,np
           unew(indp(j),ichem+ich-1)=unew(indp(j),ichem+ich-1)+mloss_spec(ich,j)
        end do
     end do
  endif

  ! NOTE: This assumes that anything with ivar>imetal+nchem is a fractional quantity
  ! Update ionisation fractions, ref, etc.
  if(i_fractions.le.nvar)then
     do j=1,np
        d = unew(indp(j),1)
        do ivar=i_fractions,nvar
           unew(indp(j),ivar) = d * fractions(j,ivar)
        end do
     end do
  endif 

#endif
  
end subroutine stellar_winds_dump
!################################################################
!################################################################
!################################################################
!################################################################
subroutine cmp_stellar_wind_props (birth_time,dteff, zstar,dfmloss, log_deloss_erg, dfmzloss, dfmloss_spec)
   use amr_commons
   use stellar_commons
   implicit none
   real(dp),intent(in)::birth_time, dteff, zstar
   real(dp)::dfmloss, dfmzloss, log_deloss_erg
   real(dp),dimension(1:nchem)::dfmloss_spec
   real(dp)::age1, age2, log_age1,log_age2,log_met
   real(dp)::ft1, ft2, fz
   real(dp)::cm1,cm2,cmz1,cmz2,ce1,ce2,dum1,dum2
   integer:: itg1, itg2, izg, ich

   ! initialise
   dfmloss = 0d0
   dfmzloss = 0d0
   log_deloss_erg = -99.

   ! convert the time to physical units
   call getStarAgeGyr(birth_time+dteff, age1)
   call getStarAgeGyr(birth_time      , age2) ! double-checked.

   log_age1    = log10(max(age1*1d9,1.d0))
   log_age2    = log10(max(age2*1d9,1.d0))
   log_met     = log10(max(zstar,z_ave*0.02))

   ! search for the time index from stellar winds library
   call binary_search(log_tSW, log_age1, nt_SW, itg1)
   call binary_search(log_tSW, log_age2, nt_SW, itg2)

   ! search for the metallicity index from stellar winds library
   call binary_search(log_zSW, log_met , nz_SW, izg )

   ! find where we are
   ft1 = (log_tSW(itg1+1) - log_age1)/(log_tSW(itg1+1)-log_tSW(itg1))
   ft2 = (log_tSW(itg2+1) - log_age2)/(log_tSW(itg2+1)-log_tSW(itg2))
   fz  = (log_zSW(izg +1) - log_met )/(log_zSW(izg +1)-log_zSW(izg ))

   ! no extrapolation
   if (ft1 < 0.0) ft1 = 0.0 
   if (ft1 > 1.0) ft1 = 1.0 
   if (ft2 < 0.0) ft2 = 0.0 
   if (ft2 > 1.0) ft2 = 1.0 
   if (fz  < 0.0) fz  = 0.0 
   if (fz  > 1.0) fz  = 1.0 

   ! if a star particle is younger than log_tSW(1), no mass loss 
   if(itg2.eq.1.and.ft2>0.999) return

   ! mass loss fraction during [birth_time, birth_time+dteff]
   dum1 = log_cmSW(itg1,izg  )*ft1 + log_cmSW(itg1+1,izg  )*(1d0-ft1)
   dum2 = log_cmSW(itg1,izg+1)*ft1 + log_cmSW(itg1+1,izg+1)*(1d0-ft1)
   cm1  = dum1*fz + dum2*(1d0-fz)
   dum1 = log_cmSW(itg2,izg  )*ft2 + log_cmSW(itg2+1,izg  )*(1d0-ft2)
   dum2 = log_cmSW(itg2,izg+1)*ft2 + log_cmSW(itg2+1,izg+1)*(1d0-ft2)
   cm2  = dum1*fz + dum2*(1d0-fz)
   dfmloss  = 10d0**cm2  - 10d0**cm1

   ! metal mass loss fraction during [birth_time, birth_time+dteff]
   dum1 = log_cmzSW(itg1,izg  )*ft1 + log_cmzSW(itg1+1,izg  )*(1d0-ft1)
   dum2 = log_cmzSW(itg1,izg+1)*ft1 + log_cmzSW(itg1+1,izg+1)*(1d0-ft1)
   cmz1 = dum1*fz + dum2*(1d0-fz)
   dum1 = log_cmzSW(itg2,izg  )*ft2 + log_cmzSW(itg2+1,izg  )*(1d0-ft2)
   dum2 = log_cmzSW(itg2,izg+1)*ft2 + log_cmzSW(itg2+1,izg+1)*(1d0-ft2)
   cmz2 = dum1*fz + dum2*(1d0-fz)
   dfmzloss = 10d0**cmz2 - 10d0**cmz1

   ! energy during [birth_time, birth_time+dteff]
   dum1 = log_ceSW(itg1,izg  )*ft1 + log_ceSW(itg1+1,izg  )*(1d0-ft1)
   dum2 = log_ceSW(itg1,izg+1)*ft1 + log_ceSW(itg1+1,izg+1)*(1d0-ft1)
   ce1  = dum1*fz + dum2*(1d0-fz)
   dum1 = log_ceSW(itg2,izg  )*ft2 + log_ceSW(itg2+1,izg  )*(1d0-ft2)
   dum2 = log_ceSW(itg2,izg+1)*ft2 + log_ceSW(itg2+1,izg+1)*(1d0-ft2)
   ce2  = dum1*fz + dum2*(1d0-fz)
   log_deloss_erg = log10(10d0**dble(ce2) - 10d0**dble(ce1) + 1d-50) 

   do ich=1,nchem
      ! mass loss fraction during [birth_time, birth_time+dteff]
      dum1 = log_cmSW_spec(ich,itg1,izg  )*ft1 + log_cmSW_spec(ich,itg1+1,izg  )*(1d0-ft1)
      dum2 = log_cmSW_spec(ich,itg1,izg+1)*ft1 + log_cmSW_spec(ich,itg1+1,izg+1)*(1d0-ft1)
      cm1  = dum1*fz + dum2*(1d0-fz)
      dum1 = log_cmSW_spec(ich,itg2,izg  )*ft2 + log_cmSW_spec(ich,itg2+1,izg  )*(1d0-ft2)
      dum2 = log_cmSW_spec(ich,itg2,izg+1)*ft2 + log_cmSW_spec(ich,itg2+1,izg+1)*(1d0-ft2)
      cm2  = dum1*fz + dum2*(1d0-fz)
      dfmloss_spec(ich) = 10d0**cm2  - 10d0**cm1
   end do

   ! For Chabrier IMF
   if(mass_loss_boost>1)then
      dfmloss = dfmloss * mass_loss_boost
      dfmzloss = dfmloss * mass_loss_boost
      do ich=1,nchem
         dfmloss_spec(ich) = dfmloss_spec(ich) * mass_loss_boost
      end do
   endif
  
!777 format(f5.2,1x,f5.2,1x,5(e15.7,1x))
!   write(*,777) log_age1, log_age2, dfmloss, dfmzloss, log_deloss_erg,ce1,ce2

end subroutine cmp_stellar_wind_props
!################################################################
!################################################################
!################################################################
!################################################################
subroutine binary_search(database,xtarget,ndata,i)
   use amr_commons,ONLY:dp
   implicit none 
   integer::i,j,k
   integer,intent(in)::ndata
   real(dp),intent(in)::database(1:ndata),xtarget

   i=1  
   j=ndata
   do   
     k=(i+j)/2
     if (xtarget<database(k)) then 
         j=k  
     else 
         i=k  
     end if
     if (i+1>=j) exit 
   end do

end subroutine binary_search


