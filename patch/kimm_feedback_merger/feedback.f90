!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! Leo
! The thermal feedback part should be very close to the main version
! (../../pm/feedback.f90)
! Only enter in thermal feedback if f_w <= 0 then id>0 (id for stars always positive)
!
! The kinetic feedback part is different (compared to ../../pm/feedback.f90), 
! since the mass-loading is done on the cell where SN explodes (no more debris particles)
! this is the kinetic feedback a la Yohan
! Only enter in kinetic feedback if f_w > 0 then if id<0 
!
! SF with id<0
!
! Be careful that some parameters/values are still hardcoded in the various subroutines below...
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Jokis SD (stellar density) patch:
! Added  st_n_tp and st_n_sn variables (oct 2013)

! Additional stuff by joki for outputting SNe statistics
module SN_stats
  use amr_parameters, ONLY: dp
  real(dp)::E_SN_want=0             ! Prescribed total SNe energy
  real(dp)::E_SN_got=0              ! Injected total SNe energy
  integer::n_SN_want=0
  integer::n_SN_stoch_got=0
  real(dp)::pr_tot=0                ! Summed probability
  integer::n_pr=0                   ! Number of probabilities taken
end module SN_stats

!################################################################
!################################################################
!################################################################
!################################################################
subroutine thermal_feedback(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine computes the thermal energy, the kinetic energy and 
  ! the metal mass dumped in the gas by stars (SNII, SNIa, winds).
  ! This routine is called every fine time step.
  !------------------------------------------------------------------------
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::t0,scale,dx_min,vsn,rdebris,ethermal
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc
  real(dp),dimension(1:3)::skip_loc
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
              if(idp(ipart).lt.0.and.tp(ipart).ne.0)then
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
              if(idp(ipart).lt.0.and.tp(ipart).ne.0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig   
              endif
              if(ip==nvector)then
                 call feedbk(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
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
     if(ip>0)call feedbk(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do 
  ! End loop over cpus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(ilevel .eq. levelmin) call output_SN_stats ! Once per coarse step    !Stoch FB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#endif

111 format('   Entering thermal_feedback for level ',I2)

end subroutine thermal_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine feedbk(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use random
  use SN_stats                                                             !Stoch FB
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine feedback. Each stellar particle
  ! dumps mass, momentum and energy in the nearest grid cell using array
  ! unew.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc
  real(kind=8)::RandNum
  real(dp)::SN_BOOST,mstar,dx_min,vol_min
  real(dp)::xxx,mmm,t0,ESN,mejecta,zloss
  real(dp)::ERAD,RAD_BOOST,tauIR,eta_sig
  real(dp)::sigma_d,delta_x,tau_factor,rad_factor
  real(dp)::dx,dx_loc,scale,birth_time,current_time
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  integer,dimension(1:nvector),save::igrid_son,ind_son
  integer,dimension(1:nvector),save::list1
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mloss,mzloss,ethermal,ekinetic,dteff
  real(dp),dimension(1:nvector),save::vol_loc
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp)::prob,de_min,mcell                                               !DSfeedb

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

  ! Minimum star particle mass
  if(m_star < 0d0)then
     mstar=n_star/(scale_nH*aexp**3)*vol_min
  else
     mstar=m_star*mass_sph
  endif

  ! Compute stochastic boost to account for target GMC mass
  SN_BOOST=MAX(mass_gmc*2d33/(scale_d*scale_l**3)/mstar,1d0)

  ! Massive star lifetime from Myr to code units
  if(use_proper_time)then
     t0=t_delay*1d6*(365.*24.*3600.)/(scale_t/aexp**2)
     current_time=texp
  else
     t0=t_delay*1d6*(365.*24.*3600.)/scale_t
     current_time=t
  endif

  ! Type II supernova specific energy from cgs to code units
  ESN=1d51/(10.*2d33)/scale_v**2

  ! Life time radiation specific energy from cgs to code units
  ERAD=1d53/(10.*2d33)/scale_v**2

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
        id(j,idim)=x(j,idim)
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
  end do

  ! Compute stellar mass loss and thermal feedback due to supernovae
  do j=1,np
     birth_time=tp(ind_part(j))
     ! Make sure that we don't count feedback twice
     if(birth_time.le.(current_time-t0))then           
        ! Stellar mass loss
        mejecta=eta_sn*mp(ind_part(j))
        mloss(j)=mloss(j)+mejecta/vol_loc(j)
        ! Thermal energy
        ethermal(j)=ethermal(j)+mejecta*ESN/vol_loc(j)
        ! Metallicity
        if(metal)then
           zloss=yield+(1d0-yield)*zp(ind_part(j))
           mzloss(j)=mzloss(j)+mejecta*zloss/vol_loc(j)
        endif
        ! Reduce star particle mass
        mp(ind_part(j))=mp(ind_part(j))-mejecta
        idp(ind_part(j))=-idp(ind_part(j))
        ! Update local density at SN event:
        st_n_SN(ind_part(j)) = uold(indp(j),1)       !-----------------!SD
        ! Boost SNII energy and depopulate accordingly
        if(SN_BOOST>1d0)then
           call ranf(localseed,RandNum)
           if(RandNum<1d0/SN_BOOST)then
              mloss(j)=SN_BOOST*mloss(j)
              mzloss(j)=SN_BOOST*mzloss(j)
              ethermal(j)=SN_BOOST*ethermal(j)
           else
              mloss(j)=0d0
              mzloss(j)=0d0
              ethermal(j)=0d0
           endif
        endif

        !BEGIN Dalla Vecchia & Schaye feedback----------------------------
        n_SN_want = n_SN_want + 1                             ! Statistics
        E_SN_want = E_SN_want + mejecta*ESN                   ! Statistics
        if(SN_dT2_min .gt. 0d0) then
           de_min = SN_dT2_min / (gamma-1d0) / scale_T2
           mcell = unew(indp(j),1) * vol_loc(j)
           prob = mejecta*ESN / de_min / (mcell+mejecta)  ! SN probability
           if(prob.le.1d0) then     ! Stoch SNe to fill in required deltaT
              n_pr = n_pr + 1                                 ! Statistics
              pr_tot = pr_tot + prob                          ! Statistics
              call ranf(localseed,RandNum)
              if(RandNum<prob) then
                 ethermal(j) = de_min * (mcell+mejecta) / vol_loc(j)
                 n_SN_stoch_got = n_SN_stoch_got+1            ! Statistics
              else
                 ethermal(j)=0d0
              endif
           endif
        endif
        E_SN_got = E_SN_got + ethermal(j) * vol_loc(j)        ! Statistics
        !END Dalla Vecchia & Schaye feedback------------------------------
        if(ethermal(j) .gt. 0d0) then 
           st_n_SN(ind_part(j)) = unew(indp(j),1)       !--------------!SD
           st_e_SN(ind_part(j)) = ethermal(j) * vol_loc(j) * scale_d  &!SD
                * scale_l**3 * scale_v**2 / 1d51        !--------------!SD
        endif

     endif
  end do

  ! Update hydro variables due to feedback

  ! For IR radiation trapping,
  ! we use a fixed length to estimate the column density of gas
  delta_x=200.*3d18
  if(metal)then
     tau_factor=kappa_IR*delta_x*scale_d/0.02
  else
     tau_factor=kappa_IR*delta_x*scale_d*z_ave
  endif
  rad_factor=ERAD/ESN
  do j=1,np

     ! Infrared photon trapping boost
     if(metal)then
        tauIR=tau_factor*max(uold(indp(j),imetal),smallr)
     else
        tauIR=tau_factor*max(uold(indp(j),1),smallr)
     endif
     if(uold(indp(j),1)*scale_nH > 10.)then
        RAD_BOOST=rad_factor*(1d0-exp(-tauIR))
     else
        RAD_BOOST=0.0
     endif

     ! Specific kinetic energy of the star
     ekinetic(j)=0.5*(vp(ind_part(j),1)**2 &
          &          +vp(ind_part(j),2)**2 &
          &          +vp(ind_part(j),3)**2)

     ! Update hydro variable in NGP cell
     unew(indp(j),1)=unew(indp(j),1)+mloss(j)
     unew(indp(j),2)=unew(indp(j),2)+mloss(j)*vp(ind_part(j),1)
     unew(indp(j),3)=unew(indp(j),3)+mloss(j)*vp(ind_part(j),2)
     unew(indp(j),4)=unew(indp(j),4)+mloss(j)*vp(ind_part(j),3)
     unew(indp(j),5)=unew(indp(j),5)+mloss(j)*ekinetic(j)+ &
          & ethermal(j)*(1d0+RAD_BOOST)
          
  end do

  ! Add metals
  if(metal)then
     do j=1,np
        unew(indp(j),imetal)=unew(indp(j),imetal)+mzloss(j)
     end do
  endif

  ! Add delayed cooling switch variable
  if(delayed_cooling)then
     do j=1,np
        unew(indp(j),idelay)=unew(indp(j),idelay)+mloss(j)
     end do
  endif

#endif
  
end subroutine feedbk

subroutine kinetic_feedback
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::nSN_tot_all
  integer,dimension(1:ncpu)::nSN_icpu_all
  real(dp),dimension(:),allocatable::mSN_all,ZSN_all
  real(dp),dimension(:,:),allocatable::xSN_all,vSN_all
#endif
  !----------------------------------------------------------------------
  ! Description: This subroutine checks SN events in cells where a
  ! star particle has been spawned.
  ! Yohan Dubois
  !----------------------------------------------------------------------
  ! local constants
  integer::icpu,igrid,jgrid,npart1,npart2,ipart,jpart,next_part
  integer::nSN,nSN_tot,info,iSN,ilevel,ivar
  integer,dimension(1:ncpu)::nSN_icpu
  logical ::ok_free
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,t0,mass_load
  integer ,dimension(:),allocatable::indSN,iSN_myid
  real(dp),dimension(:),allocatable::mSN,ZSN,vol_gas,ekBlast,mloadSN,ZloadSN
  real(dp),dimension(:,:),allocatable::xSN,vSN,dq,vloadSN
  integer ,dimension(:),allocatable::indSN_tot,itemp
  real(dp),dimension(:),allocatable::mSN_tot,ZSN_tot
  real(dp),dimension(:,:),allocatable::xSN_tot,vSN_tot
  integer::isort

  if(.not. hydro)return
  if(ndim.ne.3)return

  if(verbose)write(*,*)'Entering kinetic_feedback'
  
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Time delay for SN explosion from Myr to code units
  t0=t_delay*(1d6*365.*24.*3600.)/scale_t

  !------------------------------------------------------
  ! Gather star particles eligible for a SN event
  !------------------------------------------------------
  nSN_tot=0
  do icpu=1,ncpu
  ! Loop over cpus
     igrid=headl(icpu,levelmin)
     ! Loop over grids
     do jgrid=1,numbl(icpu,levelmin)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0        
        ! Count old enough star particles that have not exploded
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if( idp(ipart).lt.0 .and. tp(ipart).ne.0d0 .and. &
                   & tp(ipart).lt.(t-t0) )then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        nSN_tot=nSN_tot+npart2   ! Add SNe to the total
        igrid=next(igrid)   ! Go to next grid
     end do
  enddo

  nSN_icpu=0
  nSN_icpu(myid)=nSN_tot  
#ifndef WITHOUTMPI
  ! Give an array of number of SN on each cpu available to all cpus
  call MPI_ALLREDUCE(nSN_icpu,nSN_icpu_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nSN_icpu=nSN_icpu_all
#endif
  nSN_tot=sum(nSN_icpu(1:ncpu))

  if(myid==1)then
     write(*,*)'-----------------------------------------------'
     write(*,*)'Number of SN to explode=',nSN_tot
     write(*,*)'-----------------------------------------------'
  endif

  if (nSN_tot .eq. 0) return
  
  ! Allocate arrays for the position and the mass of the SN
  allocate(xSN_tot(1:nSN_tot,1:3),vSN_tot(1:nSN_tot,1:3),mSN_tot(1:nSN_tot),ZSN_tot(1:nSN_tot),itemp(1:nSN_tot))

  xSN_tot=0.;vSN_tot=0.;mSN_tot=0.;ZSN_tot=0.
  !------------------------------------------------------
  ! Give position and mass of the star to the SN array
  !------------------------------------------------------
  if(myid==1)then
     iSN=0
  else
     iSN=sum(nSN_icpu(1:myid-1))
  endif
  do icpu=1,ncpu
  ! Loop over cpus
     igrid=headl(icpu,levelmin)
     ! Loop over grids
     do jgrid=1,numbl(icpu,levelmin)
        npart1=numbp(igrid)  ! Number of particles in the grid        
        ! Count old enough star particles that have not exploded
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if( idp(ipart).lt.0 .and. tp(ipart).ne.0d0 .and. &
                   & tp(ipart).lt.(t-t0) )then
                 iSN=iSN+1
                 xSN_tot(iSN,1)=xp(ipart,1)
                 xSN_tot(iSN,2)=xp(ipart,2)
                 xSN_tot(iSN,3)=xp(ipart,3)
                 vSN_tot(iSN,1)=vp(ipart,1)
                 vSN_tot(iSN,2)=vp(ipart,2)
                 vSN_tot(iSN,3)=vp(ipart,3)
                 mSN_tot(iSN)  =eta_sn*mp(ipart)
                 if(metal)ZSN_tot(iSN)=zp(ipart)
                 ! Remove the mass ejected by the SN
                 mp(ipart) =mp(ipart)-eta_sn*mp(ipart)
                 idp(ipart)=-idp(ipart)
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
     
        igrid=next(igrid)   ! Go to next grid
     end do
  enddo

#ifndef WITHOUTMPI
  allocate(xSN_all(1:nSN_tot,1:3),vSN_all(1:nSN_tot,1:3),mSN_all(1:nSN_tot),ZSN_all(1:nSN_tot))
  call MPI_ALLREDUCE(xSN_tot,xSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vSN_tot,vSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mSN_tot,mSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ZSN_tot,ZSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  xSN_tot=xSN_all
  vSN_tot=vSN_all
  mSN_tot=mSN_all
  ZSN_tot=ZSN_all
  deallocate(xSN_all,vSN_all,mSN_all,ZSN_all)
#endif

  call getSNonmyid(itemp,nSN,xSN_tot,nSN_tot)

  ! Allocate the arrays for the position and the mass of the SN
  allocate(xSN(1:nSN,1:3),vSN(1:nSN,1:3),mSN(1:nSN),ZSN(1:nSN),iSN_myid(1:nSN))
  xSN=0d0; vSN=0d0; mSN=0d0; ZSN=0d0; iSN_myid=0

  do iSN=1,nSN
     isort=itemp(iSN)
     iSN_myid(iSN)=isort
     xSN(iSN,1)=xSN_tot(isort,1)
     xSN(iSN,2)=xSN_tot(isort,2)
     xSN(iSN,3)=xSN_tot(isort,3)
     vSN(iSN,1)=vSN_tot(isort,1)
     vSN(iSN,2)=vSN_tot(isort,2)
     vSN(iSN,3)=vSN_tot(isort,3)
     mSN(iSN)  =mSN_tot(isort)
     ZSN(iSN)  =ZSN_tot(isort)
  enddo
  deallocate(xSN_tot,vSN_tot,mSN_tot,ZSN_tot,itemp)

  allocate(vol_gas(1:nSN),dq(1:nSN,1:3),ekBlast(1:nSN),indSN(1:nSN))
  allocate(mloadSN(1:nSN),ZloadSN(1:nSN),vloadSN(1:nSN,1:3))

  ! Compute the grid discretization effects
  call average_SN(xSN,vSN,vol_gas,dq,ekBlast,indSN,nSN,nSN_tot,iSN_myid,mSN,mloadSN,ZSN,ZloadSN,vloadSN)

  ! Modify hydro quantities to account for a Sedov blast wave
  call Sedov_blast(xSN,mSN,indSN,vol_gas,dq,ekBlast,nSN,mloadSN,ZloadSN,vloadSN)

  deallocate(xSN,vSN,mSN,ZSN,iSN_myid)
  deallocate(indSN,vol_gas,dq,ekBlast)
  deallocate(mloadSN,ZloadSN,vloadSN)

  ! Update hydro quantities for split cells
  do ilevel=nlevelmax,levelmin,-1
     call upload_fine(ilevel)
     do ivar=1,nvar
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     enddo
  enddo

end subroutine kinetic_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine average_SN(xSN,vSN,vol_gas,dq,ekBlast,ind_blast,nSN,nSN_tot,iSN_myid,mSN,mloadSN,ZSN,ZloadSN,vloadSN)
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine average the hydro quantities inside the SN bubble
  ! and do the mass loading process
  !------------------------------------------------------------------------
  integer::ilevel,ncache,nSN,nSN_tot,j,iSN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,info
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dr_SN,d,u,v,w,ek,u2,v2,w2,dr_cell
  real(dp)::scale,dx,dxx,dyy,dzz,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::eint,ekk,ekk1,ekk2,heat,mload,Zload
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ,dimension(1:nSN)::ind_blast
  real(dp),dimension(1:nSN)::mSN,mloadSN,ZSN,ZloadSN,vol_gas,ekBlast
  real(dp),dimension(1:nSN,1:3)::xSN,vSN,dq,u2Blast,vloadSN
#ifndef WITHOUTMPI
  real(dp),dimension(1:nSN_tot)::vol_gas_mpi
  real(dp),dimension(1:nSN_tot)::vol_gas_all,ekBlast_all
  real(dp),dimension(1:nSN_tot)::mloadSN_mpi,mloadSN_all,ZloadSN_mpi,ZloadSN_all
  real(dp),dimension(1:nSN_tot,1:3)::dq_mpi,u2Blast_mpi,vloadSN_mpi
  real(dp),dimension(1:nSN_tot,1:3)::dq_all,u2Blast_all,vloadSN_all
#endif
  logical ,dimension(1:nvector),save::ok
  integer ,dimension(1:nSN)::iSN_myid
  integer::ind_SN
  real(dp)::rcell
  rcell=2.0d0 !Leo: bubble radius

  if(verbose)write(*,*)'Entering average_SN'

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  rmax=MAX(rcell*dx_min*scale_l/aexp,rbubble*3.08d18)
  !!write(*,*)'[debug] rmax = ',rmax,rcell,rbubble,scale_l,dx_min
  rmax=rmax/scale_l
  rmax2=rmax*rmax
  
  ! Initialize the averaged variables
  vol_gas=0.0;dq=0.0;u2Blast=0.0;ekBlast=0.0;ind_blast=-1;mloadSN=0.0;ZloadSN=0.0;vloadSN=0.0

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities) 
     dx=0.5D0**ilevel 
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim  
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 do iSN=1,nSN
                    ! Check if the cell lies within the SN radius
                    dxx=x-xSN(iSN,1)
                    dyy=y-xSN(iSN,2)
                    dzz=z-xSN(iSN,3)
                    dr_SN=dxx**2+dyy**2+dzz**2
                    dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
                    if(dr_SN.lt.rmax2)then
                       !!write(*,*)'[debug] i am in = ',iSN,dr_SN,dr_cell,rmax2,ilevel
                       vol_gas(iSN)=vol_gas(iSN)+vol_loc
                       ! Take account for grid effects on the conservation of the
                       ! normalized linear momentum
                       u=dxx/rmax
                       v=dyy/rmax
                       w=dzz/rmax
                       ! Add the local normalized linear momentum to the total linear
                       ! momentum of the blast wave (should be zero with no grid effect)
                       dq(iSN,1)=dq(iSN,1)+u*vol_loc
                       dq(iSN,2)=dq(iSN,2)+v*vol_loc
                       dq(iSN,3)=dq(iSN,3)+w*vol_loc
                       u2Blast(iSN,1)=u2Blast(iSN,1)+u*u*vol_loc
                       u2Blast(iSN,2)=u2Blast(iSN,2)+v*v*vol_loc
                       u2Blast(iSN,3)=u2Blast(iSN,3)+w*w*vol_loc
                    endif
                    if(dr_cell.le.dx_loc/2.0)then
                       ! Leo : debug
                       !!if(dr_SN.lt.rmax2)then
                       !!   write(*,*)'[debug] both conditions ok ',iSN,dr_cell,dx_loc/2.0,ilevel
                       !!else
                       !!   write(*,*)'[debug] oh oh SN in the cell but cell not in the bubble...',iSN,dr_cell,dx_loc/2.0,ilevel,dr_SN,rmax2
                       !!endif
                       ind_blast(iSN)=ind_cell(i)
                       ekBlast  (iSN)=vol_loc
                       d=uold(ind_blast(iSN),1)
                       u=uold(ind_blast(iSN),2)/d
                       v=uold(ind_blast(iSN),3)/d
                       w=uold(ind_blast(iSN),4)/d
                       ekk=0.5d0*d*(u*u+v*v+w*w)
                       eint=uold(ind_blast(iSN),5)-ekk
                       ! Mass loading factor of the Sedov explosion
                       ! Ensure that no more that 25% of the gas content is removed
                       mload=min(f_w*mSN(iSN),0.25d0*d*vol_loc)
                       mloadSN(iSN)=mSN(iSN)+mload
                       ! Update gas mass and metal content in the cell
                       if(metal)then
                          Zload=uold(ind_blast(iSN),imetal)/d
                          ZloadSN(iSN)=( mload*Zload + yield*(1-zSN(iSN))*mSN(iSN) + zSN(iSN)*mSN(iSN)) / mloadSN(iSN)
                          uold(ind_blast(iSN),imetal)=uold(ind_blast(iSN),imetal)-Zload*mload/vol_loc
                       endif
                       d=uold(ind_blast(iSN),1)-mload/vol_loc

                       uold(ind_blast(iSN),1)=d
                       uold(ind_blast(iSN),2)=d*u
                       uold(ind_blast(iSN),3)=d*v
                       uold(ind_blast(iSN),4)=d*w
                       uold(ind_blast(iSN),5)=eint+0.5d0*d*(u*u+v*v+w*w)

                       vloadSN(iSN,1)=(mSN(iSN)*vSN(iSN,1)+mload*u)/mloadSN(iSN)
                       vloadSN(iSN,2)=(mSN(iSN)*vSN(iSN,2)+mload*v)/mloadSN(iSN)
                       vloadSN(iSN,3)=(mSN(iSN)*vSN(iSN,3)+mload*w)/mloadSN(iSN)
                    endif
                 end do
              endif
           end do
           
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

  !################################################################
#ifndef WITHOUTMPI
  vol_gas_mpi=0d0; dq_mpi=0d0; u2Blast_mpi=0d0; mloadSN_mpi=0d0; ZloadSN_mpi=0d0; vloadSN_mpi=0d0
  ! Put the nSN size arrays into nSN_tot size arrays to synchronize processors
  do iSN=1,nSN
     ind_SN=iSN_myid(iSN)
     vol_gas_mpi(ind_SN)=vol_gas(iSN)
     mloadSN_mpi(ind_SN)=mloadSN(iSN)
     ZloadSN_mpi(ind_SN)=ZloadSN(iSN)
     vloadSN_mpi(ind_SN,1)=vloadSN(iSN,1)
     vloadSN_mpi(ind_SN,2)=vloadSN(iSN,2)
     vloadSN_mpi(ind_SN,3)=vloadSN(iSN,3)
     dq_mpi     (ind_SN,1)=dq     (iSN,1)
     dq_mpi     (ind_SN,2)=dq     (iSN,2)
     dq_mpi     (ind_SN,3)=dq     (iSN,3)
     u2Blast_mpi(ind_SN,1)=u2Blast(iSN,1)
     u2Blast_mpi(ind_SN,2)=u2Blast(iSN,2)
     u2Blast_mpi(ind_SN,3)=u2Blast(iSN,3)
  enddo
  call MPI_ALLREDUCE(vol_gas_mpi,vol_gas_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mloadSN_mpi,mloadSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ZloadSN_mpi,ZloadSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vloadSN_mpi,vloadSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dq_mpi     ,dq_all     ,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(u2Blast_mpi,u2Blast_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  vol_gas_mpi=vol_gas_all
  mloadSN_mpi=mloadSN_all
  ZloadSN_mpi=ZloadSN_all
  vloadSN_mpi=vloadSN_all
  dq_mpi     =dq_all
  u2Blast_mpi=u2Blast_all
  ! Put the nSN_tot size arrays into nSN size arrays
  do iSN=1,nSN
     ind_SN=iSN_myid(iSN)
     vol_gas(iSN)=vol_gas_mpi(ind_SN)
     mloadSN(iSN)=mloadSN_mpi(ind_SN)
     ZloadSN(iSN)=ZloadSN_mpi(ind_SN)
     vloadSN(iSN,1)=vloadSN_mpi(ind_SN,1)
     vloadSN(iSN,2)=vloadSN_mpi(ind_SN,2)
     vloadSN(iSN,3)=vloadSN_mpi(ind_SN,3)
     dq     (iSN,1)=dq_mpi     (ind_SN,1)
     dq     (iSN,2)=dq_mpi     (ind_SN,2)
     dq     (iSN,3)=dq_mpi     (ind_SN,3)
     u2Blast(iSN,1)=u2Blast_mpi(ind_SN,1)
     u2Blast(iSN,2)=u2Blast_mpi(ind_SN,2)
     u2Blast(iSN,3)=u2Blast_mpi(ind_SN,3)
  enddo
#endif
  !################################################################
  !!write(*,*)'[debug] in average_SN entering loop iSN',nSN
  do iSN=1,nSN
     if(vol_gas(iSN)>0d0)then
        dq(iSN,1)=dq(iSN,1)/vol_gas(iSN)
        dq(iSN,2)=dq(iSN,2)/vol_gas(iSN)
        dq(iSN,3)=dq(iSN,3)/vol_gas(iSN)
        u2Blast(iSN,1)=u2Blast(iSN,1)/vol_gas(iSN)
        u2Blast(iSN,2)=u2Blast(iSN,2)/vol_gas(iSN)
        u2Blast(iSN,3)=u2Blast(iSN,3)/vol_gas(iSN)
        u2=u2Blast(iSN,1)-dq(iSN,1)**2
        v2=u2Blast(iSN,2)-dq(iSN,2)**2
        w2=u2Blast(iSN,3)-dq(iSN,3)**2
        !!write(*,*)'[debug] in average_blast',0.5d0*(u2+v2+w2)
        ekBlast(iSN)=max(0.5d0*(u2+v2+w2),0.0d0)
     endif
     !!write(*,*)'[debug] in average_SN',iSN,vol_gas(iSN),ekBlast(iSN)
     !!write(*,*)' '
  end do

  if(verbose)write(*,*)'Exiting average_SN'

end subroutine average_SN

!################################################################
!################################################################
!################################################################
!################################################################
subroutine Sedov_blast(xSN,mSN,indSN,vol_gas,dq,ekBlast,nSN,mloadSN,ZloadSN,vloadSN)
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine merges SN using the FOF algorithm.
  !------------------------------------------------------------------------
  integer::ilevel,j,iSN,nSN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,info,ncache
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dx,dxx,dyy,dzz,dr_SN,d,u,v,w,ek,u_r,d_gas,ESN
  real(dp)::scale,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nSN)::mSN,p_gas,vol_gas,uSedov,ekBlast,mloadSN,ZloadSN
  real(dp),dimension(1:nSN,1:3)::xSN,dq,vloadSN
  integer ,dimension(1:nSN)::indSN
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nSN)::Etot,Ptot, Mtot ! joki debug
  integer(dp),dimension(1:nSN)::nCells ! joki debug
  real(dp)::rcell
  real(dp)::eps_sn
  eps_sn =1.0D0 ! Leo: Supernovae radiation efficiency (1 corresponds to 100%)
  rcell=2.0d0 !Leo: bubble radius

  if(verbose)write(*,*)'Entering Sedov_blast'

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  rmax=MAX(rcell*dx_min*scale_l/aexp,rbubble*3.08d18)
  rmax=rmax/scale_l
  rmax2=rmax*rmax
  
  ! Ejecta specific energy (accounting for dilution)
  ESN=(eps_sn*1d51/(10d0*2d33))/scale_v**2
  ! Etot=0.5*ESN*MSN, Ptot=sqrt(ESN/MSN)
  Etot=0d0 ; Ptot=0d0; Mtot=0d0; nCells=0 !joki debug

  !!!write(*,*)'[debug] in sedov_blast entering loop iSN',nSN
  do iSN=1,nSN
     !!!write(*,*)'[debug]',iSN,vol_gas(iSN),ekBlast(iSN),mloadSN(iSN)
     if(vol_gas(iSN)>0d0)then
        d_gas=mSN(iSN)/vol_gas(iSN)
        if(ekBlast(iSN)==0d0)then
           p_gas (iSN)=d_gas*ESN
           uSedov(iSN)=0d0
        else
           p_gas (iSN)=(1d0-f_ek)*d_gas*ESN
           uSedov(iSN)=sqrt(f_ek*mSN(iSN)*ESN/ekBlast(iSN)/mloadSN(iSN))
        endif
     else
        ! Leo: code crashes when ekBlast = 0. here... not clear when and how it could happen...
        ! write(*,*)'[debug]',iSN,ekBlast(iSN)
        p_gas(iSN)=mSN(iSN)*ESN/ekBlast(iSN)
     endif
  end do

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities) 
     dx=0.5D0**ilevel 
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     ! Cells center position relative to grid center position
     do ind=1,twotondim  
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim  
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 do iSN=1,nSN
                    ! Check if the cell lies within the SN radius
                    dxx=x-xSN(iSN,1)
                    dyy=y-xSN(iSN,2)
                    dzz=z-xSN(iSN,3)
                    dr_SN=dxx**2+dyy**2+dzz**2
                    if(dr_SN.lt.rmax2)then
                       d_gas=mloadSN(iSN)/vol_gas(iSN)
                       ! Compute the density and the metal density of the cell
                       uold(ind_cell(i),1)=uold(ind_cell(i),1)+d_gas
                       if(metal)uold(ind_cell(i),imetal)=uold(ind_cell(i),imetal)+d_gas*ZloadSN(iSN)
                       ! Velocity at a given dr_SN linearly interpolated between zero and uSedov
                       u=uSedov(iSN)*(dxx/rmax-dq(iSN,1))+vloadSN(iSN,1)
                       v=uSedov(iSN)*(dyy/rmax-dq(iSN,2))+vloadSN(iSN,2)
                       w=uSedov(iSN)*(dzz/rmax-dq(iSN,3))+vloadSN(iSN,3)

                       Etot(iSN) = Etot(iSN) + 0.5 * d_gas * vol_loc & !debug joki
                            * ( (u-vloadSN(iSN,1))**2 &
                            + (v-vloadSN(iSN,2))**2   &
                            + (w-vloadSN(iSN,3))**2) 
                       Ptot(iSN)=Ptot(iSN) + d_gas * vol_loc * sqrt( &
                            ( (u-vloadSN(iSN,1))**2 &
                            + (v-vloadSN(iSN,2))**2   &
                            + (w-vloadSN(iSN,3))**2)  )
                       Mtot(iSN)=Mtot(iSN)+d_gas*vol_loc !debug joki
                       nCells(iSN) = nCells(iSN)+1
                       ! Add each momentum component of the blast wave to the gas
                       uold(ind_cell(i),2)=uold(ind_cell(i),2)+d_gas*u
                       uold(ind_cell(i),3)=uold(ind_cell(i),3)+d_gas*v
                       uold(ind_cell(i),4)=uold(ind_cell(i),4)+d_gas*w
                       ! Finally update the total energy of the gas
                       uold(ind_cell(i),5)=uold(ind_cell(i),5)+0.5*d_gas*(u*u+v*v+w*w)+p_gas(iSN)
                    endif
                 end do
              endif
           end do
           
        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

  do iSN=1,nSN
     if(vol_gas(iSN)==0d0)then
        d_gas=mloadSN(iSN)/ekBlast(iSN)
        u=vloadSN(iSN,1)
        v=vloadSN(iSN,2)
        w=vloadSN(iSN,3)
        if(indSN(iSN)>0)then
           uold(indSN(iSN),1)=uold(indSN(iSN),1)+d_gas
           uold(indSN(iSN),2)=uold(indSN(iSN),2)+d_gas*u
           uold(indSN(iSN),3)=uold(indSN(iSN),3)+d_gas*v
           uold(indSN(iSN),4)=uold(indSN(iSN),4)+d_gas*w
           uold(indSN(iSN),5)=uold(indSN(iSN),5)+d_gas*0.5*(u*u+v*v+w*w)+p_gas(iSN)
           if(metal)uold(indSN(iSN),imetal)=uold(indSN(iSN),imetal)+d_gas*ZloadSN(iSN)
        endif
     endif
     !if(ncells(iSN).gt. 0) then
     !   write(*, 113) myid, iSN, nCells(iSN), Etot(iSN), f_ek*ESN*MSN(iSN) &
     !        , Ptot(iSN), sqrt(2d0*f_ek*ESN)*MSN(iSN), Mtot(iSN), MloadSN(iSN)
     !endif
  end do

  if(verbose)write(*,*)'Exiting Sedov_blast'
  ! Etot=0.5*ESN*MSN, Ptot=sqrt(2*ESN*MSN)

113 format('SNk_sstats ',i9, i9, i9, 1pe12.5, 1pe12.5, 1pe12.5, 1pe12.5, 1pe12.5, 1pe12.5)
end subroutine Sedov_blast
!################################################################
!################################################################
!################################################################
!################################################################
subroutine getSNonmyid(iSN_myid,nSN_myid,xSN,nSN)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,dimension(1:nSN)::iSN_myid
  integer::nSN_myid,ii,info,nSN
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  integer::lmin,iSN,nx_loc,ilevel,lmax,bit_length,maxdom,icpu
  integer::imin,jmin,kmin,imax,jmax,kmax,ndom,impi,i,j,k,ncpu_read
  integer,dimension(1:ncpu)::cpu_list
  logical,dimension(1:ncpu)::cpu_read
  real(dp)::scale,dx,dx_min,drSN
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp)::xxmin,yymin,zzmin,xxmax,yymax,zzmax,dmax
  real(qdp),dimension(1:8)::bounding_min,bounding_max
  real(qdp)::dkey,order_min,oneqdp=1.0
  real(dp),dimension(1:nSN,1:3)::xSN
  real(dp)::rcell
  rcell=2.0d0 !Leo: bubble radius
  
  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  drSN=2d0*MAX(rcell*dx_min*scale_l/aexp,rbubble*3.08d18)
  drSN=drSN/scale_l

  !-----------------------
  ! Map parameters
  !-----------------------
  lmax=nlevelmax
  iSN_myid=0
  ii=0
  do iSN=1,nSN

     cpu_read=.false.
     ! Compute boundaries for the SN cube of influence
     xxmin=(xSN(iSN,1)-drSN)/scale ; xxmax=(xSN(iSN,1)+drSN)/scale
     yymin=(xSN(iSN,2)-drSN)/scale ; yymax=(xSN(iSN,2)+drSN)/scale
     zzmin=(xSN(iSN,3)-drSN)/scale ; zzmax=(xSN(iSN,3)+drSN)/scale

     if(TRIM(ordering).eq.'hilbert')then
        
        dmax=max(xxmax-xxmin,yymax-yymin,zzmax-zzmin)
        do ilevel=1,lmax
           dx=0.5d0**ilevel
           if(dx.lt.dmax)exit
        end do
        lmin=ilevel
        bit_length=lmin-1
        maxdom=2**bit_length
        imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
        if(bit_length>0)then
           imin=int(xxmin*dble(maxdom))
           imax=imin+1
           jmin=int(yymin*dble(maxdom))
           jmax=jmin+1
           kmin=int(zzmin*dble(maxdom))
           kmax=kmin+1
        endif
        
        dkey=(real(2**(nlevelmax+1),kind=qdp)/real(maxdom,kind=qdp))**ndim
        ndom=1
        if(bit_length>0)ndom=8
        idom(1)=imin; idom(2)=imax
        idom(3)=imin; idom(4)=imax
        idom(5)=imin; idom(6)=imax
        idom(7)=imin; idom(8)=imax
        jdom(1)=jmin; jdom(2)=jmin
        jdom(3)=jmax; jdom(4)=jmax
        jdom(5)=jmin; jdom(6)=jmin
        jdom(7)=jmax; jdom(8)=jmax
        kdom(1)=kmin; kdom(2)=kmin
        kdom(3)=kmin; kdom(4)=kmin
        kdom(5)=kmax; kdom(6)=kmax
        kdom(7)=kmax; kdom(8)=kmax
        
        do i=1,ndom
           if(bit_length>0)then
              call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
           else
              order_min=0.0d0
           endif
           bounding_min(i)=(order_min)*dkey
           bounding_max(i)=(order_min+oneqdp)*dkey
        end do
        
        cpu_min=0; cpu_max=0
        do impi=1,ncpu
           do i=1,ndom
              if (   bound_key(impi-1).le.bounding_min(i).and.&
                   & bound_key(impi  ).gt.bounding_min(i))then
                 cpu_min(i)=impi
              endif
              if (   bound_key(impi-1).lt.bounding_max(i).and.&
                   & bound_key(impi  ).ge.bounding_max(i))then
                 cpu_max(i)=impi
              endif
           end do
        end do
        
        ncpu_read=0
        do i=1,ndom
           do j=cpu_min(i),cpu_max(i)
              if(.not. cpu_read(j))then
                 ncpu_read=ncpu_read+1
                 cpu_list(ncpu_read)=j
                 cpu_read(j)=.true.
              endif
           enddo
        enddo
     else
        ncpu_read=ncpu
        do j=1,ncpu
           cpu_list(j)=j
        end do
     end  if
     
     ! Create the index array for SN in processor myid
     do k=1,ncpu_read
        icpu=cpu_list(k)
        if(icpu==myid)then
           ii=ii+1
           iSN_myid(ii)=iSN
        endif
     enddo

  enddo
  
  ! Number of SN in processor myid
  nSN_myid=ii

end subroutine getSNonmyid

!*************************************************************************
SUBROUTINE output_SN_stats

! Print the average star formation rate in the simulation box during the
! coarse timestep
!-------------------------------------------------------------------------
  use amr_commons
  use SN_stats
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::  info
  real(dp):: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp):: E_SN_want_all, E_SN_got_all,pr_tot_all
  integer:: n_SN_want_all,n_SN_stoch_got_all,n_pr_all
!-------------------------------------------------------------------------
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(E_SN_want, E_SN_want_all,  1,        &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)          
  call MPI_ALLREDUCE(E_SN_got, E_SN_got_all,  1,        &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)          
  call MPI_ALLREDUCE(n_SN_want, n_SN_want_all,  1,        &
          MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)          
  call MPI_ALLREDUCE(n_SN_stoch_got, n_SN_stoch_got_all,  1,        &
          MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)          
  call MPI_ALLREDUCE(pr_tot, pr_tot_all,  1,        &
          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)          
  call MPI_ALLREDUCE(n_pr, n_pr_all,  1,        &
          MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)          
#endif

  if(myid .eq. 1) then
     write(*, 113) t * scale_t / (60d0*60d0*24d0*365d0*1d6)                &
                    ,E_SN_want_all * scale_d * scale_l**3 * scale_v**2 /1d51 &
                    ,E_SN_got_all * scale_d * scale_l**3 * scale_v**2 /1d51 &
                    ,pr_tot_all &
                    ,n_SN_want_all          &
                    ,n_SN_stoch_got_all   &
                    ,n_pr_all
  endif

113 format('SN_stats ',1pe12.5,1pe12.5,1pe12.5,1pe12.5,i9,i9,i9)

END SUBROUTINE output_SN_stats
