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
              if(idp(ipart).gt.0.and.tp(ipart).ne.0)then
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
              if(idp(ipart).gt.0.and.tp(ipart).ne.0)then
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

  ! Removed since this is done right after anyway (in move_particles)
  !! Check for illegal moves
  !error=.false.
  !do idim=1,ndim
  !   do j=1,np
  !      if(x(j,idim)<=0.5D0.or.x(j,idim)>=5.5D0)error=.true.
  !   end do
  !end do
  !if(error)then
  !   write(*,*)'problem in sn2'
  !   write(*,*)ilevel,ng,np
  !end if

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
  if(f_w==0)then
     do j=1,np
        birth_time=tp(ind_part(j))
        ! Make sure that we don't count feedback twice
        if(birth_time.lt.(current_time-t0).and.birth_time.ge.(current_time-t0-dteff(j)))then           
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
        endif
     end do
  endif

  ! Update hydro variables due to feedback

  ! For IR radiation trapping,
  ! we use the cell resolution to estimate the column density of gas
  delta_x=200*3d18
  if(metal)then
     tau_factor=kappa_IR*delta_x*scale_d/0.02
  else
     tau_factor=kappa_IR*delta_x*scale_d*z_ave
  endif
  rad_factor=ERAD/ESN
  do j=1,np
     ! Infrared photon trapping boost
     if(metal)then
        tauIR=tau_factor*max(unew(indp(j),imetal),smallr)
     else
        tauIR=tau_factor*max(unew(indp(j),1),smallr)
     endif
     if(unew(indp(j),1)*scale_nH > 10.)then
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
!################################################################
!################################################################
!################################################################
!################################################################
subroutine kinetic_feedback
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::nSN_tot_all
  integer,dimension(1:ncpu)::nSN_icpu_all
  real(dp),dimension(:),allocatable::mSN_all,sSN_all,ZSN_all
  real(dp),dimension(:,:),allocatable::xSN_all,vSN_all
#endif
  !----------------------------------------------------------------------
  ! This subroutine compute the kinetic feedback due to SNII and
  ! imolement this using exploding GMC particles. 
  ! This routine is called only at coarse time step.
  ! Yohan Dubois
  !----------------------------------------------------------------------
  ! local constants
  integer::ip,icpu,igrid,jgrid,npart1,npart2,ipart,jpart,next_part
  integer::nSN,nSN_loc,nSN_tot,info,iSN,ilevel,ivar
  integer,dimension(1:ncpu)::nSN_icpu
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,t0
  real(dp)::current_time
  real(dp)::scale,dx_min,vol_min,nISM,nCOM,d0,mstar
  integer::nx_loc
  integer,dimension(:),allocatable::ind_part,ind_grid
  logical,dimension(:),allocatable::ok_free
  integer ,dimension(:),allocatable::indSN
  real(dp),dimension(:),allocatable::mSN,sSN,ZSN,m_gas,vol_gas,ekBlast
  real(dp),dimension(:,:),allocatable::xSN,vSN,u_gas,dq

  if(.not. hydro)return
  if(ndim.ne.3)return

  if(verbose)write(*,*)'Entering make_sn'
  
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  ! Initial star particle mass
  mstar=n_star/(scale_nH*aexp**3)*vol_min

  ! Lifetime of Giant Molecular Clouds from Myr to code units
  ! Massive star lifetime from Myr to code units
  ! !MT!: replace 10 Myr by t_delay?
  if(use_proper_time)then
     t0=t_delay*1d6*(365.*24.*3600.)/(scale_t/aexp**2)
     current_time=texp
  else
     t0=t_delay*1d6*(365.*24.*3600.)/scale_t
     current_time=t
  endif

  !------------------------------------------------------
  ! Gather GMC particles eligible for disruption
  !------------------------------------------------------
  nSN_loc=0
  ! Loop over levels
  do icpu=1,ncpu
  ! Loop over cpus
     igrid=headl(icpu,levelmin)
     ! Loop over grids
     do jgrid=1,numbl(icpu,levelmin)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        ! Count old enough GMC particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).le.0.and. tp(ipart).lt.(current_time-t0))then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
            end do
        endif
        nSN_loc=nSN_loc+npart2   ! Add SNe to the total
        igrid=next(igrid)   ! Go to next grid
     end do
  end do
  ! End loop over levels
  nSN_icpu=0
  nSN_icpu(myid)=nSN_loc
#ifndef WITHOUTMPI
  ! Give an array of number of SN on each cpu available to all cpus
  call MPI_ALLREDUCE(nSN_icpu,nSN_icpu_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nSN_icpu=nSN_icpu_all
#endif

  nSN_tot=sum(nSN_icpu(1:ncpu))

  if (nSN_tot .eq. 0) return
  
  if(myid==1)then
     write(*,*)'-----------------------------------------------'
     write(*,*)'Number of GMC to explode=',nSN_tot
     write(*,*)'-----------------------------------------------'
  endif

  ! Allocate arrays for the position and the mass of the SN
  allocate(xSN(1:nSN_tot,1:3),vSN(1:nSN_tot,1:3))
  allocate(mSN(1:nSN_tot),sSN(1:nSN_tot),ZSN(1:nSN_tot))
  xSN=0.;vSN=0.;mSN=0.;sSN=0.;ZSN=0.
  ! Allocate arrays for particles index and parent grid
  if(nSN_loc>0)then
     allocate(ind_part(1:nSN_loc),ind_grid(1:nSN_loc),ok_free(1:nSN_loc))
  endif

  !------------------------------------------------------
  ! Store position and mass of the GMC into the SN array
  !------------------------------------------------------
  if(myid==1)then
     iSN=0
  else
     iSN=sum(nSN_icpu(1:myid-1))
  endif
  ! Loop over levels
  ip=0
  do icpu=1,ncpu
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
              if(idp(ipart).le.0.and. tp(ipart).lt.(current_time-t0))then
                 iSN=iSN+1
                 xSN(iSN,1)=xp(ipart,1)
                 xSN(iSN,2)=xp(ipart,2)
                 xSN(iSN,3)=xp(ipart,3)
                 vSN(iSN,1)=vp(ipart,1)
                 vSN(iSN,2)=vp(ipart,2)
                 vSN(iSN,3)=vp(ipart,3)
                 mSN(iSN)=mp(ipart)
                 sSN(iSN)=dble(-idp(ipart))*mstar
                 if(metal)ZSN(iSN)=zp(ipart)
                 ip=ip+1
                 ind_grid(ip)=igrid
                 ind_part(ip)=ipart
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        igrid=next(igrid)   ! Go to next grid
     end do
  end do 
  ! End loop over levels

  ! Remove GMC particle
  IF(nSN_loc>0)then
     ok_free=.true.
     call remove_list(ind_part,ind_grid,ok_free,nSN_loc)
     call add_free_cond(ind_part,ok_free,nSN_loc)
     deallocate(ind_part,ind_grid,ok_free)
  endif

#ifndef WITHOUTMPI
  allocate(xSN_all(1:nSN_tot,1:3),vSN_all(1:nSN_tot,1:3),mSN_all(1:nSN_tot),sSN_all(1:nSN_tot),ZSN_all(1:nSN_tot))
  call MPI_ALLREDUCE(xSN,xSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vSN,vSN_all,nSN_tot*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mSN,mSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(sSN,sSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ZSN,ZSN_all,nSN_tot  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  xSN=xSN_all
  vSN=vSN_all
  mSN=mSN_all
  sSN=sSN_all
  ZSN=ZSN_all
  deallocate(xSN_all,vSN_all,mSN_all,sSN_all,ZSN_all)
#endif

  nSN=nSN_tot
  allocate(m_gas(1:nSN),u_gas(1:nSN,1:3),vol_gas(1:nSN),dq(1:nSN,1:3),ekBlast(1:nSN))
  allocate(indSN(1:nSN))

  ! Compute the grid discretization effects
  call average_SN(xSN,vol_gas,dq,ekBlast,indSN,nSN)

  ! Modify hydro quantities to account for a Sedov blast wave
  call Sedov_blast(xSN,vSN,mSN,sSN,ZSN,indSN,vol_gas,dq,ekBlast,nSN)

  deallocate(xSN,vSN,mSN,sSN,ZSN,indSN,m_gas,u_gas,vol_gas,dq,ekBlast)

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
subroutine average_SN(xSN,vol_gas,dq,ekBlast,ind_blast,nSN)
  use pm_commons
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine average the hydro quantities inside the SN bubble
  !------------------------------------------------------------------------
  integer::ilevel,ncache,nSN,j,iSN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,info
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dr_SN,d,u,v,w,ek,u2,v2,w2,dr_cell
  real(dp)::scale,dx,dxx,dyy,dzz,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ,dimension(1:nSN)::ind_blast
  real(dp),dimension(1:nSN)::mSN,m_gas,vol_gas,ekBlast
  real(dp),dimension(1:nSN,1:3)::xSN,vSN,u_gas,dq,u2Blast
#ifndef WITHOUTMPI
  real(dp),dimension(1:nSN)::m_gas_all,vol_gas_all,ekBlast_all
  real(dp),dimension(1:nSN,1:3)::u_gas_all,dq_all,u2Blast_all
#endif
  logical ,dimension(1:nvector),save::ok

  if(nSN==0)return
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
  rmax=MAX(2.0d0*dx_min*scale_l/aexp,rbubble*3.08d18)
  rmax=rmax/scale_l
  rmax2=rmax*rmax

  ! Initialize the averaged variables
  vol_gas=0.0;dq=0.0;u2Blast=0.0;ekBlast=0.0;ind_blast=-1

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
                       ind_blast(iSN)=ind_cell(i)
                       ekBlast  (iSN)=vol_loc
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

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(vol_gas,vol_gas_all,nSN  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dq     ,dq_all     ,nSN*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(u2Blast,u2Blast_all,nSN*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ekBlast,ekBlast_all,nSN  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  vol_gas=vol_gas_all
  dq     =dq_all
  u2Blast=u2Blast_all
  ekBlast=ekBlast_all
#endif
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
        ekBlast(iSN)=max(0.5d0*(u2+v2+w2),0.0d0)
     endif
  end do

  if(verbose)write(*,*)'Exiting average_SN'

end subroutine average_SN
!################################################################
!################################################################
!################################################################
!################################################################
subroutine Sedov_blast(xSN,vSN,mSN,sSN,ZSN,indSN,vol_gas,dq,ekBlast,nSN)
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
  real(dp)::x,y,z,dx,dxx,dyy,dzz,dr_SN,d,u,v,w,ek,u_r,ESN
  real(dp)::scale,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nSN)::mSN,sSN,ZSN,m_gas,p_gas,d_gas,d_metal,vol_gas,uSedov,ekBlast
  real(dp),dimension(1:nSN,1:3)::xSN,vSN,u_gas,dq
  integer ,dimension(1:nSN)::indSN
  logical ,dimension(1:nvector),save::ok

  if(nSN==0)return
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
  rmax=MAX(2.0d0*dx_min*scale_l/aexp,rbubble*3.08d18)
  rmax=rmax/scale_l
  rmax2=rmax*rmax
  
  ! Supernova specific energy from cgs to code units
  ESN=(1d51/(10d0*2d33))/scale_v**2

  do iSN=1,nSN
     if(vol_gas(iSN)>0d0)then
        d_gas(iSN)=mSN(iSN)/vol_gas(iSN)
        if(metal)d_metal(iSN)=ZSN(iSN)*mSN(iSN)/vol_gas(iSN)
        if(ekBlast(iSN)==0d0)then
           p_gas(iSN)=eta_sn*sSN(iSN)*ESN/vol_gas(iSN)
           uSedov(iSN)=0d0
        else
           p_gas(iSN)=(1d0-f_ek)*eta_sn*sSN(iSN)*ESN/vol_gas(iSN)
           uSedov(iSN)=sqrt(f_ek*eta_sn*sSN(iSN)*ESN/mSN(iSN)/ekBlast(iSN))
        endif
     else
        d_gas(iSN)=mSN(iSN)/ekBlast(iSN)
        p_gas(iSN)=eta_sn*sSN(iSN)*ESN/ekBlast(iSN)
        if(metal)d_metal(iSN)=ZSN(iSN)*mSN(iSN)/ekBlast(iSN)
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
                       ! Compute the mass density in the cell
                       uold(ind_cell(i),1)=uold(ind_cell(i),1)+d_gas(iSN)
                       ! Compute the metal density in the cell
                       if(metal)uold(ind_cell(i),imetal)=uold(ind_cell(i),imetal)+d_metal(iSN)
                       ! Velocity at a given dr_SN linearly interpolated between zero and uSedov
                       u=uSedov(iSN)*(dxx/rmax-dq(iSN,1))+vSN(iSN,1)
                       v=uSedov(iSN)*(dyy/rmax-dq(iSN,2))+vSN(iSN,2)
                       w=uSedov(iSN)*(dzz/rmax-dq(iSN,3))+vSN(iSN,3)
                       ! Add each momentum component of the blast wave to the gas
                       uold(ind_cell(i),2)=uold(ind_cell(i),2)+d_gas(iSN)*u
                       uold(ind_cell(i),3)=uold(ind_cell(i),3)+d_gas(iSN)*v
                       uold(ind_cell(i),4)=uold(ind_cell(i),4)+d_gas(iSN)*w
                       ! Finally update the total energy of the gas
                       uold(ind_cell(i),5)=uold(ind_cell(i),5)+0.5*d_gas(iSN)*(u*u+v*v+w*w)+p_gas(iSN)
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
        u=vSN(iSN,1)
        v=vSN(iSN,2)
        w=vSN(iSN,3)
        if(indSN(iSN)>0)then
           uold(indSN(iSN),1)=uold(indSN(iSN),1)+d_gas(iSN)
           uold(indSN(iSN),2)=uold(indSN(iSN),2)+d_gas(iSN)*u
           uold(indSN(iSN),3)=uold(indSN(iSN),3)+d_gas(iSN)*v
           uold(indSN(iSN),4)=uold(indSN(iSN),4)+d_gas(iSN)*w
           uold(indSN(iSN),5)=uold(indSN(iSN),5)+d_gas(iSN)*0.5*(u*u+v*v+w*w)+p_gas(iSN)
           if(metal)uold(indSN(iSN),imetal)=uold(indSN(iSN),imetal)+d_metal(iSN)
        endif
     endif
  end do

  if(verbose)write(*,*)'Exiting Sedov_blast'

end subroutine Sedov_blast
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine kinetic_feedback_cell
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer,dimension(1:ncpu)::nSNcomm_icpu_mpi
  real(dp),dimension(:),allocatable::mSN_mpi,ZSN_mpi
  real(dp),dimension(:,:),allocatable::xSN_mpi,vSN_mpi
  integer, dimension(:),allocatable::lSN_mpi
#endif
  !----------------------------------------------------------------------
  ! Description: This subroutine gathers SN events on a cell-by-cell basis
  ! Here 'local' variable means things that do not require MPI comm
  ! Taysun Kimm
  !----------------------------------------------------------------------
  ! local constants
  integer::igrid,npart1,ipart,jpart,next_part
  integer::nSNp,nSNp_tot,nSNc,nSNc_tot,nSNcomm,nSNcomm_tot,nSN_glo,nSN_loc
  integer::ilevel,ind,ivar,nx_loc,iSN,nSN
  integer,dimension(1:ncpu)::nSNcomm_icpu
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,t0
  real(dp)::scale,dx,dx_loc,skip_loc(3),current_time
  integer ,dimension(:),allocatable::iSN_myid
  real(dp),dimension(:),allocatable::mSN_loc,ZSN_loc
  real(dp),dimension(:,:),allocatable::xSN_loc,vSN_loc
  real(dp),dimension(:),allocatable::mSN,ZSN
  real(dp),dimension(:,:),allocatable::xSN,vSN
  integer ,dimension(:),allocatable::itemp,lSN_loc,lSN_glo,lSN
  real(dp),dimension(:),allocatable::mSN_glo,ZSN_glo
  real(dp),dimension(:,:),allocatable::xSN_glo,vSN_glo
  real(dp),dimension(1:twotondim,1:3)::xc,p_ej
  real(dp),dimension(1:twotondim)::m_ej,Z_ej
  real(dp),dimension(1:ndim)::x0,x
  real(dp)::ttsta,ttend,mejecta
  integer::ind_grid,ind_cell,ngrid,ind_son,idim,iSN_loc,iSN_glo
  integer::isort,iskip,i,ix,iy,iz,info,nok_sn,ncpu_read
  logical::ok

  if(.not. hydro)return
  if(ndim.ne.3)return

#ifndef WITHOUTMPI
  ttsta = MPI_WTIME(info)
#endif
  
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Lifetime of Giant Molecular Clouds from Myr to code units
  ! Massive star lifetime from Myr to code units
  if(use_proper_time)then
     t0=t_delay*1d6*(365.*24.*3600.)/(scale_t/aexp**2)
     current_time=texp
  else
     t0=t_delay*1d6*(365.*24.*3600.)/scale_t
     current_time=t
  endif


  ! update the particle lists on ilevel > levelmin
  do ilevel=levelmin,nlevelmax
     call kill_tree_fine(ilevel)
     call virtual_tree_fine(ilevel)
  end do

  ! gather particle from the grid
  do ilevel=nlevelmax-1,1,-1
     call merge_tree_fine(ilevel)
  end do

  ! Scatter particle to the grid
  do ilevel=1,nlevelmax
     call make_tree_fine(ilevel)
     call kill_tree_fine(ilevel)
     call virtual_tree_fine(ilevel)
  end do

  ! gather particle
  do ilevel=nlevelmax,levelmin,-1
     call merge_tree_fine(ilevel)
  end do

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)


  !------------------------------------------------------
  ! Count cells eligible for a SN event
  !------------------------------------------------------
  nSNc     = 0  ! the total number of cells that will launch SN in myid
  nSNp     = 0  ! the total number of SN particles in myid
  nSNcomm  = 0  ! the total number of cells that need communications
  nSNc_tot = 0  ! total(nSN) for the entire cpus
  nSNp_tot = 0  ! total(nSNp) for the entire cpus
  do ilevel=levelmin,nlevelmax
     dx=0.5D0**ilevel
     dx_loc=dx*scale

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
     ngrid=active(ilevel)%ngrid
     do igrid=1,ngrid
        ind_grid=active(ilevel)%igrid(igrid)
        npart1=numbp(ind_grid) ! number of particles in the grid
        ok=.true.
        if(npart1>0)then
           do idim=1,ndim
              x0(idim) = xg(ind_grid,idim) - dx
           end do

           do ind=1,twotondim  ! Loop over cells
              iskip=ncoarse+(ind-1)*ngridmax
              ind_cell=iskip+ind_grid
              ok=(ok.and.(son(ind_cell)==0))   ! Flag leaf cells
           end do

           if(ok)then
              m_ej=0d0 
              ipart=headp(ind_grid)
              do jpart=1,npart1
                 next_part=nextp(ipart)
                 if(idp(ipart).le.0.and.tp(ipart).lt.(current_time-t0))then
                     nSNp=nSNp+1
                     ! find the cell index to get the position of it
                     ind_son=1  ! double-checked
                     do idim=1,ndim
                        i = int((xp(ipart,idim)/scale+skip_loc(idim)-x0(idim))/dx)
                        ind_son=ind_son+i*2**(idim-1)
                     end do
                     m_ej(ind_son)=1.
                     !ind_cell = ncoarse + (ind_son-1)*ngridmax + ind_grid
                 endif
                 ipart=next_part
              enddo

              nok_sn = count(m_ej>0d0)
              nSNc = nSNc + nok_sn ! add the number of cells that contain SN(e)

              if (nok_sn>0) then ! check how many SNe need MPI comm
                 do ind=1,twotondim
                    if (m_ej(ind)>0.5) then
                       x(1) = (xg(ind_grid,1)+xc(ind,1)-skip_loc(1))*scale 
                       x(2) = (xg(ind_grid,2)+xc(ind,2)-skip_loc(2))*scale 
                       x(3) = (xg(ind_grid,3)+xc(ind,3)-skip_loc(3))*scale 
                       call checkSNboundary(x,ncpu_read,ilevel)
                       if(ncpu_read>1) nSNcomm = nSNcomm + 1
                    endif
                 enddo
              endif


           endif
        endif
     enddo
  enddo
  
  nSNp_tot=nSNp
  nSNc_tot=nSNc
  nSNcomm_tot=0 
  nSNcomm_icpu=0  ! This array is for cells that need MPI communications
  nSNcomm_icpu(myid)=nSNcomm
#ifndef WITHOUTMPI
  ! Give an array of number of SN on each cpu available to all cpus
  nSNcomm_icpu_mpi=0
  nSNp_tot=0
  call MPI_ALLREDUCE(nSNc    ,nSNc_tot    ,1,   MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(nSNp    ,nSNp_tot    ,1,   MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(nSNcomm_icpu,nSNcomm_icpu_mpi,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nSNcomm_icpu=nSNcomm_icpu_mpi
#endif
  nSNcomm_tot=sum(nSNcomm_icpu(1:ncpu))  ! the total number of cells that contain SNe for all cpus

  if(myid==1.and.nSNc_tot>0)then
     write(*,*)'-----------------------------------------------'
     write(*,*)'Number of SN/SNc/SNcomm to explode=',nSNp_tot,nSNc_tot,nSNcomm_tot
     write(*,*)'-----------------------------------------------'
  endif

  if (nSNc_tot .eq. 0) return


  ! local SN cells except the ones that will go to "global". This wouldn't require MPI
  nSN_loc = nSNc - nSNcomm
  nSN_glo = nSNcomm_tot

  ! Allocate arrays for the position and the mass of the SN
  allocate(xSN_loc(1:nSN_loc,1:3),vSN_loc(1:nSN_loc,1:3))
  allocate(mSN_loc(1:nSN_loc),ZSN_loc(1:nSN_loc),lSN_loc(1:nSN_loc))
  xSN_loc=0.;vSN_loc=0.;mSN_loc=0.;ZSN_loc=0.; lSN_loc=0

  allocate(xSN_glo(1:nSN_glo,1:3),vSN_glo(1:nSN_glo,1:3))
  allocate(mSN_glo(1:nSN_glo),ZSN_glo(1:nSN_glo),lSN_glo(1:nSN_glo))
  xSN_glo=0.;vSN_glo=0.;mSN_glo=0.;ZSN_glo=0.; lSN_glo=0

  iSN_loc=0
  if(myid==1) then
     iSN_glo=0
  else
     iSN_glo=sum(nSNcomm_icpu(1:myid-1))
  endif


  !---------------------------------------------------------------------
  ! Gather ejecta m, p, Z, etc., for local/global SN events, separately
  ! here 'global' means that it requires cells belonging to other cpus 
  !---------------------------------------------------------------------
  do ilevel=levelmin,nlevelmax
     dx=0.5D0**ilevel
     dx_loc=dx*scale

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
     ngrid=active(ilevel)%ngrid
     do igrid=1,ngrid
        ind_grid=active(ilevel)%igrid(igrid)
        npart1=numbp(ind_grid) ! number of particles in the grid
        ok=.true.
        if(npart1>0)then
           do idim=1,ndim
              x0(idim) = xg(ind_grid,idim) - dx
           end do

           do ind=1,twotondim  ! Loop over cells
              iskip=ncoarse+(ind-1)*ngridmax
              ind_cell=iskip+ind_grid
              ok=(ok.and.(son(ind_cell)==0))   ! Flag leaf cells
           end do

           if(ok)then
              m_ej=0d0; p_ej=0d0; Z_ej=0d0
              ipart=headp(ind_grid)
              do jpart=1,npart1
                 next_part=nextp(ipart)
                 if(idp(ipart).le.0.and.tp(ipart).lt.(current_time-t0))then
                     nSNp=nSNp+1
                     ! find the cell index to get the position of it
                     ind_son=1  ! double-checked
                     do idim=1,ndim
                        i = int((xp(ipart,idim)/scale+skip_loc(idim)-x0(idim))/dx)
                        ind_son=ind_son+i*2**(idim-1)
                     end do

                     mejecta = mp(ipart)*eta_sn
                     m_ej(ind_son) = m_ej(ind_son) + mejecta
                     p_ej(ind_son,1) = p_ej(ind_son,1) + mejecta*vp(ipart,1)
                     p_ej(ind_son,2) = p_ej(ind_son,2) + mejecta*vp(ipart,2)
                     p_ej(ind_son,3) = p_ej(ind_son,3) + mejecta*vp(ipart,3)
                     if(metal) Z_ej(ind_son) = Z_ej(ind_son) + mejecta*zp(ipart)
                     ! Remove the mass ejected by the SN
                     mp(ipart)  = mp(ipart) - mejecta
                     idp(ipart) = -idp(ipart)
                     !ind_cell = ncoarse + (ind_son-1)*ngridmax + ind_grid
                 endif
                 ipart=next_part
              enddo

              nok_sn = count(m_ej>0d0)
              nSNc = nSNc + nok_sn ! add the number of cells that contain SN(e)

              if (nok_sn>0) then ! check how many SNe need MPI comm
                 do ind=1,twotondim
                    if (m_ej(ind)>0d0) then
                       x(1) = (xg(ind_grid,1)+xc(ind,1)-skip_loc(1))*scale 
                       x(2) = (xg(ind_grid,2)+xc(ind,2)-skip_loc(2))*scale 
                       x(3) = (xg(ind_grid,3)+xc(ind,3)-skip_loc(3))*scale 
                       call checkSNboundary(x,ncpu_read,ilevel)
                       if(ncpu_read>1)then  ! SNe that requires cells in other cpus
                          iSN_glo = iSN_glo + 1
                          xSN_glo(iSN_glo,1) = x(1)
                          xSN_glo(iSN_glo,2) = x(2)
                          xSN_glo(iSN_glo,3) = x(3)
                          mSN_glo(iSN_glo  ) = m_ej(ind)
                          vSN_glo(iSN_glo,1) = p_ej(ind,1)/m_ej(ind)
                          vSN_glo(iSN_glo,2) = p_ej(ind,2)/m_ej(ind)
                          vSN_glo(iSN_glo,3) = p_ej(ind,3)/m_ej(ind)
                          if(metal)ZSN_glo(iSN_glo)=Z_ej(ind)/m_ej(ind)
                          lSN_glo(iSN_glo  ) = ilevel
                       else                 ! SNe that are happy with cells in myid
                          iSN_loc = iSN_loc + 1
                          xSN_loc(iSN_loc,1) = x(1)
                          xSN_loc(iSN_loc,2) = x(2)
                          xSN_loc(iSN_loc,3) = x(3)
                          mSN_loc(iSN_loc  ) = m_ej(ind)
                          vSN_loc(iSN_loc,1) = p_ej(ind,1)/m_ej(ind)
                          vSN_loc(iSN_loc,2) = p_ej(ind,2)/m_ej(ind)
                          vSN_loc(iSN_loc,3) = p_ej(ind,3)/m_ej(ind)
                          if(metal)ZSN_loc(iSN_loc)=Z_ej(ind)/m_ej(ind)
                          lSN_loc(iSN_loc  ) = ilevel
                       endif
                    endif
                 enddo
              endif


           endif
        endif
     enddo
  enddo
  


#ifndef WITHOUTMPI
  !---------------------------------------------------------------------
  !  Part-1 release: global SNe
  !---------------------------------------------------------------------
  if (nSN_glo>0) then

     allocate(xSN_mpi(1:nSN_glo,3),vSN_mpi(1:nSN_glo,3))
     allocate(mSN_mpi(1:nSN_glo),ZSN_mpi(1:nSN_glo),lSN_mpi(1:nSN_glo))
     xSN_mpi=0d0;vSN_mpi=0d0;mSN_mpi=0d0;ZSN_mpi=0d0;lSN_mpi=0
     call MPI_ALLREDUCE(xSN_glo,xSN_mpi,nSN_glo*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(vSN_glo,vSN_mpi,nSN_glo*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(mSN_glo,mSN_mpi,nSN_glo  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(lSN_glo,lSN_mpi,nSN_glo  ,MPI_INTEGER         ,MPI_SUM,MPI_COMM_WORLD,info)
     if(metal)call MPI_ALLREDUCE(ZSN_glo,ZSN_mpi,nSN_glo,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     xSN_glo=xSN_mpi
     vSN_glo=vSN_mpi
     mSN_glo=mSN_mpi
     lSN_glo=lSN_mpi
     if(metal)ZSN_glo=ZSN_mpi
     deallocate(xSN_mpi,vSN_mpi,mSN_mpi,ZSN_mpi,lSN_mpi)

     allocate(itemp(1:nSN_glo))
     itemp=0
     call getSNonmyid2(itemp,nSN,xSN_glo,nSN_glo,lSN_glo)

     allocate(xSN(1:nSN,1:3),vSN(1:nSN,1:3),mSN(1:nSN),ZSN(1:nSN),iSN_myid(1:nSN),lSN(1:nSN))
     xSN=0d0; vSN=0d0; mSN=0d0; ZSN=0d0; iSN_myid=0; lSN=0

     do iSN=1,nSN
        isort = itemp(iSN)
        iSN_myid(iSN) = isort
        xSN(iSN,1) = xSN_glo(isort,1)
        xSN(iSN,2) = xSN_glo(isort,2)
        xSN(iSN,3) = xSN_glo(isort,3)
        vSN(iSN,1) = vSN_glo(isort,1)
        vSN(iSN,2) = vSN_glo(isort,2)
        vSN(iSN,3) = vSN_glo(isort,3)
        mSN(iSN  ) = mSN_glo(isort)
        ZSN(iSN  ) = ZSN_glo(isort)
        lSN(iSN  ) = lSN_glo(isort)
     end do

     ! Inject momentum
     call inject_momentum_SN(xSN,vSN,mSN,ZSN,lSN,nSN,iSN_myid,nSN_glo,.true.)

     deallocate(xSN,vSN,mSN,ZSN,iSN_myid,itemp,lSN)

  endif
#endif
  deallocate(xSN_glo,vSN_glo,mSN_glo,ZSN_glo,lSN_glo)

  !---------------------------------------------------------------------
  !  Part-2 release: local SNe
  !---------------------------------------------------------------------
  allocate(itemp(1:nSN_loc))
  itemp=0
  call getSNonmyid2(itemp,nSN,xSN_loc,nSN_loc,lSN_loc)

  allocate(xSN(1:nSN,1:3),vSN(1:nSN,1:3),mSN(1:nSN),ZSN(1:nSN),iSN_myid(1:nSN),lSN(1:nSN))
  xSN=0d0; vSN=0d0; mSN=0d0; ZSN=0d0; iSN_myid=0; lSN=0

  do iSN=1,nSN
     isort = itemp(iSN)
     iSN_myid(iSN) = isort
     xSN(iSN,1) = xSN_loc(isort,1)
     xSN(iSN,2) = xSN_loc(isort,2)
     xSN(iSN,3) = xSN_loc(isort,3)
     vSN(iSN,1) = vSN_loc(isort,1)
     vSN(iSN,2) = vSN_loc(isort,2)
     vSN(iSN,3) = vSN_loc(isort,3)
     mSN(iSN  ) = mSN_loc(isort)
     ZSN(iSN  ) = ZSN_loc(isort)
     lSN(iSN  ) = lSN_loc(isort)
  end do

  ! Inject momentum
  call inject_momentum_SN(xSN,vSN,mSN,ZSN,lSN,nSN,iSN_myid,nSN_loc,.false.)

  deallocate(xSN,vSN,mSN,ZSN,iSN_myid,itemp,lSN)
  deallocate(xSN_loc,vSN_loc,mSN_loc,ZSN_loc)

  ! Update hydro quantities for split cells
  do ilevel=nlevelmax,levelmin,-1
     call upload_fine(ilevel)
     do ivar=1,nvar
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     enddo
  enddo


#ifndef WITHOUTMPI
  ttend = MPI_WTIME(info)
  if(myid.eq.1)then
     write(*,*) 'Time elapsed in kinetic_feedback_cell [s]', sngl(ttend-ttsta)
  endif
#endif

end subroutine kinetic_feedback_cell
!################################################################
!################################################################
!################################################################
!################################################################
subroutine inject_momentum_SN(xSN,vSN,mSN,ZSN,lSN,nSN,iSN_myid,nSN_glo,global_search)
  use amr_commons
  use hydro_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  real(dp),dimension(1:nSN_glo,1:3)::vloadSN_mpi
  real(dp),dimension(1:nSN_glo,1:3)::dq_mpi
  real(dp),dimension(1:nSN_glo)::mloadSN_mpi,ZloadSN_mpi
#endif
  integer::ilevel,iSN,nSN,nSN_glo,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,ncache,ind_SN,info
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dx,dxx,dyy,dzz,drr,f_adjacent
  real(dp)::ekk,eint,eint2,u,v,w,d,vol_nbor,dload,pload(3),dloadloc
  real(dp)::scale,dx_min,dx_loc,vol_loc,rloose2,dx_sn,dd_sn
  real(dp)::scale_nh,scale_t2,scale_l,scale_d,scale_t,scale_v,scale_msun
  real(dp)::pxmi,pxma,pymi,pyma,pzmi,pzma,rmaxloc,dr_cell
  real(dp)::mload,vload,Zload,M_SNII,ESN,f_w_cell,navg
  real(dp)::f_esn,f_w_crit,num_sn,Zdepen,f_leftover,f_load
  real(dp),dimension(1:3)::skip_loc,x0
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nSN)::mSN,ZSN
  real(dp),dimension(1:nSN,1:3)::xSN,vSN
  integer, dimension(1:nSN)::lSN,ind_blast
  real(dp),dimension(1:nSN_glo)::mloadSN,ZloadSN
  real(dp),dimension(1:nSN_glo,1:3)::vloadSN,dq
  integer, dimension(1:nSN_glo)::iSN_myid
  logical::global_search
  logical,dimension(1:nvector),save::ok
!--------------------------------------------------------------
  integer::ista,iend,jsta,jend,ksta,kend,n48,icello
  real(dp),dimension(1:48,1:nSN)::rho_nbor,Z_nbor ! rho for neigboring cells
  integer ,dimension(1:48,1:nSN)::icell_nbor ! icell to speed up
  integer ,dimension(1:48,1:nSN)::lv_nbor
  integer ,dimension(1:48)::idx48
  real(dp),dimension(1:48,1:3)::v48
!--------------------------------------------------------------
! TKNOTE
!--------------------------------------------------------------
! There are recurrent numbers like '48' '56'
! The structure one can suppose is 4*4*4 of which the central 2*2*2 is a SN cell.
! The total number of cells in this structure is 64, but we decide to neglect 
! 8 vertices, making the total number to be 56. On top of that, 
! what we inject momentum is only for neibouring cells, making the relevant cells 
! to be 48. Note that we deposit 8/56 of the SN mass + mass in the SN cell. This is desirable 
! for numerical stability. If there are many stars with different ages in a cell,
! it will redistribute mass very efficiently, in which case the cell is gonna run 
! out of gas. 

  f_leftover=1d0/56d0  ! M_leftover in the SN cell = f_leftover*(m_cell(iSN) + m_ejecta)
  f_load    =1d0 - f_leftover

  ! mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5d0**nlevelmax
  f_adjacent = 1.7 ! the longest distance=1.658 if we neglect 8 edge cells  

  ! conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nh,scale_t2)
  scale_msun=(scale_d*scale_l**3d0)/2d33 

  ! In the original code, M_SN=10Msun is used. If you want to stick to this,
  ! eff_sfbk=0.657510, eta_sn=0.31675331 (Chabrier)
  ! This does not mean that the actual efficiency is 66%....
  ! Here I change M_SNII to make things more intuitive.
  ESN=1d51
  M_SNII = 15.208892     ! Mean SN progenitor mass for Chabrier
  f_esn  = 0.676 ! Blondin+(98) at t=trad (applies only for high f_w_cell)

  !              Salpeter / Kroupa / Chabirer
  ! M_SNII       18.72825  19.13473  15.20889
  ! eta          0.138997  0.213007  0.316753
  ! eta/M_SNII   0.007422  0.011132  0.020827   (= N_SN/Msun)    


  ! To reduce time spent to search for cells
  pxmi=minval(xSN(:,1))
  pxma=maxval(xSN(:,1))
  pymi=minval(xSN(:,2))
  pyma=maxval(xSN(:,2))
  pzmi=minval(xSN(:,3))
  pzma=maxval(xSN(:,3))

  ind_blast = -1

  mloadSN =0d0; ZloadSN=0d0; vloadSN=0d0; rho_nbor=0d0; Z_nbor=0d0; icell_nbor=-1


  !---------------------------------------------------------------------
  ! Gather rho,v,Z,icell for SN + adjacent cell
  !---------------------------------------------------------------------
  ! loop over levels
  do ilevel=levelmin,nlevelmax
     ! computing local volume (important for averaging hydro quantities) 
     dx=0.5d0**ilevel
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     rmaxloc = 2*dx_loc 
     rloose2 = (f_adjacent*dx_loc*2)**2d0

     ! cells center position relative to grid center position
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5d0)*dx
        xc(ind,2)=(dble(iy)-0.5d0)*dx
        xc(ind,3)=(dble(iz)-0.5d0)*dx
     end do

     ! loop over grids
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=min(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           ! To speed up
           do i=1,ngrid
              if(ok(i))then
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 if((x.lt.pxmi-rmaxloc).or.(x.gt.pxma+rmaxloc).or.&
                   &(y.lt.pymi-rmaxloc).or.(y.gt.pyma+rmaxloc).or.&
                   &(z.lt.pzmi-rmaxloc).or.(z.gt.pzma+rmaxloc)) then
                    ok(i)=.false.
                 endif
              endif
           enddo
 
           do i=1,ngrid
              if(ok(i))then
                 ! get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
                 do iSN=1,nSN
                    ind_SN = iSN_myid(iSN)
                    dxx = x - xSN(iSN,1)
                    dyy = y - xSN(iSN,2)
                    dzz = z - xSN(iSN,3)
                    drr = dxx*dxx + dyy*dyy + dzz*dzz
                    if(drr<=rloose2*2)then
                       dx_sn = 0.5d0**lSN(iSN)*scale    ! cell size at which iSN exist
                       dd_sn = (dx_loc+dx_sn)/2d0 + dx_sn/1d5 ! maximum separation at this level
                       dr_cell = max(abs(dxx),abs(dyy),abs(dzz))
                       if(dr_cell<dd_sn*2) then !if adjacent
                          x0(1) = xSN(iSN,1) - dx_sn ! left-bottom pos of the SN region (not cell)
                          x0(2) = xSN(iSN,2) - dx_sn
                          x0(3) = xSN(iSN,3) - dx_sn
                          
                          ista  = idnint(((x-dx_loc/2d0) - x0(1))/dx_sn*2d0)
                          iend  = idnint(((x+dx_loc/2d0) - x0(1))/dx_sn*2d0) -1
                          jsta  = idnint(((y-dx_loc/2d0) - x0(2))/dx_sn*2d0)
                          jend  = idnint(((y+dx_loc/2d0) - x0(2))/dx_sn*2d0) -1
                          ksta  = idnint(((z-dx_loc/2d0) - x0(3))/dx_sn*2d0)
                          kend  = idnint(((z+dx_loc/2d0) - x0(3))/dx_sn*2d0) -1

                          ! calculate special indices designed for momentum
                          ! injection. if this cell is twice larger than 
                          ! the SN  cell, the number of indices for this cell is
                          ! gonna be 4.
                          call get_idx48(ista,iend,jsta,jend,ksta,kend,idx48,n48)

                          ! remember surrounding densities
                          rho_nbor(idx48(1:n48),iSN)   =uold(ind_cell(i),1) 
                          if(metal)then
                             Z_nbor(idx48(1:n48),iSN)  =uold(ind_cell(i),imetal)/uold(ind_cell(i),1)
                          endif
                          icell_nbor(idx48(1:n48),iSN) =ind_cell(i)
                          lv_nbor(idx48(1:n48),iSN)    =ilevel
                          if(dr_cell<dx_loc/2d0)then  ! if the SN cell, not adjacent cells 
                             ind_blast(iSN)=ind_cell(i)
                             d = uold(ind_blast(iSN),1)
                             mload = f_load*d*vol_loc !min(f_w*mSN(iSN),0.90d0*d*vol_loc)
                             mloadSN(ind_SN) = mSN(iSN)*f_load + mload
                             if(metal)then
                                Zload = uold(ind_blast(iSN),imetal)/d
                                ZloadSN(ind_SN)=(mload*Zload + (yield+(1d0-yield)*ZSN(iSN))*mSN(iSN)*f_load)/mloadSN(ind_SN)
                             endif
                             u = uold(ind_blast(iSN),1)/d  !PROBLEM!!
                             v = uold(ind_blast(iSN),2)/d  !PROBLEM!!
                             w = uold(ind_blast(iSN),3)/d  !PROBLEM!!
                             vloadSN(ind_SN,1)=(f_load*mSN(iSN)*vSN(iSN,1)+mload*u)/mloadSN(ind_SN)
                             vloadSN(ind_SN,2)=(f_load*mSN(iSN)*vSN(iSN,2)+mload*v)/mloadSN(ind_SN)
                             vloadSN(ind_SN,3)=(f_load*mSN(iSN)*vSN(iSN,3)+mload*w)/mloadSN(ind_SN)
                          endif
                       endif
                        
                    end if 
                 end do
              end if
           end do

        end do ! loop over cells
     end do ! loop over grids
  end do ! loop over levels


#ifndef WITHOUTMPI
  if(global_search)then
     mloadSN_mpi=0d0; vloadSN_mpi=0d0; ZloadSN_mpi=0d0
     call MPI_ALLREDUCE(mloadSN,mloadSN_mpi,nSN_glo  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(vloadSN,vloadSN_mpi,nSN_glo*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     if(metal)call MPI_ALLREDUCE(ZloadSN,ZloadSN_mpi,nSN_glo,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     mloadSN = mloadSN_mpi
     vloadSN = vloadSN_mpi
     ZloadSN = ZloadSN_mpi
  endif
#endif


  !---------------------------------------------------------------------
  ! Calculate the local mass loading and enforce momentum conservation
  !---------------------------------------------------------------------
  dq=0d0 ; v48=0d0
  call get_v48(v48)

  do iSN=1,nSN
     ind_SN = iSN_myid(iSN)
     dx=0.5d0**lSN(iSN)
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim   
     mload=mloadSN(ind_SN)
     if(metal)Zload=ZloadSN(ind_SN)


     do i=1,48
        if(icell_nbor(i,iSN)>0)then
           ! idx48 is measured for sub volume of vol_loc/8
           f_w_cell = (mload/48d0 + rho_nbor(i,iSN)*vol_loc/8d0) / (mSN(iSN)/56d0) - 1d0
           vol_nbor = (0.5d0**lv_nbor(i,iSN)*scale)**ndim
           dloadloc = mload/48d0/vol_loc  ! double-checked
           navg     = (dloadloc+rho_nbor(i,iSN))/2d0*scale_nH
           num_sn   = (mSN(iSN)/(M_SNII/scale_msun))
           if(metal)then
              Zdepen = (dloadloc*Zload+rho_nbor(i,iSN)*Z_nbor(i,iSN))/(dloadloc+rho_nbor(i,iSN))
              Zdepen = (max(0.01,Zdepen/0.02))**(-0.28) ! From Thornton et al. (1999)
           else
              Zdepen = 1d0
           endif 
           !From Blondin et al. (1998)
           f_w_crit = max(0d0,9d2/(f_esn*M_SNII)*num_sn**(-2d0/17d0)*navg**(-4d0/17d0)*Zdepen-1d0)
           if(f_w_cell.ge.f_w_crit)then !radiative phase (momentum-conserving phase)
              vload    = dsqrt(2d0*f_esn*ESN*(1+f_w_crit)/(M_SNII*2d33))/scale_v/(1+f_w_cell)
           else
              vload    = dsqrt(2d0      *ESN/(1+f_w_cell)/(M_SNII*2d33))/scale_v
           endif

           ! dq: net residual momentum for ind_SN
           dq(ind_SN,1) = dq(ind_SN,1) + ((1+f_w_cell)*mSN(iSN)/56d0)*vload*v48(i,1) 
           dq(ind_SN,2) = dq(ind_SN,2) + ((1+f_w_cell)*mSN(iSN)/56d0)*vload*v48(i,2)
           dq(ind_SN,3) = dq(ind_SN,3) + ((1+f_w_cell)*mSN(iSN)/56d0)*vload*v48(i,3)
        endif
     enddo

     ! remove gas from SN cells
     if(ind_blast(iSN)>0)then  ! if SN cells are inside the 'myid' domain
         d=uold(ind_blast(iSN),1)
         u=uold(ind_blast(iSN),2)/d
         v=uold(ind_blast(iSN),3)/d
         w=uold(ind_blast(iSN),4)/d
         ekk=0.5d0*d*(u*u+v*v+w*w)
         eint=uold(ind_blast(iSN),5)-ekk
         if(metal)Zload=uold(ind_blast(iSN),imetal)/d

         ! remove 48/56*mcell gas
         d=d-(mload-mSN(iSN)*f_load)/vol_loc  
         uold(ind_blast(iSN),1)=d
         uold(ind_blast(iSN),2)=d*u
         uold(ind_blast(iSN),3)=d*v
         uold(ind_blast(iSN),4)=d*w
         uold(ind_blast(iSN),5)=eint+0.5*d*(u*u + v*v + w*w)

         ! add f_leftover*SN ejecta to this cell 
         d = mSN(iSN)*f_leftover/vol_loc
         u = vSN(iSN,1)
         v = vSN(iSN,2)
         w = vSN(iSN,3)
         uold(ind_blast(iSN),1)=uold(ind_blast(iSN),1)+d
         uold(ind_blast(iSN),2)=uold(ind_blast(iSN),2)+d*u
         uold(ind_blast(iSN),3)=uold(ind_blast(iSN),3)+d*v
         uold(ind_blast(iSN),4)=uold(ind_blast(iSN),4)+d*w
         uold(ind_blast(iSN),5)=uold(ind_blast(iSN),5)+0.5*d*(u*u+v*v+w*w) 
         if(metal)then
            uold(ind_blast(iSN),imetal)=uold(ind_blast(iSN),imetal) &
                                      & - Zload*(mload-mSN(iSN)*f_load)/vol_loc &
                                      & + (yield + (1d0-yield)*ZSN(iSN))*mSN(iSN)*f_leftover/vol_loc
         endif
     endif
  enddo


#ifndef WITHOUTMPI
  if(global_search) then
     dq_mpi=0d0
     call MPI_ALLREDUCE(dq,dq_mpi,nSN_glo*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     dq=dq_mpi
  endif
#endif
 

  !---------------------------------------------------------------------
  ! Inject momenta from SNe
  !---------------------------------------------------------------------
  do iSN=1,nSN
     ind_SN= iSN_myid(iSN)
     dx=0.5d0**lSN(iSN)
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim   
     mload=mloadSN(ind_SN)
     do i=1,48
        if(icell_nbor(i,iSN)>0)then
           ! idx48 is measured for sub volume of vol_loc/8
           f_w_cell = (mload/48d0 + rho_nbor(i,iSN)*vol_loc/8d0) / (mSN(iSN)/56d0) - 1d0
           dloadloc = mload/48d0/vol_loc
           navg     = (dloadloc+rho_nbor(i,iSN))/2d0*scale_nH
           num_sn   = mSN(iSN)/(M_SNII/scale_msun)
           if(metal)then
              Zdepen = (dloadloc*Zload+rho_nbor(i,iSN)*Z_nbor(i,iSN))/(dloadloc+rho_nbor(i,iSN))
              Zdepen = (max(0.01,Zdepen/0.02))**(-0.28)
           else
              Zdepen = 1d0
           endif 

           f_w_crit = max(0d0,9d2/(f_esn*M_SNII)*num_sn**(-2d0/17d0)*navg**(-4d0/17d0)*Zdepen-1d0)
           if(f_w_cell.ge.f_w_crit)then !radiative phase
              vload    = dsqrt(2d0*f_esn*ESN*(1+f_w_crit)/(M_SNII*2d33))/scale_v/(1d0+f_w_cell)
           else
              vload    = dsqrt(2d0      *ESN/(1+f_w_cell)/(M_SNII*2d33))/scale_v
           endif

           icello=icell_nbor(i,iSN)
           d=uold(icello,1)      
           u=uold(icello,2)/d 
           v=uold(icello,3)/d 
           w=uold(icello,4)/d
           ekk=0.5d0*d*(u*u+v*v+w*w)
           eint=uold(icello,5)-ekk

           vol_nbor=(0.5d0**lv_nbor(i,iSN)*scale)**ndim
           pload(1) = (dq(ind_SN,1) + (1+f_w_cell)*mSN(iSN)*vload*v48(i,1))/56d0/vol_nbor ! double-checked
           pload(2) = (dq(ind_SN,2) + (1+f_w_cell)*mSN(iSN)*vload*v48(i,2))/56d0/vol_nbor
           pload(3) = (dq(ind_SN,3) + (1+f_w_cell)*mSN(iSN)*vload*v48(i,3))/56d0/vol_nbor

           dload=mload/48d0/vol_nbor
           uold(icello,1)=d+dload
           uold(icello,2)=uold(icello,2)+pload(1)
           uold(icello,3)=uold(icello,3)+pload(2)
           uold(icello,4)=uold(icello,4)+pload(3)
 
           ! energy conservation (if momentum vector cancels out, it naturally becomes thermal)
           u=pload(1)/dload 
           v=pload(2)/dload 
           w=pload(3)/dload
           uold(icello,5)=uold(icello,5)+0.5*dload*(u*u+v*v+w*w)
           u=uold(icello,2)/(d+dload) 
           v=uold(icello,3)/(d+dload) 
           w=uold(icello,4)/(d+dload)
           eint2 = uold(icello,5)- 0.5*(d+dload)*(u*u+v*v+w*w)
!           if(ind_blast(iSN)>0) &
!           print *,'#####', icello, ind_blast(iSN), uold(icello,1),uold(ind_blast(iSN),1) 
!write(1000+myid) sngl(pload/(d+dload)*scale_v/1d5),sngl(d*scale_nH),&
!     &sngl(dload*scale_nH),sngl(Zdepen),sngl(eint2/(d+dload)*scale_T2), &
!     &sngl(f_w_cell),sngl(f_w_crit),ind_SN,nstep,lSN(iSN),sngl(dq(ind_SN,1:3)/(d+dload)*scale_v/1d5)

           if(metal)then
               uold(icello,imetal)=uold(icello,imetal)+dload*ZloadSN(ind_SN)
           endif
        endif
     end do
  end do
 
  !if(.not.global_search) call clean_stop

end subroutine inject_momentum_SN
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_idx48(ista,iend,jsta,jend,ksta,kend,idx48,n48)
   implicit none
   integer::ista,iend,jsta,jend,ksta,kend,i,j,k,kk,n48
   integer, dimension(1:48)::idx48
   integer, dimension(1:4,1:4,1:4)::indall
   logical::error

   n48=0
   idx48=0
   error=.false.
   if(ista<0)ista=0
   if(jsta<0)jsta=0
   if(ksta<0)ksta=0
   if(iend>4)error=.true.
   if(jend>4)error=.true.
   if(kend>4)error=.true.
   if(iend>3)iend=3
   if(jend>3)jend=3
   if(kend>3)kend=3


   indall=reshape((/&
         -1, 1, 2,-1,   3, 4, 5, 6,   7, 8, 9,10,  -1,11,12,-1,&
         13,14,15,16,  17,-1,-1,18,  19,-1,-1,20,  21,22,23,24,&
         25,26,27,28,  29,-1,-1,30,  31,-1,-1,32,  33,34,35,36,&
         -1,37,38,-1,  39,40,41,42,  43,44,45,46, -1,47,48,-1/),(/4,4,4/))
   
   do k=ksta,kend
   do j=jsta,jend
   do i=ista,iend
      kk = indall(i+1,j+1,k+1)
      if(kk>0)then
         n48=n48+1
         idx48(n48)=kk
      endif
   end do
   end do
   end do 
    
end subroutine get_idx48
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_v48(v48)
   use amr_commons
   implicit none 
   real(dp)::v48(1:48,3),x,y,z,r
   integer::i,j,k,ind,indall(1:64)
   logical::ok

   ind = 0 
   indall=-1
   ! centre: [0,0,0]
   do k=1,4
   do j=1,4
   do i=1,4
      ok = .true.
      if((i==1.or.i==4).and.&
         (j==1.or.j==4).and.&
         (k==1.or.k==4)) ok=.false. ! edge
      if((i==2.or.i==3).and.&
         (j==2.or.j==3).and.&
         (k==2.or.k==3)) ok=.false. ! centre
      if(ok)then
         ind=ind+1
         x = (i-1)+0.5d0 - 2 
         y = (j-1)+0.5d0 - 2 
         z = (k-1)+0.5d0 - 2 
         r = dsqrt(dble(x*x+y*y+z*z))
         v48(ind,1) = x/r 
         v48(ind,2) = y/r 
         v48(ind,3) = z/r 
         indall(i+(j-1)*4+(k-1)*4*4) = ind 
      endif
   end do
   end do
   end do
end subroutine get_v48
!################################################################
!################################################################
!################################################################
!################################################################

subroutine checkSNboundary(xSN,ncpu_read,mylevel)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::mylevel,ncpu_read
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  integer::lmin,nx_loc,ilevel,lmax,bit_length,maxdom
  integer::imin,jmin,kmin,imax,jmax,kmax,ndom,impi,i,j
  integer,dimension(1:ncpu)::cpu_list
  logical,dimension(1:ncpu)::cpu_read
  real(dp)::scale,dx,dx_loc,drSN
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp)::xxmin,yymin,zzmin,xxmax,yymax,zzmax,dmax
  real(qdp),dimension(1:8)::bounding_min,bounding_max
  real(qdp)::dkey,order_min,oneqdp=1.0
  real(dp),dimension(1:3)::xSN 

 
  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=scale*0.5D0**mylevel

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  drSN=2d0*MAX(1.5d0*dx_loc*scale_l/aexp,rbubble*3.08d18) ! adaptive radius
  drSN=drSN/scale_l

  !-----------------------
  ! Map parameters
  !-----------------------
  lmax=nlevelmax
  cpu_read=.false.
  ! Compute boundaries for the SN cube of influence
  xxmin=(xSN(1)-drSN)/scale ; xxmax=(xSN(1)+drSN)/scale
  yymin=(xSN(2)-drSN)/scale ; yymax=(xSN(2)+drSN)/scale
  zzmin=(xSN(3)-drSN)/scale ; zzmax=(xSN(3)+drSN)/scale

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
  
end subroutine checkSNboundary
!################################################################
!################################################################
!################################################################
!################################################################
subroutine getSNonmyid2(iSN_myid,nSN_myid,xSN,nSN,lvSN)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,dimension(1:nSN)::iSN_myid,lvSN
  integer::nSN_myid,ii,nSN
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
  !drSN=2d0*MAX(1.5d0*dx_min*scale_l/aexp,rbubble*3.08d18)
  !drSN=drSN/scale_l

  !-----------------------
  ! Map parameters
  !-----------------------
  lmax=nlevelmax
  iSN_myid=0
  ii=0
  do iSN=1,nSN

     cpu_read=.false.
 
     drSN = dx_min*3*2d0**(nlevelmax-lvSN(iSN))
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

end subroutine getSNonmyid2
