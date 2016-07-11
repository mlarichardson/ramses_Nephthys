!####################################################################
!####################################################################
!####################################################################
subroutine mechanical_feedback_fine(ilevel,icount)
  use pm_commons
  use amr_commons
  use mechanical_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine computes the energy liberated from stellar winds
  ! based on the input library, and inject thermal or kinetic energy to 
  ! the surrounding of the young stars.
  ! This routine is called every fine time step.
  ! ind_pos_cell: position of the cell within an oct
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::npart1,npart2,icpu,icount,idim,ip
  integer::ind,ind_son,ind_cell,ilevel,iskip,info,nSNc,nSNc_mpi,nsn_star
  integer,dimension(1:nvector),save::ind_grid,ind_pos_cell
  real(dp)::tyoung,current_time,dteff
  real(dp)::skip_loc(1:3),scale,dx,dx_loc,vol_loc,x0(1:3)
  real(dp),dimension(1:twotondim,1:ndim),save::xc
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale_msun
  real(dp),dimension(1:twotondim),save::mw8, mzw8 ! SNe
  real(dp),dimension(1:twotondim,1:3),save:: pw8  ! SNe
  real(dp),dimension(1:nvector),save::mSN, mZSN
  real(dp),dimension(1:nvector,1:3),save::pSN
  real(dp)::mejecta
  real(dp),parameter::msun2g=2d33
  real(dp),parameter::myr2s=3.1536000d+13
  real(dp)::ttsta,ttend
  logical::ok,done_star
 
  if(icount==2) return
  if(.not.hydro) return
  if(ndim.ne.3)  return
  if(numbtot(1,ilevel)==0)return
  if(nstar_tot==0)return

#ifndef WITHOUTMPI
  if(myid.eq.1) ttsta=MPI_WTIME(info)
#endif 
  nSNc=0

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = scale_l**3*scale_d/msun2g

  ! Mesh spacing in that level
  call mesh_info (ilevel,skip_loc,scale,dx,dx_loc,vol_loc,xc) 

  ! To filter out old particles and compute individual time steps
  ! NB: the time step is always the coarser level time step, since no feedback for icount=2
  if (ilevel==levelmin)then
     dteff = dtold(ilevel)
  else
     dteff = dtold(ilevel-1)
  endif

  if (use_proper_time)then
     tyoung = t_delay*myr2s/(scale_t/aexp**2) 
     current_time=texp
     dteff = dteff*aexp**2
  else
     tyoung = t_delay*myr2s/scale_t 
     current_time=t
  endif
  tyoung = current_time - tyoung

  nSN_comm=0  ! important to initialize; number of communications (not SNe)
#ifndef WITHOUTMPI
  xSN_comm=0d0;ploadSN_comm=0d0;mSN_comm=0d0
  mloadSN_comm=0d0;mZloadSN_comm=0d0;iSN_comm=0
  floadSN_comm=0d0;eloadSN_comm=0d0
#endif

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ip=0

     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0

        ! Count star particles
        if(npart1>0)then
           do idim=1,ndim
              x0(idim) = xg(igrid,idim) -dx -skip_loc(idim) 
           end do
 
           ipart=headp(igrid)
    
           mw8=0d0;mzw8=0d0;pw8=0d0
    
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ok=.false.
              if(sn2_real_delay)then
                 ! if tp is younger than t_delay
                 if (idp(ipart).le.0.and.tp(ipart).ge.tyoung) then
                    call get_number_of_sn2  (tp(ipart), zp(ipart), idp(ipart),&
                             & mp0(ipart)*scale_msun,mp(ipart)*scale_msun,nsn_star,done_star)
                    if(nsn_star>0)ok=.true.
                 endif
              else ! single SN event
                 ! if tp is older than t_delay
                 if (idp(ipart).le.0.and.tp(ipart).le.tyoung)then
                    ok=.true.
                 endif
              endif

              if(ok)then
                 ! Find the cell index to get the position of it
                 ind_son=1
                 do idim=1,ndim
                    ind = int((xp(ipart,idim)/scale - x0(idim))/dx)
                    ind_son=ind_son+ind*2**(idim-1)
                 end do 
                 iskip=ncoarse+(ind_son-1)*ngridmax
                 ind_cell=iskip+igrid
                 if(son(ind_cell)==0)then  ! leaf cell
                    if(sn2_real_delay)then
                       mejecta = M_SNII/scale_msun*nsn_star
                    else
                       mejecta = mp(ipart)*eta_sn
                    endif
                    !mp_sun = mp(ipart)*scale_msun*eta_sn
                    mw8 (ind_son)  =mw8(ind_son)+mejecta
                    pw8 (ind_son,1)=pw8(ind_son,1)+mejecta*vp(ipart,1)
                    pw8 (ind_son,2)=pw8(ind_son,2)+mejecta*vp(ipart,2)
                    pw8 (ind_son,3)=pw8(ind_son,3)+mejecta*vp(ipart,3)
                    if(metal) mzw8(ind_son) = mzw8(ind_son) + mejecta*(zp(ipart)+(1d0-zp(ipart))*yield)
    
                    mp(ipart)=mp(ipart)-mejecta

                    if(sn2_real_delay) then
                       if(done_star) idp(ipart)=-idp(ipart) ! only if all SNe exploded
                    else
                       idp(ipart)=-idp(ipart)
                    endif
                 endif
              endif
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
    
           do ind=1,twotondim
              if (abs(mw8(ind))>0d0)then
                 nSNc=nSNc+1
                 ip=ip+1
                 ind_grid(ip)=igrid
                 ind_pos_cell(ip)=ind
     
                 ! collect information 
                 mSN(ip)=mw8(ind)
                 mZSN(ip)=mzw8(ind)
                 pSN(ip,1)=pw8(ind,1)
                 pSN(ip,2)=pw8(ind,2)
                 pSN(ip,3)=pw8(ind,3)
    
                 if(ip==nvector)then
                    call mech_fine(ind_grid,ind_pos_cell,ip,ilevel,mSN,pSN,mZSN,dteff)
                    ip=0
                 endif 
              endif
           enddo
    
    
        end if
        igrid=next(igrid)   ! Go to next grid
     end do ! End loop over grids
    
     if (ip>0) then
        call mech_fine(ind_grid,ind_pos_cell,ip,ilevel,mSN,pSN,mZSN,dteff)
        ip=0
     endif

  end do ! End loop over cpus


#ifndef WITHOUTMPI
  nSNc_mpi=0
  ! Deal with the stars around the bounary of each cpu (need MPI)
  call mech_fine_mpi(ilevel)
  call MPI_ALLREDUCE(nSNc,nSNc_mpi,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nSNc = nSNc_mpi
  if(myid.eq.1.and.nSNc>0.and.log_mfb) then
     ttend=MPI_WTIME(info)
     write(*,*) '--------------------------------------'
     write(*,*) 'Time elapsed in mechanical_fine [sec]:', sngl(ttend-ttsta), nSNc, nSN_comm 
     write(*,*) '--------------------------------------'
  endif
#endif

end subroutine mechanical_feedback_fine
!################################################################
!################################################################
!################################################################
subroutine mech_fine(ind_grid,ind_pos_cell,np,ilevel,mSN,pSN,mZSN,dteff)
  use amr_commons
  use pm_commons
  use hydro_commons
  use mechanical_commons
  implicit none
  integer::np,ilevel ! actually the number of cells
  integer,dimension(1:nvector)::ind_grid,ind_pos_cell
  real(dp),dimension(1:nvector)::mSN,mZSN,floadSN
  real(dp),dimension(1:nvector)::mloadSN,mZloadSN,eloadSN
  real(dp),dimension(1:nvector,1:3)::pSN,ploadSN
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine mechanical_feedback_fine 
  !-----------------------------------------------------------------------
  integer::i,j,nwco,idim,icell,igrid,ista,iend,ilevel2,ind_cell,ncell,irad
  real(dp)::d,u,v,w,e,z,eth,ekk,Tk,d0,u0,v0,w0,dteff
  real(dp)::dx,dx_loc,scale,vol_loc,nH_cen,tsim_yr,fleftSN
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale_msun,msun2g=2d33
  real(dp)::skip_loc(1:3),Tk0,ekk0,eth0,etot0
  real(dp),dimension(1:twotondim,1:ndim),save::xc
  ! Grid based arrays
  real(dp),dimension(1:ndim,1:nvector),save::xc2
  real(dp),dimension(1:nvector,1:nSNnei), save::p_solid,ek_solid
  real(dp)::d_nei,Z_nei,dm_ejecta,vol_nei
  real(dp)::mload,vload,Zload=0d0,f_esn2,eturb
  real(dp)::num_sn,nH_nei,Zdepen=1d0,f_w_cell,f_w_crit
  real(dp)::t_rad,r_rad,r_shell,m_cen,ekk_ej
  real(dp)::uavg,vavg,wavg,ul,vl,wl,ur,vr,wr
  real(dp)::d1,d2,d3,d4,d5,d6,dtot,pvar(1:ndim+3+512)
  real(dp)::vturb,vth,Mach,sig_s2,dratio,mload_cen
  ! For stars affecting across the boundary of a cpu
  integer, dimension(1:nSNnei),save::icpuSNnei
  integer ,dimension(1:nvector,0:twondim):: ind_nbor
  logical(dp),dimension(1:nvector,1:nSNnei),save ::snowplough

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = scale_l**3*scale_d/msun2g
  tsim_yr    = dteff*scale_t/365d0/3600d0/24d0

  ! Mesh variables
  call mesh_info (ilevel,skip_loc,scale,dx,dx_loc,vol_loc,xc)

  ! Record position of each cell [0.0-1.0] regardless of boxlen
  xc2=0d0
  do i=1,np
     do idim=1,ndim
        xc2(idim,i)=xg(ind_grid(i),idim)-skip_loc(idim)+xc(ind_pos_cell(i),idim)
     end do 
  end do

  !======================================================================
  ! Determine p_solid before redistributing mass 
  !   (momentum along some solid angle or cell) 
  ! - This way is desirable when two adjacent SNe explode simulataenously.
  ! - if the neighboring cell does not belong to myid, this will be done 
  !      in mech_fine_mpi
  !======================================================================
  p_solid=0d0;ek_solid=0d0;snowplough=.false.

  do i=1,np
     ind_cell = ncoarse+ind_grid(i)+(ind_pos_cell(i)-1)*ngridmax

     ! redistribute the mass/metals to the central cell
     call get_icell_from_pos (xc2(1:3,i), ilevel+1, igrid, icell,ilevel2)

     ! Sanity Check
     if((cpu_map(father(igrid)).ne.myid).or.&
        (ilevel.ne.ilevel2).or.&
        (ind_cell.ne.icell))     then 
        print *,'>>> fatal error in mech_fine'
        print *, cpu_map(father(igrid)),myid
        print *, ilevel, ilevel2
        print *, ind_cell, icell
        stop 
     endif
 
     num_sn    = mSN(i)/(M_SNII/scale_msun) ! doesn't have to be an integer
     nH_cen    = uold(icell,1)*scale_nH
     m_cen     = uold(icell,1)*vol_loc*scale_msun

     d   = uold(icell,1)
     u   = uold(icell,2)/d
     v   = uold(icell,3)/d
     w   = uold(icell,4)/d
     e   = uold(icell,5)
     e   = e-0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
     do irad=1,nener
        e = e - uold(icell,ndim+2+irad) 
     end do
#endif
     Tk  = e/d*scale_T2*(gamma-1.0)*0.6
    
     if(Tk<0)then
        print *,'TKERR : mech fbk (pre-call): TK<0', TK,icell
     endif

     !==========================================
     ! estimate floadSN(i)
     !==========================================
     ! notice that f_LOAD / f_LEFT is a global parameter for SN ejecta themselves!!!
     floadSN(i) = f_LOAD 

     if(log_mfb)then
398     format('MFB = ',f7.3,1x,f7.3,1x,f5.1,1x,f5.3,1x,f9.5,1x,f7.3)
        write(*,398) log10(d*scale_nH),log10(Tk),sngl(num_sn),floadSN(i),1./aexp-1,log10(dx_loc*scale_l/3.08d18)
     endif

     dm_ejecta = f_LOAD*mSN(i)/dble(nSNnei)  ! per solid angle
     mload     = f_LOAD*mSN(i) + uold(icell,1)*vol_loc*floadSN(i)  ! total SN ejecta + host cell
     if(metal) Zload = (f_LOAD*mZSN(i) + uold(icell,imetal)*vol_loc*floadSN(i))/mload

    
     do j=1,nSNnei
        call get_icell_from_pos (xc2(1:3,i)+xSNnei(1:3,j)*dx, ilevel+1, igrid, icell,ilevel2)
        if(cpu_map(father(igrid)).eq.myid) then ! if belong to myid

           Z_nei = z_ave*0.02 ! For metal=.false. 
           if(ilevel>ilevel2)then ! touching level-1 cells
              d_nei     = unew(icell,1)
              if(metal) Z_nei = unew(icell,imetal)/d_nei
           else
              d_nei     = uold(icell,1)
              if(metal) Z_nei = uold(icell,imetal)/d_nei
           endif
           f_w_cell  = (mload/dble(nSNnei) + d_nei*vol_loc/8d0)/dm_ejecta - 1d0
           nH_nei    = d_nei*scale_nH
           Zdepen    = (max(0.01,Z_nei/0.02))**(expZ_SN*2d0) 

           ! transition mass loading factor (momentum conserving phase)
           f_w_crit = (A_SN/1d4)**2d0/(f_esn*M_SNII)*num_sn**((expE_SN-1d0)*2d0)*nH_nei**(expN_SN*2d0)*Zdepen - 1d0
           f_w_crit = max(0d0,f_w_crit)

           if(f_w_cell.ge.f_w_crit)then ! radiative phase
              vload = dsqrt(2d0*f_esn*ESN*(1d0+f_w_crit)/(M_SNII*msun2g))/scale_v/(1d0+f_w_cell)/f_LOAD
              snowplough(i,j)=.true.
           else ! adiabatic phase
              f_esn2 = 1d0-(1d0-f_esn)*f_w_cell/f_w_crit
              vload = dsqrt(2d0*f_esn2*ESN/(1d0+f_w_cell)/(M_SNII*msun2g))/scale_v/f_LOAD
              snowplough(i,j)=.false.
           endif
           p_solid(i,j)=(1d0+f_w_cell)*dm_ejecta*vload
           ek_solid(i,j)=ek_solid(i,j)+p_solid(i,j)*(vload*f_LOAD)/2d0 !ek=(m*v)*v/2, not (d*v)*v/2

           if(log_mfb_mega)then
             write(*,'(" MFBN nHcen=", f6.2," nHnei=", f6.2, " mej=", f6.2, " mcen=", f6.2, " vload=", f6.2, " lv2=",I3," fwcrit=", f6.2, " fwcell=", f6.2, " psol=",f6.2, " mload/48=",f6.2," mnei/8=",f6.2," mej/48=",f6.2)') &
            & log10(nH_cen),log10(nH_nei),log10(mSN(i)*scale_msun),log10(m_cen),log10(vload*scale_v/1d5),ilevel2-ilevel,&
            & log10(f_w_crit),log10(f_w_cell),log10(p_solid(i,j)*scale_msun*scale_v/1d5),log10(mload*scale_msun/48),&
            & log10(d_nei*vol_loc/8d0*scale_msun),log10(dm_ejecta*scale_msun)
            endif

        endif
        
     enddo ! loop over neighboring cells
  enddo ! loop over SN cells

  !-----------------------------------------
  ! Redistribute mass from the SN cell
  !-----------------------------------------
  do i=1,np
     icell = ncoarse+ind_grid(i)+(ind_pos_cell(i)-1)*ngridmax
     d     = uold(icell,1)
     u     = uold(icell,2)/d
     v     = uold(icell,3)/d
     w     = uold(icell,4)/d
!     e     = uold(icell,5)
!     e     = e - 0.5*d*(u**2 + v**2 + w**2) ! total pressure (thermal+non-thermal)

     mloadSN (i) = mSN (i)*f_LOAD + d*vol_loc*floadSN(i)
     if(metal)then
        z = uold(icell,imetal)/d
        mZloadSN(i) = mZSN(i)*f_LOAD + d*z*vol_loc*floadSN(i)
     endif

     ! original momentum by star + gas entrained from the SN cell
     ploadSN(i,1) = pSN(i,1)*f_LOAD + vol_loc*d*u*floadSN(i)
     ploadSN(i,2) = pSN(i,2)*f_LOAD + vol_loc*d*v*floadSN(i)
     ploadSN(i,3) = pSN(i,3)*f_LOAD + vol_loc*d*w*floadSN(i)

     ! update the hydro variable
     fleftSN = 1d0 - floadSN(i)
     uold(icell,1) = mSN(i)  /vol_loc*f_LEFT + d*fleftSN
     uold(icell,2) = pSN(i,1)/vol_loc*f_LEFT + d*u*fleftSN ! make sure this is rho*v, not v
     uold(icell,3) = pSN(i,2)/vol_loc*f_LEFT + d*v*fleftSN
     uold(icell,4) = pSN(i,3)/vol_loc*f_LEFT + d*w*fleftSN
     if(metal)then
        uold(icell,imetal) = mZSN(i)/vol_loc*f_LEFT + d*z*fleftSN
     endif

     ! original kinetic energy of the gas entrained
     eloadSN(i) = 0.5d0*d*(u**2 + v**2 + w**2)*vol_loc*floadSN(i) 

     ! by commenting out the following two lines, 
     ! we are assuming that original thermal energy is left in the SN host cell
     ! original thermal energy of the gas entrained (including the non-thermal part)
     !eloadSN(i) = eloadSN(i) + e*vol_loc*floadSN(i)

     ! add the kinetic energy from SN
     d = mSN(i)/vol_loc
     u = pSN(i,1)/mSN(i)
     v = pSN(i,2)/mSN(i)
     w = pSN(i,3)/mSN(i)
     uold(icell,5) = uold(icell,5) + 0.5d0*d*(u**2 + v**2 + w**2)*f_LEFT
 
     ! original energy of the SN particle
     eloadSN(i) = eloadSN(i) + 0.5d0*mSN(i)*(u**2 + v**2 + w**2)*f_LOAD

     ! update ek_solid     
     ek_solid(i,:) = ek_solid(i,:) + eloadSN(i)/dble(nSNnei)

  enddo  ! loop over SN cell


  !-------------------------------------------------------------
  ! Find and save stars affecting across the boundary of a cpu
  !-------------------------------------------------------------
  do i=1,np

     ind_cell = ncoarse+ind_grid(i)+(ind_pos_cell(i)-1)*ngridmax

     nwco=0; icpuSNnei=0
     do j=1,nSNnei
        call get_icell_from_pos (xc2(1:3,i)+xSNnei(1:3,j)*dx, ilevel+1, igrid, icell,ilevel2)
 
        if(cpu_map(father(igrid)).ne.myid) then ! need mpi
           nwco=nwco+1
           icpuSNnei(nwco)=cpu_map(father(igrid))
        else  ! can be handled locally
           vol_nei = vol_loc*(2d0**ndim)**(ilevel-ilevel2)
           pvar(:) = 0d0 ! temporary primitive variable
           if(ilevel>ilevel2)then ! touching level-1 cells
              pvar(1:5) = unew(icell,1:5)
              if(metal) pvar(imetal) = unew(icell,imetal)
#if NENER>0
              do irad=1,nener
                 pvar(ndim+2+irad) = unew(icell,ndim+2+irad)
              enddo
#endif
           else
              pvar(1:5) = uold(icell,1:5)
              if(metal) pvar(imetal) = uold(icell,imetal)
#if NENER>0
              do irad=1,nener
                 pvar(ndim+2+irad) = uold(icell,ndim+2+irad)
              enddo
#endif
           endif 

           d0=pvar(1)
           u0=pvar(2)/d0
           v0=pvar(3)/d0
           w0=pvar(4)/d0
           ekk0=0.5d0*d0*(u0**2+v0**2+w0**2)
           eth0=pvar(5)-ekk0
#if NENER>0
           do irad=1,nener
              eth0=eth0-pvar(ndim+2+irad)
           end do
#endif

           Tk0 = eth0/d0*scale_T2*(gamma-1)*0.6
           if(Tk0<0)then
              print *, 'TKERR: mech (post-call): Tk<0 =', Tk0
              stop
           endif 

           d= mloadSN(i  )/dble(nSNnei)/vol_nei
           u=(ploadSN(i,1)/dble(nSNnei)+p_solid(i,j)*vSNnei(1,j))/vol_nei/d
           v=(ploadSN(i,2)/dble(nSNnei)+p_solid(i,j)*vSNnei(2,j))/vol_nei/d
           w=(ploadSN(i,3)/dble(nSNnei)+p_solid(i,j)*vSNnei(3,j))/vol_nei/d
           pvar(1)=pvar(1)+d
           pvar(2)=pvar(2)+d*u
           pvar(3)=pvar(3)+d*v
           pvar(4)=pvar(4)+d*w
           ekk_ej = 0.5*d*(u**2 + v**2 + w**2)   
           etot0  = eth0+ekk0+ek_solid(i,j)/vol_nei

           ! the minimum thermal energy input floor
           d   = pvar(1)
           u   = pvar(2)/d
           v   = pvar(3)/d
           w   = pvar(4)/d
           ekk = 0.5*d*(u**2 + v**2 + w**2)
           eth = d*Tk0/(scale_T2*(gamma-1)*0.6)
           if (etot0.lt.(ekk+eth0)) then
              pvar(5) = ekk0 + eth0 + ekk_ej
           else
              pvar(5) = ekk0 + eth0 + ek_solid(i,j)/vol_nei ! additional energy from SNe+entrained gas
           endif
           eturb = pvar(5) - ekk - eth0
#if NENER>0
           do irad=1,nener
              pvar(5) = pvar(5) + pvar(ndim+2+irad)
           end do
#endif

           if(metal)then
               pvar(imetal)=pvar(imetal)+mzloadSN(i)/dble(nSNnei)/vol_nei
           end if

           ! update the hydro variable
           if(ilevel>ilevel2)then ! touching level-1 cells
              unew(icell,1:5) = pvar(1:5)
              if(metal) unew(icell,imetal) = pvar(imetal)
#if NENER>0
              do irad=1,nener
                 unew(icell,ndim+2+irad) = pvar(ndim+2+irad)
              end do
#endif
           else
              uold(icell,1:5) = pvar(1:5)
              if(metal) uold(icell,imetal) = pvar(imetal)
#if NENER>0
              do irad=1,nener
                 uold(icell,ndim+2+irad) = pvar(ndim+2+irad)
              end do
#endif
           endif 

        end if
     end do ! loop over 48 neighbors


#ifndef WITHOUTMPI 
     if(nwco>0)then  ! for SNs across different cpu
        if(nwco>1)then
           ! remove redundant cpu list for this SN cell
           call redundant_non_1d(icpuSNnei(1:nwco), nwco, nwco)
        endif
        ista=nSN_comm+1
        iend=ista+nwco-1

        if(iend>ncomm_max)then
           write(*,*) 'Error: increase ncomm_max in mechanical_fine.f90', ncomm_max, iend
           call clean_stop
        endif
        iSN_comm (ista:iend)=icpuSNnei(1:nwco)
        !lSN_comm (ista:iend)=ilevel
        mSN_comm (ista:iend )=mSN(i)
        mloadSN_comm (ista:iend  )=mloadSN(i)
        xSN_comm (1,ista:iend)=xc2(1,i)
        xSN_comm (2,ista:iend)=xc2(2,i)
        xSN_comm (3,ista:iend)=xc2(3,i)
        ploadSN_comm (1,ista:iend)=ploadSN(i,1)
        ploadSN_comm (2,ista:iend)=ploadSN(i,2)
        ploadSN_comm (3,ista:iend)=ploadSN(i,3)
        floadSN_comm (ista:iend  )=floadSN(i)
        eloadSN_comm (ista:iend  )=eloadSN(i)
        if(metal) mZloadSN_comm(ista:iend  )=mZloadSN(i)
        nSN_comm=nSN_comm+nwco

     endif
#endif

  end do ! loop over SN cell

end subroutine mech_fine
!################################################################
!################################################################
!################################################################
subroutine mech_fine_mpi(ilevel)
  use amr_commons
  use mechanical_commons
  use pm_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::i,j,info,nSN_tot,icpu,ncpu_send,ncpu_recv,ncc
  integer::ncell_recv,ncell_send,cpu2send,cpu2recv,tag,np
  integer::isend_sta,irecv_sta,irecv_end
  real(dp),dimension(:,:),allocatable::SNsend,SNrecv,p_solid,ek_solid
  integer ,dimension(:),allocatable::list2recv,list2send
  integer, dimension(:),allocatable::reqrecv,reqsend
  integer, dimension(:,:),allocatable::statrecv,statsend
  ! SN variables
  real(dp)::dx,dx_loc,scale,vol_loc
  real(dp)::mloadSN_i,zloadSN_i,ploadSN_i(1:3),mSN_i,xSN_i(1:3),fload_i
  real(dp)::f_esn2,d_nei,Z_nei,f_w_cell,f_w_crit,nH_nei
  real(dp)::num_sn,vload,Tk,Zdepen=1d0,vol_nei,dm_ejecta
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale_msun,msun2g=2d33, pvar(1:ndim+3+512),eturb,etot0
  real(dp)::skip_loc(1:3),d,u,v,w,ekk,eth,d0,u0,v0,w0,eth0,ekk0,Tk0,ekk_ej
  integer::igrid,icell,ilevel,ilevel2,irad
  real(dp),dimension(1:twotondim,1:ndim),save::xc
  logical(dp),dimension(1:nvector,1:nSNnei),save ::snowplough

  if(ndim.ne.3) return

  snowplough = .false.

  !============================================================
  ! For MPI communication
  !============================================================
  ncpu_send=0;ncpu_recv=0

  nSN_comm_cpu=0 
  nSN_comm_mpi=0
  nSN_comm_mpi(myid)=nSN_comm
  ! compute the total number of communications needed
  call MPI_ALLREDUCE(nSN_comm_mpi,nSN_comm_cpu,ncpu,&
                   & MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nSN_tot = sum(nSN_comm_cpu)
  if(nSN_tot==0) return


  allocate(icpuSN_comm    (1:nSN_tot,1:2))
  allocate(icpuSN_comm_mpi(1:nSN_tot,1:2))

  ! index for mpi variable
  if(myid==1)then
     isend_sta = 0 
  else
     isend_sta = sum(nSN_comm_cpu(1:myid-1)) 
  endif

  icpuSN_comm=0
  do i=1,nSN_comm_cpu(myid)
     icpuSN_comm(isend_sta+i,1)=myid
     icpuSN_comm(isend_sta+i,2)=iSN_comm(i)
     ! iSN_comm:   local variable
     ! icpuSN_comm:  local (but extended) variable to be passed to a mpi variable
     ! icpuSN_comm_mpi: mpi variable
  end do

  ! share the list of communications
  icpuSN_comm_mpi=0
  call MPI_ALLREDUCE(icpuSN_comm,icpuSN_comm_mpi,nSN_tot*2,&
                   & MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)

  ncell_send = nSN_comm_cpu(myid)
  ncell_recv = count(icpuSN_comm_mpi(:,2).eq.myid, 1)

  ! check if myid needs to send anything
  if(ncell_send>0)then
     allocate( SNsend (1:nvarSN,1:ncell_send)) ! x(3),m,mz,p,pr
     allocate( list2send(1:ncell_send)    )
     list2send=0;  SNsend=0d0
     list2send=icpuSN_comm_mpi(isend_sta+1:isend_sta+ncell_send,2)
     ncpu_send=1
     if(ncell_send>1) call redundant_non_1d (list2send,ncell_send,ncpu_send)
     ! ncpu_send = No. of cpus to which myid should send info 
     allocate( reqsend  (1:ncpu_send) )
     allocate( statsend (1:MPI_STATUS_SIZE,1:ncpu_send) )
     reqsend=0; statsend=0
  endif

  ! check if myid needs to receive anything
  if(ncell_recv>0)then
     allocate( SNrecv (1:nvarSN,1:ncell_recv)) ! x(3),m,mz,p,pr
     allocate( list2recv(1:ncell_recv)    )
     list2recv=0;  SNrecv=0d0
     j=0
     do i=1,nSN_tot
        if(icpuSN_comm_mpi(i,2).eq.myid)then
           j=j+1 
           list2recv(j) = icpuSN_comm_mpi(i,1)
        endif
     end do
 
     ncc = j
     if(ncc.ne.ncell_recv)then ! sanity check
        write(*,*) 'Error in mech_fine_mpi: ncc != ncell_recv',ncc,ncell_recv,myid
        call clean_stop
     endif

     ncpu_recv=1 ! No. of cpus from which myid should receive info 
     if(j>1) call redundant_non_1d(list2recv,ncc,ncpu_recv)

     allocate( reqrecv  (1:ncpu_recv) )
     allocate( statrecv (1:MPI_STATUS_SIZE,1:ncpu_recv) )
     reqrecv=0; statrecv=0
  endif

  ! prepare one variable and send
  if(ncell_send>0)then
     do icpu=1,ncpu_send
        cpu2send = list2send(icpu)
        ncc=0 ! number of SN host cells that need communications with myid=cpu2send
        do i=1,ncell_send
           j=i+isend_sta
           if(icpuSN_comm_mpi(j,2).eq.cpu2send)then
              ncc=ncc+1
              SNsend(1:3,ncc)=xSN_comm (1:3,i)
              SNsend(4  ,ncc)=mSN_comm (    i)
              SNsend(5  ,ncc)=mloadSN_comm (i)
              SNsend(6:8,ncc)=ploadSN_comm (1:3,i)
              SNsend(9  ,ncc)=floadSN_comm (i)
              SNsend(10 ,ncc)=eloadSN_comm (i)
              if(metal)SNsend(11,ncc)=mZloadSN_comm(i)
           endif
        end do ! i

        tag = myid + cpu2send + ncc
        call MPI_ISEND (SNsend(1:nvarSN,1:ncc),ncc*nvarSN,MPI_DOUBLE_PRECISION, &
                      & cpu2send-1,tag,MPI_COMM_WORLD,reqsend(icpu),info) 
     end do ! icpu

  endif ! ncell_send>0


  ! receive one large variable
  if(ncell_recv>0)then
     irecv_sta=1
     do icpu=1,ncpu_recv
        cpu2recv = list2recv(icpu)
        ncc=0 ! number of SN host cells that need communications with cpu2recv
        do i=1,nSN_tot
           if((icpuSN_comm_mpi(i,1)==cpu2recv).and.&
             &(icpuSN_comm_mpi(i,2)==myid)      )then
              ncc=ncc+1
           endif
        end do
        irecv_end=irecv_sta+ncc-1
        tag = myid + cpu2recv + ncc
        
! NB. Some mpi routines do not like the receiving buffer
        call MPI_IRECV (SNrecv(1:nvarSN,irecv_sta:irecv_end),ncc*nvarSN,MPI_DOUBLE_PRECISION,&
                     & cpu2recv-1,tag,MPI_COMM_WORLD,reqrecv(icpu),info)

        

        irecv_sta=irecv_end+1
     end do ! icpu 

  endif ! ncell_recv >0

  if(ncpu_send>0)call MPI_WAITALL(ncpu_send,reqsend,statsend,info)
  if(ncpu_recv>0)call MPI_WAITALL(ncpu_recv,reqrecv,statrecv,info)


  !============================================================
  ! inject mass/metal/momentum
  !============================================================

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = scale_l**3*scale_d/msun2g

  ! Mesh variables
  call mesh_info (ilevel,skip_loc,scale,dx,dx_loc,vol_loc,xc)

  np = ncell_recv
  if(ncell_recv>0) then
     allocate(p_solid(1:np,1:nSNnei))
     allocate(ek_solid(1:np,1:nSNnei))
     p_solid=0d0;ek_solid=0d0
  endif

  ! Compute the momentum first before redistributing mass
  do i=1,np

     xSN_i(1:3) = SNrecv(1:3,i)
     mSN_i      = SNrecv(4,i)
     mloadSN_i  = SNrecv(5,i)
     fload_i    = SNrecv(9,i)
     dm_ejecta  = f_LOAD*mSN_i/dble(nSNnei)
     num_sn     = mSN_i/(M_SNII/scale_msun)
     if(metal) ZloadSN_i = SNrecv(11,i)/SNrecv(5,i)
     ek_solid(i,:) = SNrecv(10,i)/dble(nSNnei) ! kinetic energy of the gas mass entrained from the host cell + SN
 

     do j=1,nSNnei
        call get_icell_from_pos (xSN_i+xSNnei(1:3,j)*dx, ilevel+1, igrid, icell,ilevel2)
        if(cpu_map(father(igrid)).eq.myid) then ! if belong to myid
           Z_nei = z_ave*0.02 ! for metal=.false.
           if(ilevel>ilevel2)then ! touching level-1 cells
              d_nei     = unew(icell,1)
              if(metal) Z_nei = unew(icell,imetal)/d_nei
           else
              d_nei     = uold(icell,1)
              if(metal) Z_nei = uold(icell,imetal)/d_nei
           endif
           f_w_cell = (mloadSN_i/dble(nSNnei) + d_nei*vol_loc/8d0)/dm_ejecta - 1d0
           nH_nei   = d_nei*scale_nH
           Zdepen   =(max(0.01,Z_nei/0.02))**(expZ_SN*2d0) !From Thornton+(98)

           ! transition mass loading factor (momentum conserving phase)
           f_w_crit = (A_SN/1d4)**2d0/(f_esn*M_SNII)*num_sn**((expE_SN-1d0)*2d0)*nH_nei**(expN_SN*2d0)*Zdepen - 1d0
           f_w_crit = max(0d0,f_w_crit)

           if(f_w_cell.ge.f_w_crit)then ! radiative phase
              vload = dsqrt(2d0*f_esn*ESN*(1d0+f_w_crit)/(M_SNII*msun2g))/scale_v/(1d0+f_w_cell)/f_LOAD
              snowplough(i,j)=.true.
           else
              f_esn2 = 1d0-(1d0-f_esn)*f_w_cell/f_w_crit
              vload = dsqrt(2d0*f_esn2*ESN/(1d0+f_w_cell)/(M_SNII*msun2g))/scale_v/f_LOAD
              snowplough(i,j)=.false.
           endif
           p_solid(i,j)=(1d0+f_w_cell)*dm_ejecta*vload
           ek_solid(i,j)=ek_solid(i,j)+p_solid(i,j)*(vload*f_LOAD)/2d0 !ek
        endif
       
     enddo ! loop over neighboring cells

  end do



  ! Apply the SNe to this cpu domain
  do i=1,np

     xSN_i(1:3) = SNrecv(1:3,i)
     mSN_i      = SNrecv(4,i)
     mloadSN_i  = SNrecv(5,i)
     ploadSN_i(1:3)=SNrecv(6:8,i)
     if(metal) ZloadSN_i = SNrecv(11,i)/SNrecv(5,i)

     do j=1,nSNnei
 
        
        call get_icell_from_pos (xSN_i(1:3)+xSNnei(1:3,j)*dx, ilevel+1, igrid, icell, ilevel2)
        if(cpu_map(father(igrid))==myid)then

           vol_nei = vol_loc*(2d0**ndim)**(ilevel-ilevel2)
           if(ilevel>ilevel2)then ! touching level-1 cells
              pvar(1:5) = unew(icell,1:5)
              if(metal) pvar(imetal) = unew(icell,imetal)
#if NENER>0 
              do irad=1,nener
                 pvar(ndim+2+irad) = unew(icell,ndim+2+irad)
              end do
#endif
           else
              pvar(1:5) = uold(icell,1:5)
              if(metal) pvar(imetal) = uold(icell,imetal)
#if NENER>0 
              do irad=1,nener
                 pvar(ndim+2+irad) = uold(icell,ndim+2+irad)
              end do
#endif
           endif
     
           d0=pvar(1)
           u0=pvar(2)/d0
           v0=pvar(3)/d0
           w0=pvar(4)/d0
           ekk0=0.5d0*d0*(u0**2d0 + v0**2d0 + w0**2d0)
           eth0=pvar(5)-ekk0
#if NENER>0
           do irad=1,nener
              eth0 = eth0 - pvar(ndim+2+irad)
           end do
#endif
           Tk0 = eth0/d0*scale_T2*(gamma-1)*0.6

           if(Tk0<0) then
              print *,'TKERR : part1 mpi', myid,Tk0,d0*scale_nH
              stop
           endif

           d= mloadSN_i   /dble(nSNnei)/vol_nei
           u=(ploadSN_i(1)/dble(nSNnei)+p_solid(i,j)*vSNnei(1,j))/vol_nei/d
           v=(ploadSN_i(2)/dble(nSNnei)+p_solid(i,j)*vSNnei(2,j))/vol_nei/d
           w=(ploadSN_i(3)/dble(nSNnei)+p_solid(i,j)*vSNnei(3,j))/vol_nei/d

           pvar(1)=pvar(1)+d
           pvar(2)=pvar(2)+d*u
           pvar(3)=pvar(3)+d*v
           pvar(4)=pvar(4)+d*w
           ekk_ej = 0.5*d*(u**2 + v**2 + w**2)
           etot0  = eth0+ekk0+ek_solid(i,j)/vol_nei  

           ! the minimum thermal energy input floor
           d   = pvar(1)
           u   = pvar(2)/d
           v   = pvar(3)/d
           w   = pvar(4)/d
           ekk = 0.5*d*(u**2 + v**2 + w**2)
           eth = d*Tk0/(scale_T2*(gamma-1)*0.6)

           if (etot0.lt.(ekk+eth0)) then
              pvar(5) = ekk0 + eth0 + ekk_ej
           else
              pvar(5) = ekk0 + eth0 + ek_solid(i,j)/vol_nei ! additional energy from SNe+entrained gas
           endif
           eturb = pvar(5) - ekk - eth0
#if NENER>0
           do irad=1,nener
              pvar(5) = pvar(5) + pvar(ndim+2+irad)
           end do
#endif

           if(metal)then
               pvar(imetal)=pvar(imetal)+mloadSN_i/dble(nSNnei)*ZloadSN_i/vol_nei
           end if

           ! update the hydro variable
           if(ilevel>ilevel2)then ! touching level-1 cells
              unew(icell,1:5) = pvar(1:5)
              if(metal) unew(icell,imetal) = pvar(imetal)
#if NENER>0
              do irad=1,nener
                 unew(icell,ndim+2+irad) = pvar(ndim+2+irad)
              end do
#endif
           else
              uold(icell,1:5) = pvar(1:5)
              if(metal) uold(icell,imetal) = pvar(imetal)
#if NENER>0
              do irad=1,nener
                 uold(icell,ndim+2+irad) = pvar(ndim+2+irad)
              end do
#endif
           endif
 
        endif ! if this belongs to me


     end do ! loop over neighbors
  end do ! loop over SN cells


  deallocate(icpuSN_comm_mpi, icpuSN_comm)
  if(ncell_send>0) deallocate(list2send,SNsend,reqsend,statsend)
  if(ncell_recv>0) deallocate(list2recv,SNrecv,reqrecv,statrecv)
  if(ncell_recv>0) deallocate(p_solid,ek_solid)

  nSN_comm=nSN_tot
#endif

end subroutine mech_fine_mpi
!################################################################
!################################################################
!################################################################
subroutine get_number_of_sn2(birth_time,zp_star,id_star,mass0,mass1,nsn,done_star)
  use amr_commons, ONLY:dp,M_SNII,eta_sn,sn2_real_delay
  use random
  implicit none
  real(kind=dp)::birth_time,zp_star,mass0,mass1 ! birth_time in code, mass in Msun
  integer::nsn,nsn_tot,nsn_sofar,nsn_ok,nsn_try
  integer::i,localseed,id_star   ! necessary for the random number
  real(kind=dp)::age_star,tMyr,xdum,ydum,logzpsun
  real(kind=dp)::co0,co1,co2
  ! fit to the cumulative number fraction for Kroupa IMF
  ! ~kimm/soft/starburst99/output_hires_new/snrate.pro 
  ! ~kimm/soft/starburst99/output_hires_new/fit.pro 
  ! ~kimm/soft/starburst99/output_hires_new/sc.pro 
  real(kind=dp),dimension(1:3)::coa=(/-2.677292E-01,1.392208E-01,-5.747332E-01/)
  real(kind=dp),dimension(1:3)::cob=(/4.208666E-02, 2.152643E-02, 7.893866E-02/)
  real(kind=dp),dimension(1:3)::coc=(/-8.612668E-02,-1.698731E-01,1.867337E-01/)
  real(kind=dp),external::ran1 
  logical:: done_star


  ! determine the SNrateCum curve for Zp
  logzpsun=max(min(zp_star,0.05),0.008) ! no extrapolation
  logzpsun=log10(logzpsun/0.02)

  co0 = coa(1)+coa(2)*logzpsun+coa(3)*logzpsun**2
  co1 = cob(1)+cob(2)*logzpsun+cob(3)*logzpsun**2
  co2 = coc(1)+coc(2)*logzpsun+coc(3)*logzpsun**2

  ! RateCum = co0+sqrt(co1*Myr+co2)

  ! get stellar age
  call getStarAgeGyr(birth_time,age_star)

  tMyr=age_star*1d3 !Myr
  nsn=0
  if(tMyr.le.(-co2/co1))then
     return
  endif

  ! total number of SNe
  nsn_tot   = NINT(mass0*(eta_sn/M_SNII),kind=4)
  nsn_sofar = NINT((mass0-mass1)/M_SNII ,kind=4)

  if(sn2_real_delay)then
     nsn_try=nsn_tot
  else
     nsn_try=1
  endif

  nsn_ok    = 0
  localseed = -abs(id_star) !make sure that the number is negative

  do i=1,nsn_try 
     xdum =  ran1(localseed)
     ! inverse function for y=co0+sqrt(co1*x+co2)
     ydum = ((xdum-co0)**2.-co2)/co1
     if(ydum.le.tMyr) then
        nsn_ok=nsn_ok+1
     endif
  end do

  if(.not.sn2_real_delay) nsn_ok=nsn_ok*nsn_tot

  nsn = nsn_ok - nsn_sofar
  if(nsn_tot.eq.0)then
      write(*,*) 'Fatal error: please increase the mass of your star particle'
      stop
  endif

  done_star=.false.
  if(nsn_ok.eq.nsn_tot) done_star=.true.

end subroutine get_number_of_sn2
!################################################################
!################################################################
!################################################################
!################################################################
function ran1(idum)
   implicit none
   integer:: idum,IA,IM,IQ,IR,NTAB,NDIV
   real(kind=8):: ran1,AM,EPS,RNMX
   parameter(IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,&
            &NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
   integer::j,k,iv(NTAB),iy
   save iv,iy
   data iv /NTAB*0/, iy /0/ 
   if (idum.le.0.or.iy.eq.0) then! initialize
      idum=max(-idum,1)
      do j=NTAB+8,1,-1
         k=idum/IQ
         idum=IA*(idum-k*IQ)-IR*k
         if (idum.lt.0) idum=idum+IM
         if (j.le.NTAB) iv(j)=idum
      end do
      iy=iv(1)
   end if
   k=idum/IQ
   idum=IA*(idum-k*IQ)-IR*k
   if (idum.lt.0) idum=idum+IM
   j=1+iy/NDIV
   iy=iv(j)
   iv(j)=idum
   ran1=min(AM*iy,RNMX)
   return
end
!################################################################
!################################################################
!################################################################
!################################################################
subroutine getStarAgeGyr(birth_time,age_star)
   use amr_commons
   implicit none
   real(dp)::age_star,birth_time

   if (use_proper_time)then
      call getAgeGyr    (birth_time, age_star)
   else
      call getProperTime(birth_time, age_star)
      call getAgeGyr    (age_star,   age_star)
   endif

end subroutine getStarAgeGyr
!################################################################
!################################################################
!################################################################
subroutine redundant_non_1d(list,ndata,ndata2)
   implicit none
   integer,dimension(1:ndata)::list,list2
   integer,dimension(1:ndata)::ind
   integer::ndata,i,y1,ndata2

   ! sort first
   call heapsort(list,ind,ndata) 

   list(:)=list(ind)

   ! check the redundancy
   list2(:) = list(:)

   y1=list(1)
   ndata2=1
   do i=2,ndata
       if(list(i).ne.y1)then
          ndata2=ndata2+1
          list2(ndata2)=list(i)
          y1=list(i)
       endif
   end do

   list =list2

end subroutine redundant_non_1d
!################################################################
!################################################################
!################################################################
subroutine heapsort(ain,ind,n)
   implicit none
   integer::n
   integer,dimension(1:n)::ain,aout,ind
   integer::i,j,l,ir,idum,rra

   l=n/2+1
   ir=n
   do i=1,n
      aout(i)=ain(i)                        ! Copy input array to output array
      ind(i)=i                                   ! Generate initial idum array
   end do
   if(n.eq.1) return                            ! Special for only one record
10 continue
   if(l.gt.1)then
      l=l-1
      rra=aout(l)
      idum=ind(l)
   else
      rra=aout(ir)
      idum=ind(ir)
      aout(ir)=aout(1)
      ind(ir)=ind(1)
      ir=ir-1
      if(ir.eq.1)then
        aout(1)=rra
        ind(1)=idum
        return
      endif
    endif
    i=l
    j=l+l
20  if(j.le.ir)then
       if(j.lt.ir)then
          if(aout(j).lt.aout(j+1))j=j+1
       endif
       if(rra.lt.aout(j))then
          aout(i)=aout(j)
          ind(i)=ind(j)
          i=j
          j=j+j
       else
          j=ir+1
       endif
       go to 20
    endif
    aout(i)=rra
    ind(i)=idum
    go to 10
end subroutine heapsort
!################################################################
!################################################################
!################################################################
subroutine mesh_info (ilevel,skip_loc,scale,dx,dx_loc,vol_loc,xc)
  use amr_commons
  implicit none
  integer::ilevel,ind,ix,iy,iz,nx_loc
  real(dp)::skip_loc(1:3),scale,dx,dx_loc,vol_loc
  real(dp),dimension(1:twotondim,1:ndim):: xc

  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
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

end subroutine
!################################################################
!################################################################
!################################################################
subroutine get_icell_from_pos (fpos,ilevel_max,ind_grid,ind_cell,ilevel_out)
  use amr_commons
  implicit none
  real(dp)::fpos(1:3)
  integer ::ind_grid,ind_cell
  !-------------------------------------------------------------------
  ! This routnies find the index of the leaf cell for a given position
  ! fpos: positional info, from [0.0-1.0] (scaled by scale)
  ! ilevel_max: maximum level you want to search
  ! ind_cell: index of the cell
  ! ind_grid: index of the grid that contains the cell
  ! ilevel_out: level of this cell
  ! You can check whether this grid belongs to this cpu
  !     by asking cpu_map(father(ind_grid))
  !-------------------------------------------------------------------
  integer::ilevel_max
  integer::ilevel_out,nx_loc,i,ind,idim
  real(dp)::scale,dx,fpos2(1:3)
  real(dp)::skip_loc(1:3),x0(1:3)
  logical ::not_found

  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
 
  fpos2=fpos
  if (fpos2(1).gt.boxlen) fpos2(1)=fpos2(1)-boxlen
  if (fpos2(2).gt.boxlen) fpos2(2)=fpos2(2)-boxlen
  if (fpos2(3).gt.boxlen) fpos2(3)=fpos2(3)-boxlen
  if (fpos2(1).lt.0d0) fpos2(1)=fpos2(1)+boxlen
  if (fpos2(2).lt.0d0) fpos2(2)=fpos2(2)+boxlen
  if (fpos2(3).lt.0d0) fpos2(3)=fpos2(3)+boxlen

  not_found=.true.
  ind_grid =1  ! this is level=1 grid
  ilevel_out=1
  do while (not_found)
     dx = 0.5D0**ilevel_out
     x0 = xg(ind_grid,1:ndim)-dx-skip_loc  !left corner of *this* grid [0-1]

     ind  = 1 
     do idim = 1,ndim
        i = int( (fpos2(idim) - x0(idim))/dx)
        ind = ind + i*2**(idim-1)
     end do

     ind_cell = ind_grid+ncoarse+(ind-1)*ngridmax
!     write(*,'(2(I2,1x),2(I10,1x),3(f10.8,1x))') ilevel_out,ilevel_max,ind_grid,ind_cell,fpos2
     if(son(ind_cell)==0.or.ilevel_out==ilevel_max) return

     ind_grid=son(ind_cell)
     ilevel_out=ilevel_out+1
  end do

end subroutine
!################################################################
!################################################################
!################################################################
