!####################################################################
!####################################################################
!####################################################################
subroutine mechanical_feedback_snIa_fine(ilevel,icount)
  use pm_commons
  use amr_commons
  use mechanical_commons
  use hydro_commons,only:uold
#ifdef RT
  use rt_parameters,only:group_egy
  use SED_module,only:nSEDgroups,inp_SED_table
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine computes the energy liberated from supernova Ia,
  ! and inject momentum and energy to the surroundings of the young stars.
  ! This routine is called every fine time step.
  ! ind_pos_cell: position of the cell within an oct
  ! m8,mz8,nph8: temporary variable necessary for an oct to add up 
  !              the mass, metal, etc. on a cell by cell basis
  ! mejecta: mass of the ejecta from SNe
  ! Zejecta: metallicity of the ejecta (not yield)
  ! mZSNe : total metal mass from SNe in each cell
  ! mchSNe: total mass of each chemical element in each cell
  ! nphSNe: total production rate of ionising radiation in each cell.
  !         This is necessary to estimate the Stromgren sphere
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::npart1,npart2,icpu,icount,idim,ip,ich
  integer::ind,ind_son,ind_cell,ilevel,iskip,info,nSNc,nSNc_mpi
  integer,dimension(1:nvector),save::ind_grid,ind_pos_cell
  real(dp)::nsnIa_star,mass0,nsnIa_tot,nsnIa_mpi
  real(dp)::tyoung,current_time,dteff
  real(dp)::skip_loc(1:3),scale,dx,dx_loc,vol_loc,x0(1:3)
  real(dp),dimension(1:twotondim,1:ndim),save::xc
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_msun
  real(dp),dimension(1:twotondim),save::m8, mz8, n8, nph8 ! SNe
  real(dp),dimension(1:twotondim,1:3),save:: p8  ! SNe
  real(dp),dimension(1:nvector),save::nSNe, mSNe, mZSNe, nphSNe
  real(dp),dimension(1:nvector,1:3),save::pSNe
  real(dp),dimension(1:nvector,1:nchem),save::mchSNe
  real(dp),dimension(1:twotondim,1:nchem),save::mch8 ! SNe
  real(dp)::mejecta,mejecta_ch(1:nchem),Zejecta,M_SN_var
  real(dp),parameter::msun2g=2d33
  real(dp),parameter::myr2s=3.1536000d+13
  real(dp)::ttsta,ttend
  logical::ok,done_star
#ifdef RT
  real(dp),allocatable,dimension(:)::L_star
  real(dp)::Z_star,age,L_star_ion
  integer::iph,igroup
  Z_star=z_ave
  allocate(L_star(1:nSEDgroups))
#endif
 
  if(icount==2) return
  if(.not.hydro) return
  if(ndim.ne.3)  return
  if(numbtot(1,ilevel)==0)return
  if(nstar_tot==0)return

#ifndef WITHOUTMPI
  if(myid.eq.1) ttsta=MPI_WTIME(info)
#endif 
  nSNc=0; nsnIa_tot=0d0

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = scale_l**3*scale_d/msun2g

  ! Mesh spacing in that level
  call mesh_info (ilevel,skip_loc,scale,dx,dx_loc,vol_loc,xc) 

  ! To filter out old particles and compute individual time steps
  ! NB: the time step is always the coarser level time step, since no feedback for icount=2
  if (ilevel==levelmin)then
     dteff = dtnew(ilevel)
  else
     dteff = dtnew(ilevel-1)
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

  ncomm_SN=0  ! important to initialize; number of communications (not SNe)
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
    
           m8=0d0;mz8=0d0;p8=0d0;n8=0d0;nph8=0d0;mch8=0d0
    
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ok=.false.

              if(use_initial_mass) then
                 mass0 = mp0(ipart)*scale_msun
              else
                 mass0 = mp (ipart)*scale_msun
                 if(idp(ipart)>0) mass0 = mass0 / (1d0 - eta_sn)
              endif

              nsnIa_star=0d0
              if(snIa.and.tp(ipart).ne.0)then ! no dark matter
                 call get_number_of_snIa (tp(ipart),dteff,idp(ipart),mass0,nsnIa_star)
                 if(nsnIa_star>0) ok=.true.
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

                    !----------------------------------
                    ! For Type Ia explosions
                    !----------------------------------
                    ! total ejecta mass in code units
                    mejecta = mejecta_Ia / scale_msun * nsnIa_star
                    ! number of SNIa
                    n8(ind_son) = n8(ind_son) + nsnIa_star
                    ! mass return from SNe
                    m8 (ind_son)  = m8(ind_son) + mejecta
                    ! momentum from the original star, not the one generated by SNe
                    p8 (ind_son,1) = p8(ind_son,1) + mejecta*vp(ipart,1)
                    p8 (ind_son,2) = p8(ind_son,2) + mejecta*vp(ipart,2)
                    p8 (ind_son,3) = p8(ind_son,3) + mejecta*vp(ipart,3)
                    ! metal mass return from SNe including the newly synthesised one
                    if(metal)then
                       mz8(ind_son) = mz8(ind_son) + mejecta*1.0 ! 100% metals
                    endif
                    do ich=1,nchem
                       mch8(ind_son,ich) = mch8(ind_son,ich) + mejecta*Zejecta_chem_Ia(ich)
                    end do 

                    ! subtract the mass return
                    mp(ipart)=mp(ipart)-mejecta

                 endif
              endif
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
    
           do ind=1,twotondim
              if (abs(n8(ind))>0d0)then
                 ip=ip+1
                 ind_grid(ip)=igrid
                 ind_pos_cell(ip)=ind
     
                 ! collect information 
                 nSNe(ip)=n8(ind)
                 mSNe(ip)=m8(ind)
                 mZSNe(ip)=mz8(ind)
                 pSNe(ip,1)=p8(ind,1)
                 pSNe(ip,2)=p8(ind,2)
                 pSNe(ip,3)=p8(ind,3)
                 nphSNe(ip)=nph8(ind)  ! mechanical_geen
                 do ich=1,nchem
                    mchSNe(ip,ich)=mch8(ind,ich)
                 end do
 
                 nSNc=nSNc+1
                 nsnIa_tot=nsnIa_tot+nsnIa_star
                 if(ip==nvector)then
                    call mech_fine_snIa(ind_grid,ind_pos_cell,ip,ilevel,dteff,nSNe,mSNe,pSNe,mZSNe,nphSNe,mchSNe)
                    ip=0
                 endif 
              endif
           enddo
    
    
        end if
        igrid=next(igrid)   ! Go to next grid
     end do ! End loop over grids
    
     if (ip>0) then
        call mech_fine_snIa(ind_grid,ind_pos_cell,ip,ilevel,dteff,nSNe,mSNe,pSNe,mZSNe,nphSNe,mchSNe)
        ip=0
     endif

  end do ! End loop over cpus


#ifndef WITHOUTMPI
  nSNc_mpi=0; nsnIa_mpi=0d0
  ! Deal with the stars around the bounary of each cpu (need MPI)
  call mech_fine_snIa_mpi(ilevel)
  call MPI_ALLREDUCE(nSNc,nSNc_mpi,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(nsnIa_tot,nsnIa_mpi,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  nSNc = nSNc_mpi
  nsnIa_tot = nsnIa_mpi
  if(myid.eq.1.and.nSNc>0.and.log_mfb) then
     ttend=MPI_WTIME(info)
     write(*,*) '--------------------------------------'
     write(*,*) 'Time elapsed in mechanical_fine_Ia [sec]:', sngl(ttend-ttsta), nSNc, sngl(nsnIa_tot)
     write(*,*) '--------------------------------------'
  endif
#endif

#ifdef RT
   deallocate(L_star)
#endif

end subroutine mechanical_feedback_snIa_fine
!################################################################
!################################################################
!################################################################
subroutine mech_fine_snIa(ind_grid,ind_pos_cell,np,ilevel,dteff,nSN,mSN,pSN,mZSN,nphSN,mchSN)
  use amr_commons
  use pm_commons
  use hydro_commons
  use mechanical_commons
#ifdef RT
  use rt_parameters, ONLY:iIons,nIons
#endif
  implicit none
  integer::np,ilevel ! actually the number of cells
  integer,dimension(1:nvector)::ind_grid,ind_pos_cell
  real(dp),dimension(1:nvector)::nSN,mSN,mZSN,floadSN,nphSN
  real(dp),dimension(1:nvector)::mloadSN,mZloadSN,eloadSN
  real(dp),dimension(1:nvector,1:3)::pSN,ploadSN
  real(dp),dimension(1:nvector,1:nchem)::mchSN,mchloadSN
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine mechanical_feedback_fine 
  !-----------------------------------------------------------------------
  integer,dimension(nvector)::icellvec
  integer::i,j,nwco,nwco_here,idim,icell,igrid,ista,iend,ilevel2
  integer::ind_cell,ncell,irad,ii,ich
  real(dp)::d,u,v,w,e,z,eth,ekk,Tk,d0,u0,v0,w0,dteff
  real(dp)::dx,dx_loc,scale,vol_loc,nH_cen,fleftSN
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale_msun,msun2g=2d33
  real(dp)::skip_loc(1:3),Tk0,ekk0,eth0,etot0,T2min
  real(dp),dimension(1:twotondim,1:ndim),save::xc
  ! Grid based arrays
  real(dp),dimension(1:ndim,1:nvector),save::xc2
  real(dp),dimension(1:nvector,1:nSNnei), save::p_solid,ek_solid
  real(dp)::d_nei,Z_nei,Z_neisol,dm_ejecta,vol_nei
  real(dp)::mload,vload,Zload=0d0,f_esn2
  real(dp)::num_sn,nH_nei,Zdepen=1d0,f_w_cell,f_w_crit
  real(dp)::t_rad,r_rad,r_shell,m_cen,ekk_ej
  real(dp)::uavg,vavg,wavg,ul,vl,wl,ur,vr,wr
  real(dp)::d1,d2,d3,d4,d5,d6,dtot,pvar(1:ndim+3+512)
  real(dp)::vturb,vth,Mach,sig_s2,dratio,mload_cen
#ifdef RT
  real(dp),dimension(1:nIons),save::xion
#endif
  ! For stars affecting across the boundary of a cpu
  integer, dimension(1:nSNnei),save::icpuSNnei
  integer ,dimension(1:nvector,0:twondim):: ind_nbor
  logical,dimension(1:nvector,1:nSNnei),save ::snowplough
  real(dp),dimension(1:nvector)::rStrom ! in pc
  real(dp)::dx_loc_pc,psn_tr,chi_tr,psn_thor98,psn_geen15,fthor
  real(dp)::km2cm=1d5,M_SN_var
  ! chemical abundance
  real(dp),dimension(1:nchem)::chload,z_ch

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = scale_l**3*scale_d/msun2g

  ! Mesh variables
  call mesh_info (ilevel,skip_loc,scale,dx,dx_loc,vol_loc,xc)
  dx_loc_pc = dx_loc*scale_l/3.08d18

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
 
     num_sn    = nSN(i) ! doesn't have to be an integer
     M_SN_var  = mSN(i)*scale_msun / num_sn
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
        stop
     endif

     !==========================================
     ! estimate floadSN(i)
     !==========================================
     ! notice that f_LOAD / f_LEFT is a global parameter for SN ejecta themselves!!!
     if(loading_type.eq.1)then
        ! based on Federrath & Klessen (2012)
        ! find the mach number of this cell
        vth   = sqrt(gamma*1.38d-16*Tk/1.673d-24) ! cm/s
        ncell = 1
        icellvec(1) = icell
        call getnbor(icellvec,ind_nbor,ncell,ilevel)
        u0 = uold(icell        ,2) 
        v0 = uold(icell        ,3) 
        w0 = uold(icell        ,4) 
        ul = uold(ind_nbor(1,1),2)-u0 
        ur = uold(ind_nbor(1,2),2)-u0
        vl = uold(ind_nbor(1,3),3)-v0
        vr = uold(ind_nbor(1,4),3)-v0
        wl = uold(ind_nbor(1,5),4)-w0
        wr = uold(ind_nbor(1,6),4)-w0
        vturb = sqrt(ul**2 + ur**2 + vl**2 + vr**2 + wl**2 + wr**2)
        vturb = vturb*scale_v ! cm/s
    
        Mach  = vturb/vth
 
        ! get volume-filling density for beta=infty (non-MHD gas), b=0.4 (a stochastic mixture of forcing modes)
        sig_s2 = log(1d0 + (0.4*Mach)**2d0)
        ! s = -0.5*sig_s2;   s = ln(rho/rho0)
        dratio = max(exp(-0.5*sig_s2),0.01) ! physicall M~100 would be hard to achieve
        floadSN(i) = min(dratio,f_LOAD)

     else
        dratio     = 1d0
        floadSN(i) = f_LOAD
     endif

     if(log_mfb)then
398     format('MFB = ',f7.3,1x,f7.3,1x,f5.1,1x,f5.3,1x,f9.5,1x,f7.3)
        write(*,398) log10(d*scale_nH),log10(Tk),sngl(num_sn),floadSN(i),1./aexp-1,log10(dx_loc*scale_l/3.08d18)
     endif

     dm_ejecta = f_LOAD*mSN(i)/dble(nSNnei)  ! per solid angle
     mload     = f_LOAD*mSN(i) + uold(icell,1)*vol_loc*floadSN(i)  ! total SN ejecta + host cell
     if(metal) Zload = (f_LOAD*mZSN(i) + uold(icell,imetal)*vol_loc*floadSN(i))/mload
     do ich=1,nchem
        chload(ich) = (f_LOAD*mchSN(i,ich) + uold(icell,ichem+ich-1)*vol_loc*floadSN(i))/mload
     end do 
 
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
           nH_nei    = d_nei*scale_nH*dratio
           Z_neisol  = max(0.01, Z_nei/0.02)
           Zdepen    = Z_neisol**(expZ_SN*2d0) 

           ! transition mass loading factor (momentum conserving phase)
           ! psn_tr = sqrt(2*chi_tr*Nsn*Esn*Msn*fe)
           ! chi_tr = (1+f_w_crit)
           psn_thor98   = A_SN * num_sn**(expE_SN) * nH_nei**(expN_SN) * Z_neisol**(expZ_SN)  !km/s Msun
           psn_tr       = psn_thor98 

           chi_tr   = psn_tr**2d0 / (2d0 * num_sn**2d0 * (E_SNIa/msun2g/km2cm**2d0) * M_SN_var * f_esn)
           !          (Msun*km/s)^2                 (Msun *km2/s2)              (Msun)
           f_w_crit = max(chi_tr-1d0, 0d0)
 
           !f_w_crit = (A_SN/1d4)**2d0/(f_esn*M_SNIa)*num_sn**((expE_SN-1d0)*2d0)*nH_nei**(expN_SN*2d0)*Zdepen - 1d0
           !f_w_crit = max(0d0,f_w_crit)

           if(f_w_cell.ge.f_w_crit)then ! radiative phase
              ! ptot = sqrt(2*chi_tr*Mejtot*(fe*Esntot))
              ! vload = ptot/(chi*Mejtot) = sqrt(2*chi_tr*fe*Esntot/Mejtot)/chi = sqrt(2*chi_tr*fe*Esn/Mej)/chi
              vload = dsqrt(2d0*f_esn*E_SNIa*(1d0+f_w_crit)/(M_SN_var*msun2g))/scale_v/(1d0+f_w_cell)/f_LOAD
              snowplough(i,j)=.true.
           else ! adiabatic phase
              ! ptot = sqrt(2*chi*Mejtot*(fe*Esntot))
              ! vload = ptot/(chi*Mejtot) = sqrt(2*fe*Esntot/chi/Mejtot) = sqrt(2*fe*Esn/chi/Mej)
              f_esn2 = 1d0-(1d0-f_esn)*f_w_cell/f_w_crit
              vload = dsqrt(2d0*f_esn2*E_SNIa/(1d0+f_w_cell)/(M_SN_var*msun2g))/scale_v/f_LOAD
              snowplough(i,j)=.false.
           endif
           p_solid(i,j)=(1d0+f_w_cell)*dm_ejecta*vload
           ek_solid(i,j)=ek_solid(i,j)+p_solid(i,j)*(vload*f_LOAD)/2d0 !ek=(m*v)*v/2, not (d*v)*v/2

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
     e     = uold(icell,5)
#if NENER>0
     do irad=1,nener
        e = e - uold(icell,ndim+2+irad)
     enddo
#endif
     ekk   = 0.5*d*(u**2 + v**2 + w**2)
     eth   = e - ekk  ! thermal pressure 
#ifdef RT
     do ii=0,nIons-1
        xion(ii+1) = uold(icell,iIons+ii)/d
     end do
#endif

     mloadSN (i) = mSN (i)*f_LOAD + d*vol_loc*floadSN(i)
     if(metal)then
        z = uold(icell,imetal)/d
        mZloadSN(i) = mZSN(i)*f_LOAD + d*z*vol_loc*floadSN(i)
     endif

     do ich=1,nchem
        z_ch(ich) = uold(icell,ichem+ich-1)/d
        mchloadSN(i,ich) = mchSN(i,ich)*f_LOAD + d*z_ch(ich)*vol_loc*floadSN(i)
     end do 

     ! original momentum by star + gas entrained from the SN cell
     ploadSN(i,1) = pSN(i,1)*f_LOAD + vol_loc*d*u*floadSN(i)
     ploadSN(i,2) = pSN(i,2)*f_LOAD + vol_loc*d*v*floadSN(i)
     ploadSN(i,3) = pSN(i,3)*f_LOAD + vol_loc*d*w*floadSN(i)

     ! update the hydro variable
     fleftSN = 1d0 - floadSN(i)
     uold(icell,1) = uold(icell,1)*fleftSN + mSN(i)  /vol_loc*f_LEFT 
     uold(icell,2) = uold(icell,2)*fleftSN + pSN(i,1)/vol_loc*f_LEFT  ! rho*v, not v
     uold(icell,3) = uold(icell,3)*fleftSN + pSN(i,2)/vol_loc*f_LEFT 
     uold(icell,4) = uold(icell,4)*fleftSN + pSN(i,3)/vol_loc*f_LEFT 
     if(metal)then
        uold(icell,imetal) = mZSN(i)/vol_loc*f_LEFT + d*z*fleftSN
     endif
     do ich=1,nchem
        uold(icell,ichem+ich-1) = mchSN(i,ich)/vol_loc*f_LEFT + d*z_ch(ich)*fleftSN
     end do
#ifdef RT
     do ii=0,nIons-1
        uold(icell,iIons+ii) = xion(ii+1) * uold(icell,1)
     end do
#endif

     ! original kinetic energy of the gas entrained
     eloadSN(i) = ekk*vol_loc*floadSN(i) 

     ! original thermal energy of the gas entrained (including the non-thermal part)
     eloadSN(i) = eloadSN(i) + eth*vol_loc*floadSN(i)

     ! reduce total energy as we are distributing it to the neighbours
     uold(icell,5) = uold(icell,5)*fleftSN 

     ! add the contribution from the original kinetic energy of SN particle
     d = mSN(i)/vol_loc
     u = pSN(i,1)/mSN(i)
     v = pSN(i,2)/mSN(i)
     w = pSN(i,3)/mSN(i)
     uold(icell,5) = uold(icell,5) + 0.5d0*d*(u**2 + v**2 + w**2)*f_LEFT
    
     ! add the contribution from the original kinetic energy of SN to outflow
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
              do ich=1,nchem
                 pvar(ichem+ich-1) = unew(icell, ichem+ich-1)
              end do 
#if NENER>0
              do irad=1,nener
                 pvar(ndim+2+irad) = unew(icell,ndim+2+irad)
              enddo
#endif
           else
              pvar(1:5) = uold(icell,1:5)
              if(metal) pvar(imetal) = uold(icell,imetal)
              do ich=1,nchem
                 pvar(ichem+ich-1) = uold(icell, ichem+ich-1)
              end do 
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
           ! For stability
           Tk0 =eth0/d0*scale_T2*(gamma-1)
           T2min=T2_star*(d0*scale_nH/n_star)**g_star
           if(Tk0<T2min)then
              eth0=T2min*d0/scale_T2/(gamma-1)
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
           etot0  = eth0+ekk0+ek_solid(i,j)/vol_nei ! additional energy from SNe+entrained gas

           ! the minimum thermal energy input floor
           d   = pvar(1)
           u   = pvar(2)/d
           v   = pvar(3)/d
           w   = pvar(4)/d
           ekk = 0.5*d*(u**2 + v**2 + w**2)

           pvar(5) = max(etot0, ekk+eth0)

           Tk = (pvar(5)-ekk)/d*scale_T2*(gamma-1)*0.6
           if(Tk<0)then
              print *, 'TKERR: mech (post-call): Tk<0 =',sngl(Tk),sngl(d*scale_nH),sngl(Tk0),ilevel2
              stop
           endif 

#if NENER>0
           do irad=1,nener
              pvar(5) = pvar(5) + pvar(ndim+2+irad)
           end do
#endif

           if(metal)then
               pvar(imetal)=pvar(imetal)+mzloadSN(i)/dble(nSNnei)/vol_nei
           end if
           do ich=1,nchem
               pvar(ichem+ich-1)=pvar(ichem+ich-1)+mchloadSN(i,ich)/dble(nSNnei)/vol_nei
           end do

           ! update the hydro variable
           if(ilevel>ilevel2)then ! touching level-1 cells
#ifdef RT
              do ii=0,nIons-1 ! get the ionisation fraction
                 xion(ii+1)=unew(icell,iIons+ii)/unew(icell,1)
              end do
#endif
              unew(icell,1:5) = pvar(1:5)
              if(metal) unew(icell,imetal) = pvar(imetal)
              do ich=1,nchem
                 unew(icell,ichem+ich-1) = pvar(ichem+ich-1)
              end do
#ifdef RT
              do ii=0,nIons-1 ! put it back with the updated density
                 unew(icell,iIons+ii)=unew(icell,1)*xion(ii+1)
              end do
#endif
#if NENER>0
              do irad=1,nener
                 unew(icell,ndim+2+irad) = pvar(ndim+2+irad)
              end do
#endif

           else

#ifdef RT
              do ii=0,nIons-1 ! get the ionisation fraction
                 xion(ii+1)=uold(icell,iIons+ii)/uold(icell,1)
              end do
#endif
              uold(icell,1:5) = pvar(1:5)
              if(metal) uold(icell,imetal) = pvar(imetal)
              do ich=1,nchem
                 uold(icell,ichem+ich-1) = pvar(ichem+ich-1)
              end do
#ifdef RT
              do ii=0,nIons-1 ! put it back with the updated density
                 uold(icell,iIons+ii)=uold(icell,1)*xion(ii+1)
              end do
#endif
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
           nwco_here=nwco
           ! remove redundant cpu list for this SN cell
           call redundant_non_1d(icpuSNnei(1:nwco), nwco_here, nwco)
        endif
        ista=ncomm_SN+1
        iend=ista+nwco-1

        if(iend>ncomm_max)then
           write(*,*) 'Error: increase ncomm_max in mechanical_fine.f90', ncomm_max, iend
           call clean_stop
        endif
        iSN_comm (ista:iend)=icpuSNnei(1:nwco)
        nSN_comm (ista:iend )=nSN(i)
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
        do ich=1,nchem
            mchloadSN_comm(ista:iend,ich)=mchloadSN(i,ich)
        end do
        if(mechanical_geen.and.rt) rSt_comm (ista:iend)=rStrom(i)
        ncomm_SN=ncomm_SN+nwco
     endif
#endif

  end do ! loop over SN cell

end subroutine mech_fine_snIa
!################################################################
!################################################################
!################################################################
subroutine mech_fine_snIa_mpi(ilevel)
  use amr_commons
  use mechanical_commons
  use pm_commons
  use hydro_commons
#ifdef RT
  use rt_parameters, ONLY:iIons,nIons
#endif
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
  real(dp)::f_esn2,d_nei,Z_nei,Z_neisol,f_w_cell,f_w_crit,nH_nei,dratio
  real(dp)::num_sn,vload,Tk,Zdepen=1d0,vol_nei,dm_ejecta
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale_msun,msun2g=2d33, pvar(1:ndim+3+512),etot0
  real(dp)::skip_loc(1:3),d,u,v,w,ekk,eth,d0,u0,v0,w0,eth0,ekk0,Tk0,ekk_ej,T2min
  integer::igrid,icell,ilevel,ilevel2,irad,ii,ich
  real(dp),dimension(1:twotondim,1:ndim),save::xc
  logical,allocatable,dimension(:,:)::snowplough
!  logical,dimension(1:nvector,1:nSNnei),save ::snowplough
#ifdef RT
  real(dp),dimension(1:nIons),save::xion
#endif
  real(dp),dimension(1:nchem),save::chloadSN_i
  real(dp)::rSt_i,dx_loc_pc,psn_tr,chi_tr,psn_thor98,psn_geen15,fthor
  real(dp)::km2cm=1d5,M_SN_var

  if(ndim.ne.3) return

  !============================================================
  ! For MPI communication
  !============================================================
  ncpu_send=0;ncpu_recv=0

  ncomm_SN_cpu=0 
  ncomm_SN_mpi=0
  ncomm_SN_mpi(myid)=ncomm_SN
  ! compute the total number of communications needed
  call MPI_ALLREDUCE(ncomm_SN_mpi,ncomm_SN_cpu,ncpu,&
                   & MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nSN_tot = sum(ncomm_SN_cpu)
  if(nSN_tot==0) return

  allocate(icpuSN_comm    (1:nSN_tot,1:2))
  allocate(icpuSN_comm_mpi(1:nSN_tot,1:2))

  ! index for mpi variable
  if(myid==1)then
     isend_sta = 0 
  else
     isend_sta = sum(ncomm_SN_cpu(1:myid-1)) 
  endif

  icpuSN_comm=0
  do i=1,ncomm_SN_cpu(myid)
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

  ncell_send = ncomm_SN_cpu(myid)
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
              SNsend(11 ,ncc)=nSN_comm     (i)
              if(metal)SNsend(12,ncc)=mZloadSN_comm(i)
              if(mechanical_geen.and.rt)SNsend(13,ncc)=rSt_comm(i)
              do ich=1,nchem
                 SNsend(13+ich,ncc)=mchloadSN_comm(i,ich)
              end do
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
  dx_loc_pc = dx_loc*scale_l/3.08d18

  np = ncell_recv
  if(ncell_recv>0) then
     allocate(p_solid(1:np,1:nSNnei))
     allocate(ek_solid(1:np,1:nSNnei))
     allocate(snowplough(1:np,1:nSNnei))
     p_solid=0d0;ek_solid=0d0;snowplough=.false.
  endif

  ! Compute the momentum first before redistributing mass
  do i=1,np

     xSN_i(1:3) = SNrecv(1:3,i)
     mSN_i      = SNrecv(4,i)
     mloadSN_i  = SNrecv(5,i)
     fload_i    = SNrecv(9,i)
     dm_ejecta  = f_LOAD*mSN_i/dble(nSNnei)
     num_sn     = SNrecv(11,i)
     ek_solid(i,:) = SNrecv(10,i)/dble(nSNnei) ! kinetic energy of the gas mass entrained from the host cell + SN
     if(metal) ZloadSN_i = SNrecv(12,i)/SNrecv(5,i)
     if(mechanical_geen.and.rt) rSt_i = SNrecv(13,i) ! Stromgren sphere in pc 

     M_SN_var = mSN_i*scale_msun/num_sn

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
           if(loading_type.eq.1.and.fload_i<f_LOAD)then
              dratio = fload_i
           else
              dratio = 1d0
           endif
           f_w_cell = (mloadSN_i/dble(nSNnei) + d_nei*vol_loc/8d0)/dm_ejecta - 1d0
           nH_nei   = d_nei*scale_nH*dratio
           Z_neisol = max(0.01,Z_nei/0.02)
           Zdepen   = Z_neisol**(expZ_SN*2d0) !From Thornton+(98)

           ! transition mass loading factor (momentum conserving phase)
           ! psn_tr = sqrt(2*chi_tr*Nsn*Esn*Msn*fe)
           ! chi_tr = (1+f_w_crit)
           psn_thor98   = A_SN * num_sn**(expE_SN) * nH_nei**(expN_SN) * Z_neisol**(expZ_SN)  !km/s Msun
           psn_tr       = psn_thor98 

           chi_tr   = psn_tr**2d0 / (2d0 * num_sn**2d0 * (E_SNIa/msun2g/km2cm**2d0) * M_SN_var * f_esn)
           !          (Msun*km/s)^2                 (Msun *km2/s2)              (Msun)
           f_w_crit = max(chi_tr-1d0, 0d0)
 
           !f_w_crit = (A_SN/1d4)**2d0/(f_esn*M_SNIa)*num_sn**((expE_SN-1d0)*2d0)*nH_nei**(expN_SN*2d0)*Zdepen - 1d0
           !f_w_crit = max(0d0,f_w_crit)

           if(f_w_cell.ge.f_w_crit)then ! radiative phase
              vload = dsqrt(2d0*f_esn*E_SNIa*(1d0+f_w_crit)/(M_SN_var*msun2g))/(1d0+f_w_cell)/scale_v/f_LOAD
              snowplough(i,j)=.true.
           else ! adiabatic phase
              f_esn2 = 1d0-(1d0-f_esn)*f_w_cell/f_w_crit
              vload = dsqrt(2d0*f_esn2*E_SNIa/(1d0+f_w_cell)/(M_SN_var*msun2g))/scale_v/f_LOAD
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
     if(metal) ZloadSN_i = SNrecv(12,i)/mloadSN_i
     do ich=1,nchem
        chloadSN_i(ich) = SNrecv(13+ich,i)/mloadSN_i
     end do

     do j=1,nSNnei
 
        call get_icell_from_pos (xSN_i(1:3)+xSNnei(1:3,j)*dx, ilevel+1, igrid, icell, ilevel2)
        if(cpu_map(father(igrid))==myid)then

           vol_nei = vol_loc*(2d0**ndim)**(ilevel-ilevel2)
           if(ilevel>ilevel2)then ! touching level-1 cells
              pvar(1:5) = unew(icell,1:5)
              if(metal) pvar(imetal) = unew(icell,imetal)
              do ich=1,nchem
                 pvar(ichem+ich-1) = unew(icell,ichem+ich-1)
              end do
#if NENER>0 
              do irad=1,nener
                 pvar(ndim+2+irad) = unew(icell,ndim+2+irad)
              end do
#endif
           else
              pvar(1:5) = uold(icell,1:5)
              if(metal) pvar(imetal) = uold(icell,imetal)
              do ich=1,nchem
                 pvar(ichem+ich-1) = uold(icell,ichem+ich-1)
              end do
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
           Tk0 = eth0/d0*scale_T2*(gamma-1)
           T2min = T2_star*(d0*scale_nH/n_star)**g_star
           if(Tk0<T2min)then
              eth0 = T2min*d0/scale_T2/(gamma-1)
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
           etot0  = eth0+ekk0+ek_solid(i,j)/vol_nei  ! additional energy from SNe+entrained gas

           ! the minimum thermal energy input floor
           d   = pvar(1)
           u   = pvar(2)/d
           v   = pvar(3)/d
           w   = pvar(4)/d
           ekk = 0.5*d*(u**2 + v**2 + w**2)

           pvar(5) = max(etot0, ekk+eth0)

#if NENER>0
           do irad=1,nener
              pvar(5) = pvar(5) + pvar(ndim+2+irad)
           end do
#endif

           if(metal)then
               pvar(imetal)=pvar(imetal)+mloadSN_i/dble(nSNnei)*ZloadSN_i/vol_nei
           end if
           do ich=1,nchem
               pvar(ichem+ich-1)=pvar(ichem+ich-1)+mloadSN_i/dble(nSNnei)*chloadSN_i(ich)/vol_nei
           end do

           ! update the hydro variable
           if(ilevel>ilevel2)then ! touching level-1 cells
#ifdef RT
              do ii=0,nIons-1 ! get the ionisation fraction
                 xion(ii+1)=unew(icell,iIons+ii)/unew(icell,1)
              end do
#endif
              unew(icell,1:5) = pvar(1:5)
              if(metal) unew(icell,imetal) = pvar(imetal)
              do ich=1,nchem
                 unew(icell,ichem+ich-1) = pvar(ichem+ich-1)
              end do
#ifdef RT
              do ii=0,nIons-1 ! put it back with the updated density
                 unew(icell,iIons+ii)=unew(icell,1)*xion(ii+1)
              end do
#endif
#if NENER>0
              do irad=1,nener
                 unew(icell,ndim+2+irad) = pvar(ndim+2+irad)
              end do
#endif

           else
#ifdef RT
              do ii=0,nIons-1 ! get the ionisation fraction
                 xion(ii+1)=uold(icell,iIons+ii)/uold(icell,1)
              end do
#endif
              uold(icell,1:5) = pvar(1:5)
              if(metal) uold(icell,imetal) = pvar(imetal)
              do ich=1,nchem
                 uold(icell,ichem+ich-1) = pvar(ichem+ich-1)
              end do
#ifdef RT
              do ii=0,nIons-1 ! put it back with the updated density
                 uold(icell,iIons+ii)=uold(icell,1)*xion(ii+1)
              end do
#endif
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
  if(ncell_recv>0) deallocate(p_solid,ek_solid,snowplough)

  ncomm_SN=nSN_tot
#endif

end subroutine mech_fine_snIa_mpi
!################################################################
!################################################################
!################################################################
subroutine get_number_of_snIa (birth_time, dteff, id_star, mass0, nsnIa )
  use amr_commons
  implicit none
  integer ::id_star
  real(dp)::mass0 ! initial mass of a star particle in Msun
  real(dp)::nsnIa ! number of snIa
  real(dp)::birth_time,dteff
!-------------------------------------------------------------
!	Use the inverse method to generate random numbers for snIa
!-------------------------------------------------------------
  real(dp)::A_DTD,t_ini,t_fin,xdum,ydum,age1,age2
  integer ::localseed,i,nsnIa_tot
  real(dp),external::ran1


!  DTD = A_DTD* t^-1
!  1 = int A_DTD / t dt
!    = A_DTD*(log_e (t_f) - log_e (t_i))
!  A_DTD = 1d0/(alog(t_fin) - alog(t_ini))

!  n_sn(<t) = A_DTD*(log (t) - log(t_ini)) ; 0 < n_sn < 1
!  log(t) = n_sn/A_DTD + log(t_ini)
!  t = exp(n_sn/A_DTD + log(t_ini))


  ! get stellar age
  call getStarAgeGyr(birth_time+dteff, age1)
  call getStarAgeGyr(birth_time      , age2)

  ! convert Gyr -> Myr
  age1 = age1*1d3
  age2 = age2*1d3

  t_ini = 50d6
  t_fin = 1.37d10
  A_DTD = 1d0 / (log(t_fin) - log(t_ini))

  nsnIa = 0d0
  if(age2*1d6.lt.t_ini) then
     return
  endif

  nsnIa_tot = NINT(mass0 * A_snIa)
  localseed = -ABS(id_star)
  do i=1,nsnIa_tot
     xdum = ran1(localseed)
     ydum = exp(xdum / A_DTD + log(t_ini))/1d6
     if(ydum.ge.age1.and.ydum.le.age2) nsnIa = nsnIa + 1
  end do 

end subroutine get_number_of_snIa
