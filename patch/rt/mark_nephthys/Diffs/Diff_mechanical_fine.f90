1,4c1,4
< ! This replaced mechanical_fine_1510.f90 on 14 Nov 2015
< ! Minor fixed to thermal injection and accounting for NENER and 
< ! ionisation fractions
< ! Taysun added chem (Jun 2016)
---
> ! This replaced mechanical_fine_1510.f90 on 14 Nov 2015.
> ! Minor fixes to thermal injection and accounting for NENER and
> ! ionisation fractions.
> 
13,16d12
< #ifdef RT
<   use rt_parameters,only:group_egy
<   use SED_module,only:nSEDgroups,inp_SED_table
< #endif
22,23c18,20
<   ! This routine computes the energy liberated from supernova II ,
<   ! and inject momentum and energy to the surroundings of the young stars.
---
>   ! This routine computes the energy liberated from stellar winds
>   ! based on the input library, and inject thermal or kinetic energy to 
>   ! the surrounding of the young stars.
26,33d22
<   ! m8,mz8,nph8: temporary variable necessary for an oct to add up 
<   !              the mass, metal, etc. on a cell by cell basis
<   ! mejecta: mass of the ejecta from SNe
<   ! Zejecta: metallicity of the ejecta (not yield)
<   ! mZSNe : total metal mass from SNe in each cell
<   ! mchSNe: total mass of each chemical element in each cell
<   ! nphSNe: total production rate of ionising radiation in each cell.
<   !         This is necessary to estimate the Stromgren sphere
36,38c25,26
<   integer::npart1,npart2,icpu,icount,idim,ip,ich
<   integer::ind,ind_son,ind_cell,ilevel,iskip,info
<   integer::nSNc,nSNc_mpi
---
>   integer::npart1,npart2,icpu,icount,idim,ip
>   integer::ind,ind_son,ind_cell,ilevel,iskip,info,nSNc,nSNc_mpi,nsn_star
40d27
<   real(dp)::nsnII_star,mass0,nsnII_tot,nsnII_mpi
44,51c31,37
<   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_msun
<   real(dp),dimension(1:twotondim),save::m8, mz8, n8, nph8 ! SNe
<   real(dp),dimension(1:twotondim,1:3),save:: p8  ! SNe
<   real(dp),dimension(1:nvector),save::nSNe, mSNe, mZSNe, nphSNe
<   real(dp),dimension(1:nvector,1:3),save::pSNe
<   real(dp),dimension(1:nvector,1:nchem),save::mchSNe
<   real(dp),dimension(1:twotondim,1:nchem),save::mch8 ! SNe
<   real(dp)::mejecta,mejecta_ch(1:nchem),Zejecta,mfrac_snII,M_SN_var,snII_freq
---
>   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
>   real(dp)::scale_msun
>   real(dp),dimension(1:twotondim),save::mw8, mzw8 ! SNe
>   real(dp),dimension(1:twotondim,1:3),save:: pw8  ! SNe
>   real(dp),dimension(1:nvector),save::mSN, mZSN
>   real(dp),dimension(1:nvector,1:3),save::pSN
>   real(dp)::mejecta
56,62d41
< #ifdef RT
<   real(dp),allocatable,dimension(:)::L_star
<   real(dp)::Z_star,age,L_star_ion
<   integer::iph,igroup
<   Z_star=z_ave
<   allocate(L_star(1:nSEDgroups))
< #endif
73c52
<   nSNc=0; nsnII_tot=0d0
---
>   nSNc=0
85c64
<      dteff = dtnew(ilevel)
---
>      dteff = dtold(ilevel)
87c66
<      dteff = dtnew(ilevel-1)
---
>      dteff = dtold(ilevel-1)
100c79
<   ncomm_SN=0  ! important to initialize; number of communications (not SNe)
---
>   nSN_comm=0  ! important to initialize; number of communications (not SNe)
107,109d85
<   ! Type II Supernova frequency per Msun
<   snII_freq = eta_sn / M_SNII
< 
128c104
<            m8=0d0;mz8=0d0;p8=0d0;n8=0d0;nph8=0d0;mch8=0d0
---
>            mw8=0d0;mzw8=0d0;pw8=0d0
135,163c111,121
< 
<               if(use_initial_mass) then
<                  mass0 = mp0(ipart)*scale_msun
<               else
<                  mass0 = mp (ipart)*scale_msun
<                  if(idp(ipart)>0) mass0 = mass0 / (1d0 - eta_sn)
<               endif
< 
< #ifdef POP3
<               if(pop3 .and. (zp(ipart).lt.Zcrit_pop3) ) then
<                  ok=.false. ! this will be done elsewhere
<                  done_star=.false.
<               else
< #endif
<                  nsnII_star=0d0
<                  if(sn2_real_delay)then
<                     ! if tp is younger than t_delay
<                     if (idp(ipart).le.0.and.tp(ipart).ge.tyoung) then
<                        call get_number_of_sn2  (tp(ipart), dteff, zp(ipart), idp(ipart),&
<                                 & mass0, nsnII_star, done_star)
<                        if(nsnII_star>0)ok=.true.
<                     endif
<                  else ! single SN event
<                     ! if tp is older than t_delay
<                     if (idp(ipart).le.0.and.tp(ipart).le.tyoung)then
<                        ok=.true. 
<                        ! number of sn is not necessarily an integer
<                        nsnII_star = mass0*eta_sn/M_SNII
<                     endif
---
>               if(sn2_real_delay)then
>                  ! if tp is younger than t_delay
>                  if (idp(ipart).le.0.and.tp(ipart).ge.tyoung) then
>                     call get_number_of_sn2  (tp(ipart), zp(ipart), idp(ipart),&
>                              & mp0(ipart)*scale_msun,mp(ipart)*scale_msun,nsn_star,done_star)
>                     if(nsn_star>0)ok=.true.
>                  endif
>               else ! single SN event
>                  ! if tp is older than t_delay
>                  if (idp(ipart).le.0.and.tp(ipart).le.tyoung)then
>                     ok=.true.
165d122
< #ifdef POP3
167d123
< #endif
179,191c135,138
< 
<                     !----------------------------------
<                     ! For Type II explosions
<                     !----------------------------------
<                     M_SN_var = M_SNII
<                     if(metal) then
<                        if(variable_yield_SNII)then
<                           call SNII_yield (zp(ipart), mfrac_snII, Zejecta, Zejecta_chem_II)
<                           ! Adjust M_SNII mass not to double-count the mass loss from massive stars
<                           M_SN_var = mass_loss_boost * mfrac_snII / snII_freq ! ex) 0.1 / 0.01 = 10 Msun
<                        else
<                           Zejecta = zp(ipart)+(1d0-zp(ipart))*yield
<                        endif
---
>                     if(sn2_real_delay)then
>                        mejecta = M_SNII/scale_msun*nsn_star
>                     else
>                        mejecta = mp(ipart)*eta_sn
192a140,145
>                     !mp_sun = mp(ipart)*scale_msun*eta_sn
>                     mw8 (ind_son)  =mw8(ind_son)+mejecta
>                     pw8 (ind_son,1)=pw8(ind_son,1)+mejecta*vp(ipart,1)
>                     pw8 (ind_son,2)=pw8(ind_son,2)+mejecta*vp(ipart,2)
>                     pw8 (ind_son,3)=pw8(ind_son,3)+mejecta*vp(ipart,3)
>                     if(metal) mzw8(ind_son) = mzw8(ind_son) + mejecta*(zp(ipart)+(1d0-zp(ipart))*yield)
194,207c147,150
<                     ! total ejecta mass in code units
<                     mejecta = M_SN_var/scale_msun*nsnII_star
< 
<                     ! number of SNII
<                     n8(ind_son) = n8(ind_son) + nsnII_star
<                     ! mass return from SNe
<                     m8 (ind_son)  = m8(ind_son) + mejecta
<                     ! momentum from the original star, not the one generated by SNe
<                     p8 (ind_son,1) = p8(ind_son,1) + mejecta*vp(ipart,1)
<                     p8 (ind_son,2) = p8(ind_son,2) + mejecta*vp(ipart,2)
<                     p8 (ind_son,3) = p8(ind_son,3) + mejecta*vp(ipart,3)
<                     ! metal mass return from SNe including the newly synthesised one
<                     if(metal)then
<                        mz8(ind_son) = mz8(ind_son) + mejecta*Zejecta
---
>                     ! Record-keeping of local density at SN events (joki):
>                     if(write_stellar_densities) then
>                        st_n_SN(ipart) = uold(ind_cell,1)                                     
>                        st_e_SN(ipart) = eta_sn/(10.*2d33) * mp(ipart) * scale_d * scale_l**3 
209,213d151
<                     do ich=1,nchem
<                        mch8(ind_son,ich) = mch8(ind_son,ich) + mejecta*Zejecta_chem_II(ich)
<                     end do 
< 
<                     ! subtract the mass return
216d153
<                     ! mark if we are done with this particle
222,247d158
< 
<                     ! Record-keeping of local density at SN events (joki):
<                     if(write_stellar_densities) then
<                        st_n_SN(ipart) = uold(ind_cell,1)                                     
<                        st_e_SN(ipart) = nsnII_star  ! eta_sn/(10.*2d33) * mp(ipart) * scale_d * scale_l**3 ! fixed by Taysun 
<                     endif
<  
< 
< #ifdef RT
<                     ! Enhanced momentum due to pre-processing of the ISM due to radiation
<                     if(rt.and.mechanical_geen) then
<                        ! Let's count the total number of ionising photons per sec 
<                        call getAgeGyr(tp(ipart), age)
<                        if(metal) Z_star=zp(ipart)
<                        Z_star=max(Z_star,10.d-5)
< 
<                        ! compute the number of ionising photons from SED
<                        call inp_SED_table(age, Z_star, 1, .false., L_star) ! L_star = [# s-1 Msun-1]
<                        L_star_ion = 0d0
<                        do igroup=1,nSEDgroups
<                           if(group_egy(igroup).ge.13.6) L_star_ion = L_star_ion + L_star(igroup)
<                        end do
<                        nph8 (ind_son)=nph8(ind_son) + mass0*L_star_ion ! [# s-1]
<                     endif
< #endif
< 
255c166,167
<               if (abs(n8(ind))>0d0)then
---
>               if (abs(mw8(ind))>0d0)then
>                  nSNc=nSNc+1
261,275c173,178
<                  nSNe(ip)=n8(ind)
<                  mSNe(ip)=m8(ind)
<                  mZSNe(ip)=mz8(ind)
<                  pSNe(ip,1)=p8(ind,1)
<                  pSNe(ip,2)=p8(ind,2)
<                  pSNe(ip,3)=p8(ind,3)
<                  nphSNe(ip)=nph8(ind)  ! mechanical_geen
<                  do ich=1,nchem
<                     mchSNe(ip,ich)=mch8(ind,ich)
<                  end do
<    
<                  ! statistics
<                  nSNc=nSNc+1
<                  nsnII_tot = nsnII_tot + nsnII_star
<  
---
>                  mSN(ip)=mw8(ind)
>                  mZSN(ip)=mzw8(ind)
>                  pSN(ip,1)=pw8(ind,1)
>                  pSN(ip,2)=pw8(ind,2)
>                  pSN(ip,3)=pw8(ind,3)
>     
277c180
<                     call mech_fine(ind_grid,ind_pos_cell,ip,ilevel,dteff,nSNe,mSNe,pSNe,mZSNe,nphSNe,mchSNe)
---
>                     call mech_fine(ind_grid,ind_pos_cell,ip,ilevel,mSN,pSN,mZSN,dteff)
289c192
<         call mech_fine(ind_grid,ind_pos_cell,ip,ilevel,dteff,nSNe,mSNe,pSNe,mZSNe,nphSNe,mchSNe)
---
>         call mech_fine(ind_grid,ind_pos_cell,ip,ilevel,mSN,pSN,mZSN,dteff)
297c200
<   nSNc_mpi=0; nsnII_mpi=0d0
---
>   nSNc_mpi=0
301d203
<   call MPI_ALLREDUCE(nsnII_tot,nsnII_mpi,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
303d204
<   nsnII_tot = nSNII_mpi
307c208
<      write(*,*) 'Time elapsed in mechanical_fine [sec]:', sngl(ttend-ttsta), nSNc, sngl(nsnII_tot) 
---
>      write(*,*) 'Time elapsed in mechanical_fine [sec]:', sngl(ttend-ttsta), nSNc, nSN_comm 
316c217
< subroutine mech_fine(ind_grid,ind_pos_cell,np,ilevel,dteff,nSN,mSN,pSN,mZSN,nphSN,mchSN)
---
> subroutine mech_fine(ind_grid,ind_pos_cell,np,ilevel,mSN,pSN,mZSN,dteff)
327c228
<   real(dp),dimension(1:nvector)::nSN,mSN,mZSN,floadSN,nphSN
---
>   real(dp),dimension(1:nvector)::mSN,mZSN,floadSN
330d230
<   real(dp),dimension(1:nvector,1:nchem)::mchSN,mchloadSN
334,335c234
<   integer::i,j,nwco,nwco_here,idim,icell,igrid,ista,iend,ilevel2
<   integer::ind_cell,ncell,irad,ii,ich
---
>   integer::i,j,nwco,nwco_here,idim,icell,igrid,ista,iend,ilevel2,ind_cell,ncell,irad,ii
337c236
<   real(dp)::dx,dx_loc,scale,vol_loc,nH_cen,fleftSN
---
>   real(dp)::dx,dx_loc,scale,vol_loc,nH_cen,tsim_yr,fleftSN
340c239
<   real(dp)::skip_loc(1:3),Tk0,ekk0,eth0,etot0,T2min
---
>   real(dp)::skip_loc(1:3),Tk0,ekk0,eth0,etot0
345,346c244,245
<   real(dp)::d_nei,Z_nei,Z_neisol,dm_ejecta,vol_nei
<   real(dp)::mload,vload,Zload=0d0,f_esn2
---
>   real(dp)::d_nei,Z_nei,dm_ejecta,vol_nei
>   real(dp)::mload,vload,Zload=0d0,f_esn2,eturb
359,363d257
<   real(dp),dimension(1:nvector)::rStrom ! in pc
<   real(dp)::dx_loc_pc,psn_tr,chi_tr,psn_thor98,psn_geen15,fthor
<   real(dp)::km2cm=1d5,M_SN_var,boost_geen_ad=0d0,p_hydro,vload_rad,f_wrt_snow
<   ! chemical abundance
<   real(dp),dimension(1:nchem)::chload,z_ch
367a262
>   tsim_yr    = dteff*scale_t/365d0/3600d0/24d0
371d265
<   dx_loc_pc = dx_loc*scale_l/3.08d18
407,408c301
<      num_sn    = nSN(i) ! doesn't have to be an integer
<      M_SN_var = mSN(i)*scale_msun / num_sn
---
>      num_sn    = mSN(i)/(M_SNII/scale_msun) ! doesn't have to be an integer
427d319
<         stop
434,469c326
<      if(loading_type.eq.1)then
<         ! based on Federrath & Klessen (2012)
<         ! find the mach number of this cell
<         vth   = sqrt(gamma*1.38d-16*Tk/1.673d-24) ! cm/s
<         ncell = 1
<         call getnbor(icell,ind_nbor,ncell,ilevel)
<         u0 = uold(icell        ,2) 
<         v0 = uold(icell        ,3) 
<         w0 = uold(icell        ,4) 
<         ul = uold(ind_nbor(1,1),2)-u0 
<         ur = uold(ind_nbor(1,2),2)-u0
<         vl = uold(ind_nbor(1,3),3)-v0
<         vr = uold(ind_nbor(1,4),3)-v0
<         wl = uold(ind_nbor(1,5),4)-w0
<         wr = uold(ind_nbor(1,6),4)-w0
<         vturb = sqrt(ul**2 + ur**2 + vl**2 + vr**2 + wl**2 + wr**2)
<         vturb = vturb*scale_v ! cm/s
<     
<         Mach  = vturb/vth
<  
<         ! get volume-filling density for beta=infty (non-MHD gas), b=0.4 (a stochastic mixture of forcing modes)
<         sig_s2 = log(1d0 + (0.4*Mach)**2d0)
<         ! s = -0.5*sig_s2;   s = ln(rho/rho0)
<         dratio = max(exp(-0.5*sig_s2),0.01) ! physicall M~100 would be hard to achieve
<         floadSN(i) = min(dratio,f_LOAD)
< 
<      else
<         dratio     = 1d0
<         floadSN(i) = f_LOAD
<      endif
< 
<      !==========================================
<      ! estimate Stromgren sphere (relevant to RHD simulations only)
<      ! (mechanical_geen=.true)
<      !==========================================
<      if(mechanical_geen.and.rt) rStrom(i) = (3d0*nphSN(i)/4./3.141592/2.6d-13/nH_cen**2d0)**(1d0/3d0)/3.08d18 ! [pc] 
---
>      floadSN(i) = f_LOAD 
479,482c336,337
<      do ich=1,nchem
<         chload(ich) = (f_LOAD*mchSN(i,ich) + uold(icell,ichem+ich-1)*vol_loc*floadSN(i))/mload
<      end do 
<  
---
> 
>     
495d349
< 
497,499c351,352
<            nH_nei    = d_nei*scale_nH*dratio
<            Z_neisol  = max(0.01, Z_nei/0.02)
<            Zdepen    = Z_neisol**(expZ_SN*2d0) 
---
>            nH_nei    = d_nei*scale_nH
>            Zdepen    = (max(0.01,Z_nei/0.02))**(expZ_SN*2d0) 
502,516c355,356
<            ! psn_tr = sqrt(2*chi_tr*Nsn*Esn*Msn*fe)
<            ! chi_tr = (1+f_w_crit)
<            psn_thor98   = A_SN * num_sn**(expE_SN) * nH_nei**(expN_SN) * Z_neisol**(expZ_SN)  !km/s Msun
<            psn_tr       = psn_thor98 
<            if(mechanical_geen)then 
<               ! For snowplough phase, psn_tr will do the job
<               psn_geen15 = A_SN_Geen * num_sn**(expE_SN)* Z_neisol**(expZ_SN)  !km/s Msun
< 
<               if(rt)then
<                  fthor   = exp(-dx_loc_pc/rStrom(i))
<                  psn_tr  = psn_thor98*fthor + psn_geen15*(1d0-fthor)
<               else
<                  psn_tr  = psn_geen15
<               endif
<               psn_tr = max(psn_tr, psn_thor98)
---
>            f_w_crit = (A_SN/1d4)**2d0/(f_esn*M_SNII)*num_sn**((expE_SN-1d0)*2d0)*nH_nei**(expN_SN*2d0)*Zdepen - 1d0
>            f_w_crit = max(0d0,f_w_crit)
518,532d357
<               ! For adiabatic phase
<               ! psn_tr =  A_SN * (E51 * boost_geen)**expE_SN_boost * nH_nei**(expN_SN_boost) * Z_neisol**(expZ_SN)
<               !        =  p_hydro * boost_geen**expE_SN_boost
<               p_hydro = A_SN * num_sn**(expN_SN_boost) * nH_nei**(expN_SN_boost) * Z_neisol**(expZ_SN)
<               boost_geen_ad = (psn_tr / p_hydro)**(1d0/expE_SN_boost)
<               boost_geen_ad = max(boost_geen_ad-1d0,0.0)
<            endif
< 
<            chi_tr   = psn_tr**2d0 / (2d0 * num_sn * (E_SNII/msun2g/km2cm**2d0) * M_SN_var * f_ESN)
<            !          (Msun*km/s)^2                 (Msun *km2/s2)              (Msun)
<            f_w_crit = max(chi_tr-1d0, 0d0)
< 
<            !f_w_crit = (A_SN/1d4)**2d0/(f_ESN*M_SNII)*num_sn**((expE_SN-1d0)*2d0)*nH_nei**(expN_SN*2d0)*Zdepen - 1d0
<            !f_w_crit = max(0d0,f_w_crit)
<            vload_rad = dsqrt(2d0*f_ESN*E_SNII*(1d0+f_w_crit)/(M_SN_var*msun2g))/scale_v/(1d0+f_w_cell)/f_LOAD/f_PCAN
534,536c359
<               ! ptot = sqrt(2*chi_tr*Mejtot*(fe*Esntot))
<               ! vload = ptot/(chi*Mejtot) = sqrt(2*chi_tr*fe*Esntot/Mejtot)/chi = sqrt(2*chi_tr*fe*Esn/Mej)/chi
<               vload = vload_rad
---
>               vload = dsqrt(2d0*f_esn*E_SNII*(1d0+f_w_crit)/(M_SNII*msun2g))/scale_v/(1d0+f_w_cell)/f_LOAD
539,550c362,363
<               ! ptot = sqrt(2*chi*Mejtot*(fe*Esntot))
<               ! vload = ptot/(chi*Mejtot) = sqrt(2*fe*Esntot/chi/Mejtot) = sqrt(2*fe*Esn/chi/Mej)
<               f_esn2 = 1d0-(1d0-f_ESN)*f_w_cell/f_w_crit ! to smoothly correct the adibatic to the radiative phase
<               vload = dsqrt(2d0*f_esn2*E_SNII/(1d0+f_w_cell)/(M_SN_var*msun2g))/scale_v/f_LOAD
<               if(mechanical_geen) then
<                  !f_wrt_snow = (f_esn2 - f_ESN)/(1d0-f_ESN)
<                  f_wrt_snow = 2d0-2d0/(1d0+exp(-f_w_cell/f_w_crit/0.3)) ! 0.3 is obtained by calibrating
<                  vload = vload * dsqrt(1d0 + boost_geen_ad*f_wrt_snow)
<                  ! NB. this sometimes give too much momentum because expE_SN_boost != expE_SN. A limiter is needed
<               endif
<               ! safety device: limit the maximum velocity so that it does not exceed p_{SN,final}
<               if(vload>vload_rad) vload = vload_rad
---
>               f_esn2 = 1d0-(1d0-f_esn)*f_w_cell/f_w_crit
>               vload = dsqrt(2d0*f_esn2*E_SNII/(1d0+f_w_cell)/(M_SNII*msun2g))/scale_v/f_LOAD
597,601d409
<      do ich=1,nchem
<         z_ch(ich) = uold(icell,ichem+ich-1)/d
<         mchloadSN(i,ich) = mchSN(i,ich)*f_LOAD + d*z_ch(ich)*vol_loc*floadSN(i)
<      end do 
< 
609,612c417,420
<      uold(icell,1) = uold(icell,1)*fleftSN + mSN(i)  /vol_loc*f_LEFT 
<      uold(icell,2) = uold(icell,2)*fleftSN + pSN(i,1)/vol_loc*f_LEFT  ! rho*v, not v
<      uold(icell,3) = uold(icell,3)*fleftSN + pSN(i,2)/vol_loc*f_LEFT 
<      uold(icell,4) = uold(icell,4)*fleftSN + pSN(i,3)/vol_loc*f_LEFT 
---
>      uold(icell,1) = mSN(i)  /vol_loc*f_LEFT + d*fleftSN
>      uold(icell,2) = pSN(i,1)/vol_loc*f_LEFT + d*u*fleftSN ! make sure this is rho*v, not v
>      uold(icell,3) = pSN(i,2)/vol_loc*f_LEFT + d*v*fleftSN
>      uold(icell,4) = pSN(i,3)/vol_loc*f_LEFT + d*w*fleftSN
616,618d423
<      do ich=1,nchem
<         uold(icell,ichem+ich-1) = mchSN(i,ich)/vol_loc*f_LEFT + d*z_ch(ich)*fleftSN
<      end do
632c437
<      uold(icell,5) = uold(icell,5)*fleftSN 
---
>      uold(icell,5) = uold(icell,5) - (ekk+eth)*floadSN(i)
634c439
<      ! add the contribution from the original kinetic energy of SN particle
---
>      ! add the contribution from the original kinetic energy of SN
640c445
<     
---
>  
670,672d474
<               do ich=1,nchem
<                  pvar(ichem+ich-1) = unew(icell, ichem+ich-1)
<               end do 
681,683d482
<               do ich=1,nchem
<                  pvar(ichem+ich-1) = uold(icell, ichem+ich-1)
<               end do 
696a496
>            Tk0 =eth0/d0*scale_T2*(gamma-1)*0.6 !not used!
702,708d501
<            ! For stability
<            Tk0 =eth0/d0*scale_T2*(gamma-1)
<            T2min=T2_star*(d0*scale_nH/n_star)**g_star
<            if(Tk0<T2min)then
<               eth0=T2min*d0/scale_T2/(gamma-1)
<            endif 
< 
729a523,524
>            eturb = pvar(5) - ekk - eth0
> 
731,732c526,527
<            if(Tk<0)then
<               print *, 'TKERR: mech (post-call): Tk<0 =',sngl(Tk),sngl(d*scale_nH),sngl(Tk0),ilevel2
---
>            if(Tk<0.)then
>               print *, 'TKERR: mech (post-call): Tk<0.01 =', Tk
745,747d539
<            do ich=1,nchem
<                pvar(ichem+ich-1)=pvar(ichem+ich-1)+mchloadSN(i,ich)/dble(nSNnei)/vol_nei
<            end do
758,760d549
<               do ich=1,nchem
<                  unew(icell,ichem+ich-1) = pvar(ichem+ich-1)
<               end do
781,783d569
<               do ich=1,nchem
<                  uold(icell,ichem+ich-1) = pvar(ichem+ich-1)
<               end do
804d589
<            nwco_here=nwco
806c591,596
<            call redundant_non_1d(icpuSNnei(1:nwco), nwco_here, nwco)
---
>            !call redundant_non_1d(icpuSNnei(1:nwco), nwco, nwco)
>            if(nwco>1)then
>               nwco_here=nwco
>               ! remove redundant cpu list for this SN cell
>               call redundant_non_1d(icpuSNnei(1:nwco_here), nwco_here, nwco)
>            endif
808c598
<         ista=ncomm_SN+1
---
>         ista=nSN_comm+1
816c606
<         nSN_comm (ista:iend )=nSN(i)
---
>         !lSN_comm (ista:iend)=ilevel
828,832c618,619
<         do ich=1,nchem
<             mchloadSN_comm(ista:iend,ich)=mchloadSN(i,ich)
<         end do
<         if(mechanical_geen.and.rt) rSt_comm (ista:iend)=rStrom(i)
<         ncomm_SN=ncomm_SN+nwco
---
>         nSN_comm=nSN_comm+nwco
> 
863c650
<   real(dp)::f_esn2,d_nei,Z_nei,Z_neisol,f_w_cell,f_w_crit,nH_nei,dratio
---
>   real(dp)::f_esn2,d_nei,Z_nei,f_w_cell,f_w_crit,nH_nei
866,868c653,655
<   real(dp)::scale_msun,msun2g=2d33, pvar(1:ndim+3+512),etot0
<   real(dp)::skip_loc(1:3),d,u,v,w,ekk,eth,d0,u0,v0,w0,eth0,ekk0,Tk0,ekk_ej,T2min
<   integer::igrid,icell,ilevel,ilevel2,irad,ii,ich
---
>   real(dp)::scale_msun,msun2g=2d33, pvar(1:ndim+3+512),eturb,etot0
>   real(dp)::skip_loc(1:3),d,u,v,w,ekk,eth,d0,u0,v0,w0,eth0,ekk0,Tk0,ekk_ej
>   integer::igrid,icell,ilevel,ilevel2,irad,ii
874,876d660
<   real(dp),dimension(1:nchem),save::chloadSN_i
<   real(dp)::rSt_i,dx_loc_pc,psn_tr,chi_tr,psn_thor98,psn_geen15,fthor
<   real(dp)::km2cm=1d5,M_SN_var,boost_geen_ad=0d0,p_hydro,vload_rad,f_wrt_snow
887,889c671,673
<   ncomm_SN_cpu=0 
<   ncomm_SN_mpi=0
<   ncomm_SN_mpi(myid)=ncomm_SN
---
>   nSN_comm_cpu=0 
>   nSN_comm_mpi=0
>   nSN_comm_mpi(myid)=nSN_comm
891c675
<   call MPI_ALLREDUCE(ncomm_SN_mpi,ncomm_SN_cpu,ncpu,&
---
>   call MPI_ALLREDUCE(nSN_comm_mpi,nSN_comm_cpu,ncpu,&
893c677
<   nSN_tot = sum(ncomm_SN_cpu)
---
>   nSN_tot = sum(nSN_comm_cpu)
903c687
<      isend_sta = sum(ncomm_SN_cpu(1:myid-1)) 
---
>      isend_sta = sum(nSN_comm_cpu(1:myid-1)) 
907c691
<   do i=1,ncomm_SN_cpu(myid)
---
>   do i=1,nSN_comm_cpu(myid)
920c704
<   ncell_send = ncomm_SN_cpu(myid)
---
>   ncell_send = nSN_comm_cpu(myid)
979,984c763
<               SNsend(11 ,ncc)=nSN_comm     (i)
<               if(metal)SNsend(12,ncc)=mZloadSN_comm(i)
<               if(mechanical_geen.and.rt)SNsend(13,ncc)=rSt_comm(i)
<               do ich=1,nchem
<                  SNsend(13+ich,ncc)=mchloadSN_comm(i,ich)
<               end do
---
>               if(metal)SNsend(11,ncc)=mZloadSN_comm(i)
1010a790
> ! NB. Some mpi routines do not like the receiving buffer
1013a794,795
>         
> 
1030c812
<   
---
> 
1033d814
<   dx_loc_pc = dx_loc*scale_l/3.08d18
1050c831,832
<      num_sn     = SNrecv(11,i)
---
>      num_sn     = mSN_i/(M_SNII/scale_msun)
>      if(metal) ZloadSN_i = SNrecv(11,i)/SNrecv(5,i)
1052,1055c834
<      if(metal) ZloadSN_i = SNrecv(12,i)/SNrecv(5,i)
<      if(mechanical_geen.and.rt) rSt_i = SNrecv(13,i) ! Stromgren sphere in pc 
< 
<      M_SN_var = mSN_i*scale_msun/num_sn
---
>  
1068,1072d846
<            if(loading_type.eq.1.and.fload_i<f_LOAD)then
<               dratio = fload_i
<            else
<               dratio = 1d0
<            endif
1074,1076c848,849
<            nH_nei   = d_nei*scale_nH*dratio
<            Z_neisol = max(0.01,Z_nei/0.02)
<            Zdepen   = Z_neisol**(expZ_SN*2d0) !From Thornton+(98)
---
>            nH_nei   = d_nei*scale_nH
>            Zdepen   =(max(0.01,Z_nei/0.02))**(expZ_SN*2d0) !From Thornton+(98)
1079,1106c852,853
<            ! psn_tr = sqrt(2*chi_tr*Nsn*Esn*Msn*fe)
<            ! chi_tr = (1+f_w_crit)
<            psn_thor98   = A_SN * num_sn**(expE_SN) * nH_nei**(expN_SN) * Z_neisol**(expZ_SN)  !km/s Msun
<            psn_tr       = psn_thor98 
<            if(mechanical_geen)then
<               ! For snowplough phase, psn_tr will do the job
<               psn_geen15= A_SN_Geen * num_sn**(expE_SN) * Z_neisol**(expZ_SN)  !km/s Msun
<               if(rt)then
<                  fthor   = exp(-dx_loc_pc/rSt_i)
<                  psn_tr  = psn_thor98*fthor + psn_geen15*(1d0-fthor)
<               else
<                  psn_tr  = psn_geen15
<               endif
<               psn_tr = max(psn_tr, psn_thor98)
< 
<               ! For adiabatic phase
<               ! psn_tr =  A_SN * (E51 * boost_geen)**expE_SN_boost * nH_nei**(expN_SN_boost) * Z_neisol**(expZ_SN)
<               !        =  p_hydro * boost_geen**expE_SN_boost
<               p_hydro = A_SN * num_sn**(expN_SN_boost) * nH_nei**(expN_SN_boost) * Z_neisol**(expZ_SN)
<               boost_geen_ad = (psn_tr / p_hydro)**(1d0/expE_SN_boost)
<               boost_geen_ad = max(boost_geen_ad-1d0,0.0)
<            endif
<            chi_tr   = psn_tr**2d0 / (2d0 * num_sn * (E_SNII/msun2g/km2cm**2d0) * M_SN_var * f_ESN)
<            !          (Msun*km/s)^2                 (Msun *km2/s2)              (Msun)
<            f_w_crit = max(chi_tr-1d0, 0d0)
<  
<            !f_w_crit = (A_SN/1d4)**2d0/(f_ESN*M_SNII)*num_sn**((expE_SN-1d0)*2d0)*nH_nei**(expN_SN*2d0)*Zdepen - 1d0
<            !f_w_crit = max(0d0,f_w_crit)
---
>            f_w_crit = (A_SN/1d4)**2d0/(f_esn*M_SNII)*num_sn**((expE_SN-1d0)*2d0)*nH_nei**(expN_SN*2d0)*Zdepen - 1d0
>            f_w_crit = max(0d0,f_w_crit)
1108d854
<            vload_rad = dsqrt(2d0*f_ESN*E_SNII*(1d0+f_w_crit)/(M_SN_var*msun2g))/(1d0+f_w_cell)/scale_v/f_LOAD/f_PCAN
1110c856
<               vload = vload_rad
---
>               vload = dsqrt(2d0*f_esn*E_SNII*(1d0+f_w_crit)/(M_SNII*msun2g))/scale_v/(1d0+f_w_cell)/f_LOAD
1112,1122c858,860
<            else ! adiabatic phase
<               f_esn2 = 1d0-(1d0-f_ESN)*f_w_cell/f_w_crit
<               vload = dsqrt(2d0*f_esn2*E_SNII/(1d0+f_w_cell)/(M_SN_var*msun2g))/scale_v/f_LOAD
<               if(mechanical_geen) then
<                  f_wrt_snow = 2d0-2d0/(1d0+exp(-f_w_cell/f_w_crit/0.30))
<                  boost_geen_ad = boost_geen_ad * f_wrt_snow
<                  vload = vload * dsqrt(1d0+boost_geen_ad) ! f_boost
<                  ! NB. this sometimes give too much momentum because expE_SN_boost != expE_SN. A limiter is needed
<               endif
<               ! safety device: limit the maximum velocity so that it does not exceed p_{SN,final}
<               if(vload>vload_rad) vload = vload_rad
---
>            else
>               f_esn2 = 1d0-(1d0-f_esn)*f_w_cell/f_w_crit
>               vload = dsqrt(2d0*f_esn2*E_SNII/(1d0+f_w_cell)/(M_SNII*msun2g))/scale_v/f_LOAD
1142,1145c880
<      if(metal) ZloadSN_i = SNrecv(12,i)/mloadSN_i
<      do ich=1,nchem
<         chloadSN_i(ich) = SNrecv(13+ich,i)/mloadSN_i
<      end do
---
>      if(metal) ZloadSN_i = SNrecv(11,i)/SNrecv(5,i)
1148a884
>         
1156,1158d891
<               do ich=1,nchem
<                  pvar(ichem+ich-1) = unew(icell,ichem+ich-1)
<               end do
1167,1169d899
<               do ich=1,nchem
<                  pvar(ichem+ich-1) = uold(icell,ichem+ich-1)
<               end do
1188,1192c918,923
<            Tk0 = eth0/d0*scale_T2*(gamma-1)
<            T2min = T2_star*(d0*scale_nH/n_star)**g_star
<            if(Tk0<T2min)then
<               eth0 = T2min*d0/scale_T2/(gamma-1)
<            endif 
---
>            Tk0 = eth0/d0*scale_T2*(gamma-1)*0.6
> 
>            if(Tk0<0) then
>               print *,'TKERR : part1 mpi', myid,Tk0,d0*scale_nH
>               stop
>            endif
1214c945
< 
---
>            eturb = pvar(5) - ekk - eth0
1224,1226d954
<            do ich=1,nchem
<                pvar(ichem+ich-1)=pvar(ichem+ich-1)+mloadSN_i/dble(nSNnei)*chloadSN_i(ich)/vol_nei
<            end do
1237,1239d964
<               do ich=1,nchem
<                  unew(icell,ichem+ich-1) = pvar(ichem+ich-1)
<               end do
1259,1261d983
<               do ich=1,nchem
<                  uold(icell,ichem+ich-1) = pvar(ichem+ich-1)
<               end do
1286c1008
<   ncomm_SN=nSN_tot
---
>   nSN_comm=nSN_tot
1293c1015
< subroutine get_number_of_sn2(birth_time,dteff,zp_star,id_star,mass0,nsn,done_star)
---
> subroutine get_number_of_sn2(birth_time,zp_star,id_star,mass0,mass1,nsn,done_star)
1297,1299c1019,1020
<   real(kind=dp)::birth_time,zp_star,mass0,dteff ! birth_time in code, mass in Msun
<   real(kind=dp)::nsn
<   integer::nsn_tot,nsn_sofar
---
>   real(kind=dp)::birth_time,zp_star,mass0,mass1 ! birth_time in code, mass in Msun
>   integer::nsn,nsn_tot,nsn_sofar,nsn_ok,nsn_try
1301c1022
<   real(kind=dp)::age1,age2,tMyr,xdum,ydum,logzpsun
---
>   real(kind=dp)::age_star,tMyr,xdum,ydum,logzpsun
1325,1330c1046
<   call getStarAgeGyr(birth_time+dteff, age1)
<   call getStarAgeGyr(birth_time      , age2)
< 
<   ! convert Gyr -> Myr
<   age1 = age1*1d3
<   age2 = age2*1d3
---
>   call getStarAgeGyr(birth_time,age_star)
1332,1333c1048,1050
<   nsn=0d0; done_star=.false.
<   if(age2.le.(-co2/co1))then
---
>   tMyr=age_star*1d3 !Myr
>   nsn=0
>   if(tMyr.le.(-co2/co1))then
1339,1341c1056,1061
<   if(nsn_tot.eq.0)then
<       write(*,*) 'Fatal error: please increase the mass of your star particle'
<       stop
---
>   nsn_sofar = NINT((mass0-mass1)/M_SNII ,kind=4)
> 
>   if(sn2_real_delay)then
>      nsn_try=nsn_tot
>   else
>      nsn_try=1
1344c1064
<   nsn_sofar = 0
---
>   nsn_ok    = 0
1347c1067
<   do i=1,nsn_tot
---
>   do i=1,nsn_try 
1351,1355c1071,1072
<      if(ydum.ge.age1.and.ydum.le.age2) then
<         nsn=nsn+1
<      endif
<      if(ydum.le.age2)then
<         nsn_sofar=nsn_sofar+1
---
>      if(ydum.le.tMyr) then
>         nsn_ok=nsn_ok+1
1359c1076,1085
<   if(nsn_sofar.eq.nsn_tot) done_star=.true.
---
>   if(.not.sn2_real_delay) nsn_ok=nsn_ok*nsn_tot
> 
>   nsn = nsn_ok - nsn_sofar
>   if(nsn_tot.eq.0)then
>       write(*,*) 'Fatal error: please increase the mass of your star particle'
>       stop
>   endif
> 
>   done_star=.false.
>   if(nsn_ok.eq.nsn_tot) done_star=.true.
1524c1250
< end subroutine mesh_info
---
> end subroutine
1585,1667c1311
< end subroutine get_icell_from_pos
< !################################################################
< !################################################################
< !################################################################
< subroutine SNII_yield (zp_star, ej_m, ej_Z, ej_chem)
<   use amr_commons, ONLY:dp,nchem,chem_list
<   use hydro_parameters, ONLY:ichem 
<   implicit none
<   real(dp)::zp_star,ej_m,ej_Z,ej_chem(1:nchem)
< !-----------------------------------------------------------------
< ! Notice that even if the name is 'yield', 
< ! the return value is actually a metallicity fraction for simplicity
< ! These numbers are based on the outputs from Starburst99 
< !                   (i.e. essentially Woosley & Weaver 95)
< !-----------------------------------------------------------------
<   real(dp),dimension(1:5)::log_SNII_m, log_Zgrid, log_SNII_Z
<   real(dp),dimension(1:5)::log_SNII_H,log_SNII_He,log_SNII_C,log_SNII_N,log_SNII_O
<   real(dp),dimension(1:5)::log_SNII_Mg,log_SNII_Si,log_SNII_S,log_SNII_Fe, dum1d
<   real(dp)::log_Zstar,fz
<   integer::nz_SN=5, izg, ich
<   character(len=2)::element_name
< 
<   ! These are the numbers calculated from Starburst99 (Kroupa with 50Msun cut-off)
<   ! (check library/make_stellar_winds.pro)
<   log_SNII_m = (/-0.85591807,-0.93501857,-0.96138483,-1.0083450,-1.0544419/)
<   log_Zgrid = (/-3.3979400,-2.3979400,-2.0969100,-1.6989700,-1.3010300/)
<   log_SNII_Z=(/-0.99530662,-0.98262223,-0.96673581,-0.94018599,-0.93181853/)
<   log_SNII_H=(/-0.28525316, -0.28988675, -0.29588988, -0.30967822, -0.31088065/)
<   log_SNII_He=(/-0.41974152, -0.41688929, -0.41331511, -0.40330181, -0.40426739/)
<   log_SNII_C=(/-1.9739731, -1.9726015, -1.9692265, -1.9626259, -1.9635311/)
<   log_SNII_N=(/-3.7647616, -3.1325173, -2.8259748, -2.4260355, -2.1260417/)
<   log_SNII_O=(/-1.1596291, -1.1491227, -1.1344617, -1.1213435, -1.1202089/)
<   log_SNII_Mg=(/-2.4897201, -2.4368979, -2.4043654, -2.3706062, -2.3682933/)
<   log_SNII_Si=(/-2.2157073, -2.1895132, -2.1518758, -2.0431845, -2.0444829/)
<   log_SNII_S=(/-2.5492508, -2.5045016, -2.4482936, -2.2964020, -2.2988656/)
<   log_SNII_Fe=(/-2.0502141, -2.0702598, -2.1074876, -2.2126987, -2.3480877/)
< 
<   ! search for the metallicity index
<   log_Zstar = log10(zp_star)
<   call binary_search(log_Zgrid, log_Zstar, nz_SN, izg )
< 
<   fz  = (log_Zgrid(izg+1) - log_Zstar )/( log_Zgrid(izg+1) - log_Zgrid(izg) )
<   ! no extraploation
<   if (fz  < 0.0) fz  = 0.0
<   if (fz  > 1.0) fz  = 1.0
< 
< 
<   ej_m = log_SNII_m(izg)*fz + log_SNII_m(izg+1)*(1d0-fz)
<   ej_m = 10d0**ej_m
< 
<   ej_Z = log_SNII_Z(izg)*fz + log_SNII_Z(izg+1)*(1d0-fz)
<   ej_Z = 10d0**ej_Z
<  
<   do ich=1,nchem
<       element_name=chem_list(ich)
<       select case (element_name)
<          case ('H ')
<             dum1d = log_SNII_H
<          case ('He')
<             dum1d = log_SNII_He
<          case ('C ')
<             dum1d = log_SNII_C
<          case ('N ')
<             dum1d = log_SNII_N
<          case ('O ')
<             dum1d = log_SNII_O
<          case ('Mg')
<             dum1d = log_SNII_Mg
<          case ('Si')
<             dum1d = log_SNII_Si
<          case ('S ')
<             dum1d = log_SNII_S
<          case ('Fe')
<             dum1d = log_SNII_Fe
<          case default
<             dum1d = 0d0
<       end select
< 
<      ej_chem(ich) = dum1d(izg)*fz + dum1d(izg+1)*(1d0-fz)
<      ej_chem(ich) = 10d0**ej_chem(ich)
<   end do
< 
< end subroutine SNII_yield
---
> end subroutine
