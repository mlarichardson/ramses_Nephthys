6,9d5
< #ifdef RT
<   use rt_parameters, ONLY:iIons
<   use rt_cooling_module, ONLY:iIRtrapVar
< #endif
42,45c38
<        & ,interpol_var,interpol_type,sink_refine,&
<        & ,refine_extra,lmax_extra_refine,delta_lvl_extra_refine &
<        & ,ivar_refine_extra,var_cut_refine_extra,refine_extra_or
< 
---
>        & ,interpol_var,interpol_type,sink_refine
49c42
<        & ,d_bound,u_bound,v_bound,w_bound,p_bound
---
>        & ,d_bound,u_bound,v_bound,w_bound,p_bound,no_inflow
55c48
<        & ,units_density,units_time,units_length,neq_chem,ir_feedback,ir_eff,t_diss &
---
>        & ,units_density,units_time,units_length,neq_chem,ir_feedback,ir_eff,t_diss,t_sne &
59,61c52
<        & ,A_SN, expN_SN, E_SNII, E_SNIa, SF_kick_kms, write_stellar_densities      &
<        & ,use_initial_mass, stellar_winds, stellar_winds_file &
<        & ,chem_list, snIa, A_snIa, variable_yield_SNII, mechanical_geen, mass_loss_boost
---
>        & ,A_SN, expN_SN, E_SNII, SF_kick_kms, write_stellar_densities
136d126
<   if(stellar_winds )        use_initial_mass=.true. ! for stellar winds
147,150c137,139
<   else if(TRIM(star_imf).eq.'chabrier')then ! Chabrier05, not Chabrier03
<      M_SNII=19.134730 
<      eta_sn=0.317      
<      if(variable_yield_SNII.and.mass_loss_boost.le.1.0) mass_loss_boost=1.5
---
>   else if(TRIM(star_imf).eq.'chabrier')then
>      M_SNII=19.134730 !TODO
>      eta_sn=0.313706  !TODO
154a144
>   if(myid==1) write(*,*) '>>>TKNOTE: M_SNII, eta_sn, t_delay =',sngl(M_SNII), sngl(eta_sn), sngl(t_delay)
309,322d298
<   !  Deal with Extra Refine parameters
<   !-----------------------------------
<   if (lmax_extra_refine .gt. 0 .and. delta_lvl_extra_refine.lt.0) then
<     delta_lvl_extra_refine = nlevelmax - lmax_extra_refine
<   endif
< 
<   if (delta_lvl_extra_refine .lt. 0) then
<     if(myid==1)write(*,*)'Error in the namelist'
<     if(myid==1)write(*,*)'Delta_lvl_refine must be > 0 or'
<     if(myid==1)write(*,*)'  lmax_extra refine must be < lmax'
<     nml_ok=.false.
<   endif
< 
<   !-----------------------------------
325,329c301,304
<                                         !  nener   1   0
<   imetal=nener+ndim+3                   !  imetal  7   6
<   idelay=imetal                         !  idelay  7   6
<   if(metal)idelay=imetal+1              !  idelay  8   7
<   ixion=idelay                          !  ixion   8   7
---
>   imetal=nener+ndim+3
>   idelay=imetal
>   if(metal)idelay=imetal+1
>   ixion=idelay
331,332c306,307
<   ichem=ixion                           !  ichem   8   7
<   if(aton)ichem=ixion+1                 
---
>   ichem=ixion
>   if(aton)ichem=ixion+1
335,359d309
< #ifdef RT
<   ! This is necessary if the code is compiled with RT, but rt=.false.
<   ! Otherwise, xions in mechanical feedback tries to keep uold(iIons) the same,
<   ! which will lead to no change in metallicities
<   if (.not.rt) iIons=ichem+nchem        !  iIons       7 (for nchem=0)
<   ! if rt=.true., iIons will be updated in rt_init.f90
< #endif
< 
<   if(myid==1) then
<      write(*,*) '>>>NOTE: M_SNII    =',sngl(M_SNII)
<      write(*,*) '>>>NOTE: eta_sn    =',sngl(eta_sn)
<      write(*,*) '>>>NOTE: t_delay   =',sngl(t_delay)
<      write(*,*) '>>>NOTE: NSN/100M  =',sngl(eta_sn/M_SNII*100)
<      write(*,*) '>>>NOTE: imetal    =', imetal
<      do i=1,nchem
<         write(*,*) '>>>NOTE: ichem ('//chem_list(i)//')=', ichem+i-1
<      end do
< #ifdef RT
<      if(.not.rt)then
<      write(*,*) '>>>NOTE: iIons  = ', iIons
<      write(*,*) '>>>NOTE: iIRtrapVar  = ', iIRtrapVar
<      endif
< #endif
<      write(*,*)
<   endif
