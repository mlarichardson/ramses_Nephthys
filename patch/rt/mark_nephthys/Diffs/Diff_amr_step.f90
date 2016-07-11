8a9,11
>   use UV_module
>   use coolrates_module, only: update_coolrates_tables
>   use rt_cooling_module, only: update_UVrates
30a34
>                                call timer('refine','start')
86a91
>                                call timer('loadbalance','start')
107a113
>                                call timer('sinks','start')
112a119
>                                call timer('particles','start')
119a127
>                                call timer('io','start')
149a158
>                                call timer('io','start')
159c168,170
<      if (hydro .and. star .and. eta_sn>0 .and. f_w>0) then
---
>                                call timer('feedback','start')
>      if (hydro .and. star .and. eta_sn>0 .and. f_w>0                &
>              .and.(mechanical_feedback.eq.0)) then
169,194d179
<   ! Mechanical feedback from stars
<   if(hydro.and.star.and.eta_sn>0.and.mechanical_feedback) then
< 
<       call mechanical_feedback_fine(ilevel,icount)
<       if (snIa) call mechanical_feedback_snIa_fine(ilevel,icount)
< 
<       if(.not.poisson)then
< #ifdef SOLVERmhd
<         do ivar=1,nvar+3
< #else
<         do ivar=1,nvar
< #endif
< 
<            call make_virtual_fine_dp(uold(1,ivar),ilevel)
< #ifdef SOLVERmhd
<         end do
< #else
<         end do
< #endif
<      end if
< 
<   endif
< 
< 
< 
< 
198a184
>                                call timer('poisson','start')
208a195
>                                call timer('particles','start')
217a205
>                                call timer('poisson','start')
239a228,234
>      ! Thermal feedback from stars
>                                call timer('feedback','start')
>      if(hydro.and.star.and.eta_sn>0 &
>        .and. (mechanical_feedback.eq.2)) then
>         call mechanical_feedback_fine(ilevel,icount)
>      endif
> 
241a237
>                                call timer('particles','start')
245a242
>                                call timer('poisson','start')
264a262
>                                call timer('sinks','start')
272a271
>                                call timer('radiative transfer','start')
278a278
>                                call timer('courant','start')
284a285
>                                call timer('hydro - set unew','start')
288a290
>                                call timer('radiative transfer','start')
316c318,320
<   if(hydro.and.star.and.f_w==0.0) then
---
>   if(hydro.and.star.and.f_w==0.0 &
>        .and.(mechanical_feedback.eq.0)) then
>                                call timer('feedback','start')
320,330d323
<   ! Stellar winds
<   if(hydro.and.star.and.stellar_winds) then
<      call stellar_winds_fine(ilevel)
<   endif
< 
< 
< #ifdef RT
<   ! Add stellar radiation sources
<   if(rt.and.rt_star) call star_RT_feedback(ilevel,dtnew(ilevel))
< #endif
<   
332a326
>                                call timer('sinks','start')
336,342d329
<   !---------------
<   ! Move particles
<   !---------------
<   if(pic)then
<      call move_fine(ilevel) ! Only remaining particles
<   end if
< 
348a336
>                                call timer('hydro - godunov','start')
351a340
>                                call timer('hydro - rev ghostzones','start')
368a358
>                                call timer('hydro - set uold','start')
371,385d360
<      ! ! Density threshold or Bondi accretion onto sink particle
<      ! if(sink)then
<      !    !this is a trick to temporarily solve the issue with sink accretion 
<      !    !from ghost zones. Only an option for simulations without dark matter.
<      !    if (.not. cosmo)then
<      !       call make_tree_fine(ilevel)
<      !       call virtual_tree_fine(ilevel)
<      !       ! assuming all sink cloud parts sit on levelmax 
<      !       ! it's better to compute the accretion_rate based on
<      !       ! the updated values
<      !       call collect_acczone_avg(ilevel)
<      !    end if
<      !    call grow_sink(ilevel,.false.)
<      ! end if
< 
387a363
>                                call timer('poisson','start')
390a367
>                                call timer('hydro upload fine','start')
393a371
>  
394a373,375
>   !---------------------
>   ! Do RT/Chemistry step
>   !---------------------
396,414c377,384
<   !---------------
<   ! Radiation step
<   !---------------
<   if(rt)then
<      ! Hyperbolic solver
<      if(rt_advect) call rt_godunov_fine(ilevel,dtnew(ilevel))
< 
<      call add_rt_sources(ilevel,dtnew(ilevel))
< 
<      ! Reverse update boundaries
<      do ivar=1,nrtvar
<         call make_virtual_reverse_dp(rtunew(1,ivar),ilevel)
<      end do
< 
<      ! Set rtuold equal to rtunew
<      call rt_set_uold(ilevel)
< 
<      ! Restriction operator
<      call rt_upload_fine(ilevel)
---
>   if(rt .and. rt_advect) then  
>                                call timer('radiative transfer','start')
>      call rt_step(ilevel)
>   else
>      ! Still need a chemistry call if RT is defined but not
>      ! actually doing radiative transfer (i.e. rt==false):
>                                call timer('cooling','start')
>      if(neq_chem.or.cooling.or.T2_star>0.0)call cooling_fine(ilevel)
415a386,399
>   ! Regular updates and book-keeping:
>   if(ilevel==levelmin) then
>                                call timer('radiative transfer','start')
>      if(cosmo) call update_rt_c
>      if(cosmo .and. haardt_madau) call update_UVrates(aexp)
>      if(cosmo .and. rt_isDiffuseUVsrc) call update_UVsrc
>                                call timer('cooling','start')
>      if(cosmo) call update_coolrates_tables(dble(aexp))
>                                call timer('radiative transfer','start')
>      if(ilevel==levelmin) call output_rt_stats
>   endif
> #else
>                                call timer('cooling','start')
>   if(neq_chem.or.cooling.or.T2_star>0.0)call cooling_fine(ilevel)
418,421c402,408
<   !-------------------------------
<   ! Source term in leaf cells only
<   !-------------------------------
<   if(neq_chem.or.cooling.or.T2_star>0.0)call cooling_fine(ilevel)
---
>   !---------------
>   ! Move particles
>   !---------------
>   if(pic)then
>                                call timer('particles','start')
>      call move_fine(ilevel) ! Only remaining particles
>   end if
425a413
>                                call timer('feedback','start')
431a420
>                                call timer('hydro - ghostzones','start')
445,452d433
< #ifdef RT
<   if(rt)then
<      do ivar=1,nrtvar
<         call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
<      end do
<      if(simple_boundary)call rt_make_boundary_hydro(ilevel)
<   end if
< #endif
457a439
>                                call timer('hydro - diffusion','start')
465a448
>                                call timer('flag','start')
471a455
>                                call timer('particles','start')
478a463
>                                call timer('aton','start')
483a469
>                                call timer('sinks','start')
511a498,501
> !##########################################################################
> !##########################################################################
> !##########################################################################
> !##########################################################################
512a503,516
> #ifdef RT
> subroutine rt_step(ilevel)
>   use amr_parameters, only: dp
>   use amr_commons,    only: levelmin, t, dtnew, myid
>   use rt_parameters, only: rt_isDiffuseUVsrc
>   use rt_cooling_module, only: update_UVrates
>   use rt_hydro_commons
>   use UV_module
>   use SED_module,     only: star_RT_feedback
>   implicit none
> #ifndef WITHOUTMPI
>   include 'mpif.h'
> #endif
>   integer, intent(in) :: ilevel
513a518,540
> !--------------------------------------------------------------------------
> !  Radiative transfer and chemistry step. Either do one step on ilevel,
> !  with radiation field updates in coarser level neighbours, or, if
> !  rt_nsubsteps>1, do many substeps in ilevel only, using Dirichlet
> !  boundary conditions for the level boundaries. 
> !--------------------------------------------------------------------------
> 
>   real(dp) :: dt_hydro, t_left, dt_rt, t_save
>   integer  :: i_substep, ivar
> 
>   dt_hydro = dtnew(ilevel)                   ! Store hydro timestep length
>   t_left = dt_hydro
>   ! We shift the time backwards one hydro-dt, to get evolution of stellar
>   ! ages within the hydro timestep, in the case of rt subcycling:
>   t_save=t ; t=t-t_left
>   
>   i_substep = 0
>   do while (t_left > 0)                      !                RT sub-cycle
>      i_substep = i_substep + 1
>      call get_rt_courant_coarse(dt_rt)
>      ! Temporarily change timestep length to rt step:
>      dtnew(ilevel) = MIN(t_left, dt_rt/2.0**(ilevel-levelmin))
>      t = t + dtnew(ilevel) ! Shift the time forwards one dt_rt
514a542,581
>      ! If (myid==1) write(*,900) dt_hydro, dtnew(ilevel), i_substep, ilevel    
>      if (i_substep > 1) call rt_set_unew(ilevel)
> 
>      if(rt_star) call star_RT_feedback(ilevel,dtnew(ilevel))
> 
>      ! Hyperbolic solver
>      if(rt_advect) call rt_godunov_fine(ilevel,dtnew(ilevel))
> 
>      call add_rt_sources(ilevel,dtnew(ilevel))
> 
>      ! Reverse update boundaries
>      do ivar=1,nrtvar
>         call make_virtual_reverse_dp(rtunew(1,ivar),ilevel)
>      end do
> 
>      ! Set rtuold equal to rtunew
>      call rt_set_uold(ilevel)
> 
>      if(neq_chem.or.cooling.or.T2_star>0.0)call cooling_fine(ilevel)
>      
>      do ivar=1,nrtvar
>         call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
>      end do
>      if(simple_boundary)call rt_make_boundary_hydro(ilevel)
> 
>      t_left = t_left - dtnew(ilevel)
>   end do                                   !          End RT subcycle loop
>   dtnew(ilevel) = dt_hydro                 ! Restore hydro timestep length
>   t = t_save       ! Restore original time (otherwise tiny roundoff error)
>   
>   ! Restriction operator to update coarser level split cells
>   call rt_upload_fine(ilevel)
> 
>   if (myid==1 .and. rt_nsubcycle .gt. 1) write(*,901) ilevel, i_substep
> 
> 900 format (' dt_hydro=', 1pe12.3, ' dt_rt=', 1pe12.3, ' i_sub=', I5, ' level=', I5)
> 901 format (' Performed level', I3, ' RT-step with ', I5, ' subcycles')
>   
> end subroutine rt_step
> #endif
