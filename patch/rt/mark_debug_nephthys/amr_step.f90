recursive subroutine amr_step(ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
#ifdef RT
  use rt_hydro_commons
  use SED_module
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,icount,ilmax
  double precision::lmax_real
  !-------------------------------------------------------------------!
  ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
  ! Each routine is called using a specific order, don't change it,   !
  ! unless you check all consequences first                           !
  !-------------------------------------------------------------------!
  integer::i,idim,ivar
  logical::ok_defrag
  logical,save::first_step=.true.

  if(numbtot(1,ilevel)==0)return

  if(verbose)write(*,999)icount,ilevel

  !-------------------------------------------
  ! Make new refinements and update boundaries
  !-------------------------------------------
  if(levelmin.lt.nlevelmax .and..not. static)then
     if(ilevel==levelmin.or.icount>1)then
        do i=ilevel,nlevelmax
           if(i>levelmin)then

              !--------------------------
              ! Build communicators
              !--------------------------
              call build_comm(i)

              !--------------------------
              ! Update boundaries
              !--------------------------
              call make_virtual_fine_int(cpu_map(1),i)
              if(hydro)then
#ifdef SOLVERmhd
                 do ivar=1,nvar+3
#else
                 do ivar=1,nvar
#endif
                    call make_virtual_fine_dp(uold(1,ivar),i)
#ifdef SOLVERmhd
                 end do
#else
                 end do
#endif
                 if(simple_boundary)call make_boundary_hydro(i)
              end if
#ifdef RT
              if(rt)then
                 do ivar=1,nrtvar
                    call make_virtual_fine_dp(rtuold(1,ivar),i)
                 end do
                 if(simple_boundary)call rt_make_boundary_hydro(i)
              end if
#endif
              if(poisson)then
                 call make_virtual_fine_dp(phi(1),i)
                 do idim=1,ndim
                    call make_virtual_fine_dp(f(1,idim),i)
                 end do
                 if(simple_boundary)call make_boundary_force(i)
              end if
           end if

           !--------------------------
           ! Refine grids
           !--------------------------
           call refine_fine(i)
        end do
     end if
  end if

  !--------------------------
  ! Load balance
  !--------------------------
  ok_defrag=.false.
  if(levelmin.lt.nlevelmax)then
     if(ilevel==levelmin)then
        if(nremap>0)then
           ! Skip first load balance because it has been performed before file dump
           if(nrestart>0.and.first_step)then
              first_step=.false.
           else
              if(MOD(nstep_coarse,nremap)==0)then
                 call load_balance
                 call defrag
                 ok_defrag=.true.
              endif
           end if
        end if
     endif
  end if

  !-----------------
  ! Update sink cloud particle properties
  !-----------------
  if(sink)call update_cloud(ilevel)

  !-----------------
  ! Particle leakage
  !-----------------
  if(pic)call make_tree_fine(ilevel)
  call check_igrid(ilevel,19)
  
  !------------------------
  ! Output results to files
  !------------------------
  if(ilevel==levelmin)then
     if(mod(nstep_coarse,foutput)==0.or.aexp>=aout(iout).or.t>=tout(iout).or.output_now.EQV..true.)then
        if(.not.ok_defrag)then
           call defrag
        endif

        call dump_all

        ! Run the clumpfinder, (produce output, don't keep arrays alive on output)
        if(clumpfind .and. ndim==3) call clump_finder(.true.,.false.)

        ! Dump lightcone
        if(lightcone) call output_cone()

        if (output_now.EQV..true.) then
          output_now=.false.
        endif

     endif

     ! Important can't be done in sink routines because it must be done after dump all
     if(sink)acc_rate=0.

  endif

  !----------------------------
  ! Output frame to movie dump (without synced levels)
  !----------------------------
  ! call check_igrid(ilevel,18)
  if(movie) then
     if(imov.le.imovout)then ! ifort returns error for next statement if looking
                             ! beyond what is allocated as an array amovout/tmovout
        if(aexp>=amovout(imov).or.t>=tmovout(imov))then
           call output_frame()
        endif
     endif
  end if

  !-----------------------------------------------------------
  ! Put here all stuffs that are done only at coarse time step
  !-----------------------------------------------------------
  ! call check_igrid(ilevel,17)
  if(ilevel==levelmin)then
     if (hydro .and. star .and. eta_sn>0 .and. f_w>0) then
        !----------------------------------------------------
        ! Kinetic feedback
        !----------------------------------------------------
        call kinetic_feedback 
        
     endif
     
  endif


  !--------------------
  ! Poisson source term
  !--------------------
  ! call check_igrid(ilevel,16)
  if(poisson)then
     !save old potential for time-extrapolation at level boundaries
     call save_phi_old(ilevel)
     call rho_fine(ilevel,icount)
  endif

  !-------------------------------------------
  ! Sort particles between ilevel and ilevel+1
  !-------------------------------------------
  ! call check_igrid(ilevel,15)
  if(pic)then
     ! Remove particles to finer levels
     call kill_tree_fine(ilevel)
     ! Update boundary conditions for remaining particles
     call virtual_tree_fine(ilevel)
  end if

  !---------------
  ! Gravity update
  !---------------
  ! call check_igrid(ilevel,14)
  if(poisson)then
 
     ! Remove gravity source term with half time step and old force
     if(hydro)then
        call synchro_hydro_fine(ilevel,-0.5*dtnew(ilevel))
     endif
     
     ! call check_igrid(ilevel,13)
     ! Compute gravitational potential
     if(ilevel>levelmin)then
        if(ilevel .ge. cg_levelmin) then
           call phi_fine_cg(ilevel,icount)
        else
           call multigrid_fine(ilevel,icount)
        end if
     else
        call multigrid_fine(levelmin,icount)
     end if
     !when there is no old potential...
     if (nstep==0)call save_phi_old(ilevel)
     ! call check_igrid(ilevel,12)

     ! Compute gravitational acceleration
     call force_fine(ilevel,icount)
     ! call check_igrid(ilevel,11)

     ! Mechanical feedback from stars
     if(hydro.and.star.and.eta_sn>0.and.mechanical_feedback) then

        call mechanical_feedback_fine(ilevel,icount)
        if (snIa) call mechanical_feedback_snIa_fine(ilevel,icount)
     endif

     ! Synchronize remaining particles for gravity
     if(pic)then
        call synchro_fine(ilevel)
     end if
     ! call check_igrid(ilevel,10)

     if(hydro)then

        ! Add gravity source term with half time step and new force
        call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))
        ! call check_igrid(ilevel,9)

        ! Update boundaries
#ifdef SOLVERmhd
        do ivar=1,nvar+3
#else
        do ivar=1,nvar
#endif
           call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
        end do
#else
        end do
#endif
        if(simple_boundary)call make_boundary_hydro(ilevel)
        ! call check_igrid(ilevel,8)
        
        ! Compute Bondi-Hoyle accretion parameters
        if(sink)call collect_acczone_avg(ilevel)
        ! call check_igrid(ilevel,7)

     end if
  end if

#ifdef RT
  ! Turn on RT in case of rt_stars and first stars just created:
  ! Update photon packages according to star particles
  if(rt .and. rt_star) call update_star_RT_feedback(ilevel)
#endif

  ! call check_igrid(ilevel,6)
  !----------------------
  ! Compute new time step
  !----------------------
  call newdt_fine(ilevel)
  if(ilevel>levelmin)then
     dtnew(ilevel)=MIN(dtnew(ilevel-1)/real(nsubcycle(ilevel-1)),dtnew(ilevel))
  end if

  ! call check_igrid(ilevel,5)
  ! Set unew equal to uold
  if(hydro)call set_unew(ilevel)
  ! call check_igrid(ilevel,4)
#ifdef RT
  ! Set rtunew equal to rtuold
  if(rt)call rt_set_unew(ilevel)
#endif

  !---------------------------
  ! Recursive call to amr_step
  !---------------------------
  ! call check_igrid(ilevel,3)
  if(ilevel<nlevelmax)then
     if(numbtot(1,ilevel+1)>0)then
        if(nsubcycle(ilevel)==2)then
           call amr_step(ilevel+1,1)
           call amr_step(ilevel+1,2)
        else
           call amr_step(ilevel+1,1)
        endif
     else 
        ! Otherwise, update time and finer level time-step
        dtold(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        dtnew(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        call update_time(ilevel)
        if(sink)call update_sink(ilevel)
     end if
  else
     call update_time(ilevel)
     if(sink)call update_sink(ilevel)
  end if

  ! call check_igrid(ilevel,2)
  ! Thermal feedback from stars (also call if no feedback, for bookkeeping)
  if(hydro.and.star.and.f_w==0.0) then
     call thermal_feedback(ilevel)
  endif

  ! Check neighbours
  ilmax = int( log(0.5d0**0.333333d0/aexp)/log(2.d0) + 1.d0 )
  lmax_real=nlevelmax-ilmax

  if (ilevel .ge. 9) call check_neighbours_unew_fine(ilevel)

  ! Stellar winds
  call check_igrid(ilevel,1)
  if(hydro.and.star.and.stellar_winds) then
     call stellar_winds_fine(ilevel)
  endif
  call check_igrid(ilevel,0)

#ifdef RT
  ! Add stellar radiation sources
  if(rt.and.rt_star) call star_RT_feedback(ilevel,dtnew(ilevel))
#endif
  
  ! call check_igrid(ilevel,-1)
  ! Density threshold or Bondi accretion onto sink particle
  if(sink)then
     call grow_sink(ilevel,.false.)
  end if
  
  !---------------
  ! Move particles
  !---------------
  if(pic)then
     call move_fine(ilevel) ! Only remaining particles
  end if

  !-----------
  ! Hydro step
  !-----------
  if(hydro)then

     ! Hyperbolic solver
     call godunov_fine(ilevel)

     ! Reverse update boundaries
#ifdef SOLVERmhd
     do ivar=1,nvar+3
#else
     do ivar=1,nvar
#endif
        call make_virtual_reverse_dp(unew(1,ivar),ilevel)
#ifdef SOLVERmhd
     end do
#else
     end do
#endif
     if(pressure_fix)then
        call make_virtual_reverse_dp(enew(1),ilevel)
        call make_virtual_reverse_dp(divu(1),ilevel)
     endif

     ! Set uold equal to unew
     call set_uold(ilevel)

     ! ! Density threshold or Bondi accretion onto sink particle
     ! if(sink)then
     !    !this is a trick to temporarily solve the issue with sink accretion 
     !    !from ghost zones. Only an option for simulations without dark matter.
     !    if (.not. cosmo)then
     !       call make_tree_fine(ilevel)
     !       call virtual_tree_fine(ilevel)
     !       ! assuming all sink cloud parts sit on levelmax 
     !       ! it's better to compute the accretion_rate based on
     !       ! the updated values
     !       call collect_acczone_avg(ilevel)
     !    end if
     !    call grow_sink(ilevel,.false.)
     ! end if

     ! Add gravity source term with half time step and old force
     ! in order to complete the time step 
     if(poisson)call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel))

     ! Restriction operator
     call upload_fine(ilevel)

  endif

#ifdef RT
  !---------------
  ! Radiation step
  !---------------
  if(rt)then
     ! Hyperbolic solver
     if(rt_advect) call rt_godunov_fine(ilevel,dtnew(ilevel))

     call add_rt_sources(ilevel,dtnew(ilevel))

     ! Reverse update boundaries
     do ivar=1,nrtvar
        call make_virtual_reverse_dp(rtunew(1,ivar),ilevel)
     end do

     ! Set rtuold equal to rtunew
     call rt_set_uold(ilevel)

     ! Restriction operator
     call rt_upload_fine(ilevel)
  endif
#endif
  
  !-------------------------------
  ! Source term in leaf cells only
  !-------------------------------
  if(neq_chem.or.cooling.or.T2_star>0.0)call cooling_fine(ilevel)

  !----------------------------------
  ! Star formation in leaf cells only
  !----------------------------------
  if(hydro.and.star)call star_formation(ilevel)

  !---------------------------------------
  ! Update physical and virtual boundaries
  !---------------------------------------
  if(hydro)then
#ifdef SOLVERmhd
     do ivar=1,nvar+3
#else
     do ivar=1,nvar
#endif
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
     end do
#else
     end do
#endif
     if(simple_boundary)call make_boundary_hydro(ilevel)
  endif
#ifdef RT
  if(rt)then
     do ivar=1,nrtvar
        call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
     end do
     if(simple_boundary)call rt_make_boundary_hydro(ilevel)
  end if
#endif

#ifdef SOLVERmhd
  ! Magnetic diffusion step
 if(hydro)then
     if(eta_mag>0d0.and.ilevel==levelmin)then
        call diffusion
     endif
  end if
#endif

  !-----------------------
  ! Compute refinement map
  !-----------------------
  if(.not.static) call flag_fine(ilevel,icount)


  !----------------------------
  ! Merge finer level particles
  !----------------------------
  if(pic)call merge_tree_fine(ilevel)

  !---------------
  ! Radiation step
  !---------------
#ifdef ATON
  if(aton.and.ilevel==levelmin)then
     call rad_step(dtnew(ilevel))
  endif
#endif

  if(sink)then
     !-------------------------------
     ! Update coarser level sink velocity
     !-------------------------------
     if(ilevel>levelmin)then
        vsold(1:nsink,1:ndim,ilevel-1)=vsnew(1:nsink,1:ndim,ilevel-1)
        if(nsubcycle(ilevel-1)==1)vsnew(1:nsink,1:ndim,ilevel-1)=vsnew(1:nsink,1:ndim,ilevel)
        if(icount==2)vsnew(1:nsink,1:ndim,ilevel-1)= &
             (vsold(1:nsink,1:ndim,ilevel)*dtold(ilevel)+vsnew(1:nsink,1:ndim,ilevel)*dtnew(ilevel))/ &
             (dtold(ilevel)+dtnew(ilevel))
     end if
     !---------------
     ! Sink production
     !---------------
     if(ilevel==levelmin)call create_sink
  end if

  !-------------------------------
  ! Update coarser level time-step
  !-------------------------------
  if(ilevel>levelmin)then
     if(nsubcycle(ilevel-1)==1)dtnew(ilevel-1)=dtnew(ilevel)
     if(icount==2)dtnew(ilevel-1)=dtold(ilevel)+dtnew(ilevel)
  end if

999 format(' Entering amr_step',i1,' for level',i2)

end subroutine amr_step



