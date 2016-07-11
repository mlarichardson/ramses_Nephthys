129c129
<   real(dp)::f_w    =0.0D0    ! Supernovae mass loading factor
---
>   real(dp)::f_w    =0.0D0     ! Supernovae mass loading factor
137a138
>   real(dp)::t_sne =10.0D0     ! Supernova blast time
183a185
>   character(len=5),dimension(0:NVAR+6)::movie_vars_txt=''
185a188
>   character(len=5),dimension(0:NVAR+2)::movie_vars_txt=''
202,212d204
<   ! second refinement criteria / level
<   integer::refine_extra=0 ! Number of extra variables to use for refinement
<   integer::lmax_extra_refine=-1 ! Maximum level of refinement for the second variables.
<                                 ! THIS MUST BE LESS THAN NLEVELMAX currently
<   integer::delta_lvl_extra_refine=-1 ! Or, you can put in a offset from max level.
<   integer,parameter::MAXNUMBER_EXTRAREFINE=10 ! Max Number of exra variables for refine
<   integer,dimension(MAXNUMBER_EXTRAREFINE)::ivar_refine_extra=-1 ! variable types
<   real(dp),dimension(MAXNUMBER_EXTRAREFINE)::var_cut_refine_extra=-1.d0 ! thresholds
<   ! NOTE: Set threshold to negative to be an upper threshold.
<   logical::refine_extra_or=.true.
< 
247c239
<   logical                           ::no_inflow
---
>   logical                           ::no_inflow=.false.
266,268c258,259
<   real(dp)::M_SNII=10d0       ! mean progenitor mass of the TypeII SNe
<   real(dp)::E_SNII=1d51       ! typical energy from Type II SNe
<   real(dp)::E_SNIa=1d51       ! typical energy from Type Ia SNe 
---
> !  real(dp)::eff_sfbk=1.0D0
>   real(dp)::E_SNII=1d51        ! different from ESN used in feedback.f90 
273,275c264,265
<   ! Activate mechanical feedback:
<   logical ::mechanical_feedback=.false.
<   logical ::mechanical_geen=.false.
---
>   ! Activate mechanical feedback (cannot use with star_particle_winds):
>   integer ::mechanical_feedback=0
278,279c268
<   real(dp)::A_SN=2.5d5
<   real(dp)::A_SN_Geen=5d5 
---
>   real(dp)::A_SN=3d5
285a275
>   real(dp)::M_SNII=10d0       ! Mean progenitor mass of the TypeII SNe
293,307d282
< 
<   ! Stellar winds stuff:
< #ifndef NCHEM
<   integer,parameter::nchem=0 ! number of chemical elements (max: 8)
< #else
<   integer,parameter::nchem=NCHEM
< #endif
<   logical::stellar_winds=.false.
<   character(LEN=256)::stellar_winds_file='/home/kimm/soft/lib/swind_krp_pagb.dat'
<   character(LEN=2),dimension(1:8):: chem_list=(/'H ','O ','Fe','Mg','C ','N ','Si','S '/)
<   ! SN Type Ia
<   logical ::snIa=.false.
<   real(dp)::A_snIa=0.05
<   logical ::variable_yield_SNII=.false.  ! TypeII yields are computed according to metallicities based on starburst99
<   real(dp)::mass_loss_boost=1d0
