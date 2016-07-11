118a119
>   real(dp)::n_turb =0.1D0     ! Star formation density threshold for checking turbulent support, in H/cc
200a202,211
>   ! second refinement criteria / level
>   integer::refine_extra=0 ! Number of extra variables to use for refinement
>   integer::lmax_extra_refine=-1 ! Maximum level of refinement for the second variables.
>                                 ! THIS MUST BE LESS THAN NLEVELMAX currently
>   integer::delta_lvl_extra_refine=-1 ! Or, you can put in a offset from max level.
>   integer,parameter::MAXNUMBER_EXTRAREFINE=10 ! Max Number of exra variables for refine
>   integer,dimension(MAXNUMBER_EXTRAREFINE)::ivar_refine_extra=-1 ! variable types
>   real(dp),dimension(MAXNUMBER_EXTRAREFINE)::var_cut_refine_extra=-1.d0 ! thresholds
>   ! NOTE: Set threshold to negative to be an upper threshold.
>   logical::refine_extra_or=.true.
