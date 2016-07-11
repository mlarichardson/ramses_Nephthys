38c38,41
<        & ,interpol_var,interpol_type,sink_refine
---
>        & ,interpol_var,interpol_type,sink_refine &
>        & ,refine_extra,lmax_extra_refine,delta_lvl_extra_refine &
>        & ,ivar_refine_extra,var_cut_refine_extra,refine_extra_or
> 
298a302,315
>   !  Deal with Extra Refine parameters
>   !-----------------------------------
>   if (lmax_extra_refine .gt. 0 .and. delta_lvl_extra_refine.lt.0) then
>     delta_lvl_extra_refine = nlevelmax - lmax_extra_refine
>   endif
> 
>   if (delta_lvl_extra_refine .lt. 0) then
>     if(myid==1)write(*,*)'Error in the namelist'
>     if(myid==1)write(*,*)'Delta_lvl_refine must be > 0 or'
>     if(myid==1)write(*,*)'  lmax_extra refine must be < lmax'
>     nml_ok=.false.
>   endif
> 
>   !-----------------------------------
