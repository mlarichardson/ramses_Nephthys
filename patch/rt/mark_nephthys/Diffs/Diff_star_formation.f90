39c39
<   integer  :: ntot,ntot_all,info,nstar_corrected,ich,ncell
---
>   integer  :: ntot,ntot_all,info,nstar_corrected,iche,ncell
62a63
> #ifdef NCHEM
63a65
> #endif
701,702c703,705
<               do ich=1,nchem
<                  chem1(ich) = uold(ind_cell_new(i),ichem+ich-1)
---
> #ifdef NCHEM
>               do iche=1,nchem
>                  chem1(iche) = uold(ind_cell_new(i),ichem+iche-1)
703a707
> #endif
750,751c754,756
<               do ich=1,nchem
<                  chp(ind_part(i),ich) = chem1(ich)  ! Initial chemical abudance
---
> #ifdef NCHEM
>               do iche=1,nchem
>                  chp(ind_part(i),iche) = chem1(iche)  ! Initial chemical abudance
752a758
> #endif
