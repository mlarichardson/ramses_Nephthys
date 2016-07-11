4d3
< ! Taysun added chp (Jun 2016) 
856c855
<   integer::i,idim,ich
---
>   integer::i,idim
914,919d912
<      do ich=1,nchem
<         do i=1,np
<            reception(icpu,ilevel)%up(ind_com(i),current_property)=chp(ind_part(i),ich)
<         end do
<         current_property = current_property+1
<      end do
