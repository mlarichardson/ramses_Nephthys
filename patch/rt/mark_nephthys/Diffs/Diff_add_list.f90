4d3
< ! Added  chp (Taysun, Jun 2016)
50c49
<   integer::np,ich
---
>   integer::np
84,88d82
<      do ich=1,nchem
<         do j=1,np
<            chp(ind_part(j),ich)=0.0
<         end do
<      end do
119c113
<   integer::np,ich
---
>   integer::np
162,168d155
<      do ich=1,nchem
<         do j=1,np
<            if(ok(j))then
<              chp(ind_part(j),ich)=0.0
<            endif
<         end do
<      end do
