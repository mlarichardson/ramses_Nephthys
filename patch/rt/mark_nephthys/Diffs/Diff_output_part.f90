4d3
< ! - writing of chp (Taysun, Jun 2016)
8a8,10
> #ifndef WITHOUTMPI
>   include 'mpif.h'  
> #endif 
11c13
<   integer::i,idim,ilun,ipart,ich
---
>   integer::i,idim,ilun,ipart
18a21,22
>   integer,parameter::tag=1122
>   integer::dummy_io,info2
21a26,36
>   ! Wait for the token
> #ifndef WITHOUTMPI
>   if(IOGROUPSIZE>0) then
>      if (mod(myid-1,IOGROUPSIZE)/=0) then
>         call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
>              & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
>      end if
>   endif
> #endif
> 
> 
146c161
<      endif ! if write_stellar_densities
---
>      endif
159,160d173
< 
<      ! Write initial mass
161a175
>         ! Write initial mass
171,184d184
< 
<      ! Write chemical abundance
<      do ich=1,nchem
<         ipart=0
<         do i=1,npartmax
<            if(levelp(i)>0)then
<               ipart=ipart+1
<               xdp(ipart)=chp(i,ich)
<            end if
<         end do
<         write(ilun)xdp
<      end do
< 
< 
187a188,199
> 
>   ! Send the token
> #ifndef WITHOUTMPI
>   if(IOGROUPSIZE>0) then
>      if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
>         dummy_io=1
>         call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
>              & MPI_COMM_WORLD,info2)
>      end if
>   endif
> #endif
>   
