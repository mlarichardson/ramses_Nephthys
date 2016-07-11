2d1
< ! Patch change: initialisation of chp
24c23
<   integer::buf_count,indglob,npart_new,ich
---
>   integer::buf_count,indglob,npart_new
59c58,60
<   character(LEN=5)::nchar
---
>   character(LEN=5)::nchar,ncharcpu
>   integer,parameter::tagg=1109,tagg2=1110,tagg3=1111
>   integer::dummy_io,info2
91c92
<      if(metal)then
---
>         if(metal)then
95,98d95
<      if(nchem>0)then
<         allocate(chp(npartmax,nchem))
<         chp=0.0
<      endif
113c110,117
<      fileloc='output_'//TRIM(nchar)//'/part_'//TRIM(nchar)//'.out'
---
> 
>      if(IOGROUPSIZEREP>0)then
>         call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
>         fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/part_'//TRIM(nchar)//'.out'
>      else
>         fileloc='output_'//TRIM(nchar)//'/part_'//TRIM(nchar)//'.out'
>      endif
> 
115a120,129
>      ! Wait for the token                                                                                                                                                                    
> #ifndef WITHOUTMPI
>      if(IOGROUPSIZE>0) then
>         if (mod(myid-1,IOGROUPSIZE)/=0) then
>            call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tagg,&
>                 & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
>         end if
>      endif
> #endif
> 
184,188d197
<         do ich=1,nchem
<            ! Read chemical abundance
<            read(ilun)xdp
<            chp(1:npart2,ich)=xdp
<         end do
195a205,216
> 
>      ! Send the token      
> #ifndef WITHOUTMPI
>      if(IOGROUPSIZE>0) then
>         if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
>            dummy_io=1
>            call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tagg, &
>                 & MPI_COMM_WORLD,info2)
>         end if
>      endif
> #endif
> 
348a370,378
>                  ! Wait for the token                                                                                                                                                        
> #ifndef WITHOUTMPI
>                  if(IOGROUPSIZE>0) then
>                     if (mod(myid-1,IOGROUPSIZE)/=0) then
>                        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tagg2,&
>                             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
>                     end if
>                  endif
> #endif
361a392,402
>                  ! Send the token                                                                                                                                                            
> #ifndef WITHOUTMPI
>                  if(IOGROUPSIZE>0) then
>                     if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
>                        dummy_io=1
>                        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tagg2, &
>                             & MPI_COMM_WORLD,info2)
>                     end if
>                  endif
> #endif
> 
660,662d700
<               do ich=1,nchem
<                  chp(ipart,ich)=0d0
<               end do
