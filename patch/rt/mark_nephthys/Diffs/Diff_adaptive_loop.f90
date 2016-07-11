16,17c16,18
<   integer::ilevel,idim,ivar,info
<   real(kind=8)::tt1,tt2
---
>   integer(kind=8)::n_step
>   integer::ilevel,idim,ivar,info,tot_pt
>   real(kind=8)::tt1,tt2,muspt,muspt_this_step
21c22
<   tt1=MPI_WTIME(info)
---
>   tt1=MPI_WTIME()
35,36c36,37
<   if(mechanical_feedback) call init_mechanical
<   if(stellar_winds) call init_stellar_winds 
---
>   !if(mechanical_feedback.eq.2.or.feedback_refine) call init_mechanical
>   if(mechanical_feedback.eq.2) call init_mechanical
52c53,55
<   tt2=MPI_WTIME(info)
---
>   muspt=0.
>   tot_pt=-1
>   tt2=MPI_WTIME()
67a71
>                                call timer('coarse levels','start')
70c74
<      tt1=MPI_WTIME(info)
---
>      tt1=MPI_WTIME()
122a127
>                                call timer('coarse levels','start')
172c177
<      tt2=MPI_WTIME(info)
---
>      tt2=MPI_WTIME()
177c182,192
<            write(*,*)'Time elapsed since last coarse step:',tt2-tt1
---
>            if (tot_pt==0) muspt=0. ! dont count first timestep
>            n_step = int(numbtot(1,levelmin),kind=8)*twotondim
>            do ilevel=levelmin+1,nlevelmax
>              n_step = n_step + int(numbtot(1,ilevel),kind=8)*product(nsubcycle(levelmin:ilevel-1))*(twotondim-1)
>            enddo
>            muspt_this_step = (tt2-tt1)*1e6/n_step*ncpu
>            muspt = muspt + muspt_this_step
>            tot_pt = tot_pt + 1
>            write(*,'(a,f8.2,a,f12.2,a,f12.2,a)')' Time elapsed since last coarse step:',tt2-tt1 &
>           ,' s',muspt_this_step,' mus/pt'  &
>           ,muspt / max(tot_pt,1), ' mus/pt (av)'
