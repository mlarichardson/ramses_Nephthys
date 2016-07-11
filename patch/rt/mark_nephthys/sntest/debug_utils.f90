!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine check_gas_conservation(gprop)
   use amr_commons
   use hydro_commons
   implicit none
#ifndef WITHOUTMPI
   include 'mpif.h'
   real(dp),dimension(1:ncpu,1:5)::gprop_tot,gprop_mpi
#endif
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   real(dp)::scale,dx,dx_loc,vol_loc
   real(dp)::d,u,v,w
   real(dp)::gprop(1:5)
   integer::nx_loc,ilevel,i,igrid,ind,ngrid,iskip,info,ncache
   integer,dimension(1:nvector),save::ind_cell,ind_grid
   logical,dimension(1:nvector)::ok

   ! Mesh spacing in that level
   nx_loc=(icoarse_max-icoarse_min+1)
   scale=boxlen/dble(nx_loc)

   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   gprop=0d0 

   do ilevel=levelmin,nlevelmax
      dx=0.5D0**ilevel
      dx_loc=dx*scale
      vol_loc=dx_loc**ndim
      ! Loop over grids
      ncache=active(ilevel)%ngrid
      do igrid=1,ncache,nvector
         ngrid=MIN(nvector,ncache-igrid+1)
         do i=1,ngrid
            ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
         end do
   
         ! Loop over cells
         do ind=1,twotondim  
            iskip=ncoarse+(ind-1)*ngridmax
            do i=1,ngrid
               ind_cell(i)=iskip+ind_grid(i)
            end do
   
            ! Flag leaf cells
            do i=1,ngrid
               ok(i)=son(ind_cell(i))==0
            end do
   
            do i=1,ngrid
               if(ok(i))then
                  d=uold(ind_cell(i),1)
                  if(d>0)then
                     u=uold(ind_cell(i),2)/d
                     v=uold(ind_cell(i),3)/d
                     w=uold(ind_cell(i),4)/d
                     gprop(1)=gprop(1)+d*vol_loc
                     gprop(2)=gprop(2)+uold(ind_cell(i),2)*vol_loc
                     gprop(3)=gprop(3)+uold(ind_cell(i),3)*vol_loc
                     gprop(4)=gprop(4)+uold(ind_cell(i),4)*vol_loc
                     gprop(5)=gprop(5)+(uold(ind_cell(i),5)-0.5d0*d*vol_loc*(u*u+v*v+w*w))
!if (d>1.31) print *, myid, ind_cell(i), d
                  endif
               endif
            enddo
   
         enddo ! Loop over cells

      enddo ! Loop over grids
   enddo ! Loop over levels

#ifndef WITHOUTMPI
   gprop_tot=0d0;gprop_mpi=0d0
   do i=1,5
     gprop_tot(myid,i)=gprop(i)
   enddo
   call MPI_ALLREDUCE(gprop_tot,gprop_mpi,ncpu*5,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
   do i=1,5
     gprop(i)=sum(gprop_mpi(:,i))
   enddo
#endif

end subroutine check_gas_conservation 
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine check_star_conservation(sprop)
   use amr_commons
   use pm_commons
   implicit none
#ifndef WITHOUTMPI
   include 'mpif.h'
   real(dp),dimension(1:ncpu,1:4)::sprop_tot,sprop_mpi
#endif
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   real(dp)::scale,dx,dx_loc,vol_loc
   real(dp)::d,u,v,w
   real(dp)::sprop(1:4)
   integer::nx_loc,icpu,jgrid,jpart,ipart,igrid,npart1,next_part,info,i
   integer,dimension(1:nvector),save::ind_cell,ind_grid
   logical,dimension(1:nvector)::ok


   ! Mesh spacing in that level
   nx_loc=(icoarse_max-icoarse_min+1)
   scale=boxlen/dble(nx_loc)

   ! Conversion factor from user units to cgs units
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   sprop=0
   ! Count the number of star particles that can *potentially* undergo the blast phase
   do icpu=1,ncpu
   ! Loop over cpus
       igrid=headl(icpu,levelmin)
       ! Loop over grids
       do jgrid=1,numbl(icpu,levelmin)
          npart1=numbp(igrid)  ! Number of particles in the grid
          ! Count star particles younger than t_ctw
          if(npart1>0)then
             ipart=headp(igrid)
             ! Loop over particles
             do jpart=1,npart1
                ! Save next particle   <--- Very important !!!
                next_part=nextp(ipart)
                if( tp(ipart).ne.0) then
                   sprop(1)=sprop(1)+mp(ipart)
                   sprop(2)=sprop(2)+mp(ipart)*vp(ipart,1)
                   sprop(3)=sprop(3)+mp(ipart)*vp(ipart,2)
                   sprop(4)=sprop(4)+mp(ipart)*vp(ipart,3)
                endif
                ipart=next_part  ! Go to next particle
             end do
          endif
          igrid=next(igrid)   ! Go to next grid
      end do
   enddo

#ifndef WITHOUTMPI
   sprop_tot=0d0;sprop_mpi=0d0
   do i=1,4
     sprop_tot(myid,i)=sprop(i)
   enddo
   call MPI_ALLREDUCE(sprop_tot,sprop_mpi,ncpu*4,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
   do i=1,4
     sprop(i)=sum(sprop_mpi(:,i))
   enddo
#endif

end subroutine check_star_conservation
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine check_gas_status(identifier)
   use amr_commons
   use hydro_commons
   implicit none
#ifndef WITHOUTMPI
   include 'mpif.h'
#endif
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   real(dp)::scale,dx,dx_loc,vol_loc
   real(dp)::d,u,v,w,eth,ekk,Tk,u2,Z
   real(dp)::gprop(1:5)
   integer::nx_loc,ilevel,i,igrid,ind,ngrid,iskip,info,ncache,identifier
   integer,dimension(1:nvector),save::ind_cell,ind_grid
   logical,dimension(1:nvector)::ok
   logical::nokay
   ! Mesh spacing in that level
   nx_loc=(icoarse_max-icoarse_min+1)
   scale=boxlen/dble(nx_loc)

   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   gprop=0d0 

   do ilevel=levelmin,nlevelmax
      dx=0.5D0**ilevel
      dx_loc=dx*scale
      vol_loc=dx_loc**ndim
      ! Loop over grids
      ncache=active(ilevel)%ngrid
      do igrid=1,ncache,nvector
         ngrid=MIN(nvector,ncache-igrid+1)
         do i=1,ngrid
            ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
         end do
   
         ! Loop over cells
         do ind=1,twotondim  
            iskip=ncoarse+(ind-1)*ngridmax
            do i=1,ngrid
               ind_cell(i)=iskip+ind_grid(i)
            end do
   
            ! Flag leaf cells
            do i=1,ngrid
               ok(i)=son(ind_cell(i))==0
            end do
   
            do i=1,ngrid
               if(ok(i))then
                  d=uold(ind_cell(i),1)
                  u=uold(ind_cell(i),2)/d
                  v=uold(ind_cell(i),3)/d
                  w=uold(ind_cell(i),4)/d
                  ekk=d*(u*u+v*v+w*w)/2d0
                  eth=uold(ind_cell(i),5)-ekk
                  Tk =eth*scale_T2*(gamma-1)*0.6/d
                  u2=sqrt(u*u+v*v+w*w)*scale_v/1d5

                  nokay=.false.
                  if(d*scale_nH>1d4)then
                     write(*,100) log10(d*scale_nH),log10(Tk),myid,ind_cell(i),identifier
                     nokay=.true.
                  endif
                  !if(nokay) call clean_stop
               endif
            enddo
   
         enddo ! Loop over cells

      enddo ! Loop over grids
   enddo ! Loop over levels

100    format(" Den Err log|nH|=",F5.1,F5.1," myid=",I3," icell=",I8," info=",I4)
101    format(" Vel Err log|v|=",f5.1," myid=",I3," icell=",I8," info=",I4, 3(1x,f5.1))
102    format(" Met Err log|Z|=",f5.1," Z =",e10.3," log nH=",f5.1," myid=",I3," icell=",I8," info=",I4)
103    format(" Tem Err log|T|=",f5.1," log T=",f5.1," log nH=",f5.1," myid=",I3," icell=",I8," info=",I2," ilevel=",I2)
end subroutine check_gas_status
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine check_gas_status_new(identifier)
   use amr_commons
   use hydro_commons
   implicit none
#ifndef WITHOUTMPI
   include 'mpif.h'
#endif
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   real(dp)::scale,dx,dx_loc,vol_loc
   real(dp)::d,u,v,w,eth,ekk,Tk,u2,Z
   real(dp)::gprop(1:5)
   integer::nx_loc,ilevel,i,igrid,ind,ngrid,iskip,info,ncache,identifier
   integer,dimension(1:nvector),save::ind_cell,ind_grid
   logical,dimension(1:nvector)::ok
   logical::nokay

   ! Mesh spacing in that level
   nx_loc=(icoarse_max-icoarse_min+1)
   scale=boxlen/dble(nx_loc)

   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   gprop=0d0 

   do ilevel=levelmin,nlevelmax
      dx=0.5D0**ilevel
      dx_loc=dx*scale
      vol_loc=dx_loc**ndim
      ! Loop over grids
      ncache=active(ilevel)%ngrid
      do igrid=1,ncache,nvector
         ngrid=MIN(nvector,ncache-igrid+1)
         do i=1,ngrid
            ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
         end do
   
         ! Loop over cells
         do ind=1,twotondim  
            iskip=ncoarse+(ind-1)*ngridmax
            do i=1,ngrid
               ind_cell(i)=iskip+ind_grid(i)
            end do
   
            ! Flag leaf cells
            do i=1,ngrid
               ok(i)=son(ind_cell(i))==0
            end do
   
            do i=1,ngrid
               if(ok(i))then
                  d=unew(ind_cell(i),1)
                  u=unew(ind_cell(i),2)/d
                  v=unew(ind_cell(i),3)/d
                  w=unew(ind_cell(i),4)/d
                  ekk=d*(u*u+v*v+w*w)/2d0
                  eth=unew(ind_cell(i),5)-ekk
                  Tk =eth*scale_T2*(gamma-1)*0.6/d
                  u2=sqrt(u*u+v*v+w*w)*scale_v/1d5
                  nokay=.false.
                  if(d.ne.d.or.d<1d-10)then
                     write(*,100) log10(abs(d)*scale_nH),myid,ind_cell(i),identifier
                     nokay=.true.
                  endif
                  if(u2>1d4.or.u2.ne.u2)then
                     write(*,101) log10(u2), myid, ind_cell(i), identifier,&
                   & log10(abs(u)),log10(abs(v)),log10(abs(w))
                     nokay=.true.
                  endif
                  if(metal)then
                     Z=unew(ind_cell(i),imetal)/d
                     if(Z<0.or.Z.ne.Z)then
                        write(*,102) log10(abs(Z)),Z,log10(d*scale_nH),myid,ind_cell(i),identifier
                        print *, 'imetal =',imetal
                        nokay=.true.
                     endif
                  endif
                  if(Tk>1d10.or.Tk<0.or.Tk.ne.Tk)then
                     write(*,103) log10(abs(Tk)), log10(Tk), log10(d*scale_nH), myid, &
                    & ind_cell(i), identifier, ilevel
                     nokay=.true.
                  endif

                  !if(nokay) call clean_stop
               endif
            enddo
   
         enddo ! Loop over cells

      enddo ! Loop over grids
   enddo ! Loop over levels

100    format(" Den Err log|nH|=",F5.1," myid=",I3," icell=",I8," info=",I4)
101    format(" Vel Err log|v|=",f5.1," myid=",I3," icell=",I8," info=",I4, 3(1x,f5.1))
102    format(" Met Err log|Z|=",f5.1," Z =",e10.3," log nH=",f5.1," myid=",I3," icell=",I8," info=",I4)
103    format(" Tem Err log|T|=",f5.1," log T=",f5.1," log nH=",f5.1," myid=",I3," icell=",I8," info=",I2," ilevel=",I2)
end subroutine check_gas_status_new
!####################################################################
!####################################################################
!####################################################################
!####################################################################
subroutine check_gas_max_temperature(identifier)
   use amr_commons
   use hydro_commons
   implicit none
#ifndef WITHOUTMPI
   include 'mpif.h'
#endif
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   real(dp)::scale,dx,dx_loc,vol_loc
   real(dp)::d,u,v,w,Tk,Tkmax,ekk,eth
   integer::nx_loc,ilevel,i,igrid,ind,ngrid,iskip,info,ncache,identifier
   integer,dimension(1:nvector),save::ind_cell,ind_grid
   logical,dimension(1:nvector)::ok

   ! Mesh spacing in that level
   nx_loc=(icoarse_max-icoarse_min+1)
   scale=boxlen/dble(nx_loc)

   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   do ilevel=levelmin,nlevelmax
      dx=0.5D0**ilevel
      dx_loc=dx*scale
      vol_loc=dx_loc**ndim
      ! Loop over grids
      ncache=active(ilevel)%ngrid
      do igrid=1,ncache,nvector
         ngrid=MIN(nvector,ncache-igrid+1)
         do i=1,ngrid
            ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
         end do
   
         ! Loop over cells
         do ind=1,twotondim  
            iskip=ncoarse+(ind-1)*ngridmax
            do i=1,ngrid
               ind_cell(i)=iskip+ind_grid(i)
            end do
   
            ! Flag leaf cells
            do i=1,ngrid
               ok(i)=son(ind_cell(i))==0
            end do
   
            do i=1,ngrid
               if(ok(i))then
                  d=uold(ind_cell(i),1)
                  u=uold(ind_cell(i),2)/d
                  v=uold(ind_cell(i),3)/d
                  w=uold(ind_cell(i),4)/d
                  ekk=d*(u*u+v*v+w*w)/2d0
                  eth=uold(ind_cell(i),5)-ekk
                  Tk =eth*scale_T2*(gamma-1)*0.6/d
                  if(Tk>Tkmax) Tkmax=Tk 
               endif
            enddo
   
         enddo ! Loop over cells

      enddo ! Loop over grids
   enddo ! Loop over levels

   if(Tkmax>1d7)then
      write(*,*)'TK',myid,Tkmax
   endif

end subroutine check_gas_max_temperature

