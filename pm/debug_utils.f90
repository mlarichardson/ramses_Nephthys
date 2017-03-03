!################################################################
!################################################################
!################################################################
subroutine check_igrid(ilevel,identifier)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine computes the thermal energy, the kinetic energy and 
  ! the metal mass dumped in the gas by stars (winds).
  ! This routine is called every fine time step.
  ! NOTICE:
  ! 1) double-check that only stars have tp.ne.0
  ! 2) set use_initial_mass=.true.
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,npart2,icpu,identifier,ierr
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  logical::error
  real(dp)::dx,scale
  integer::nx_loc,icell_out,igrid_out,ilevel2

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  if(ndim.ne.3) return

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  error=.false. 

  ! Gather star particles only.
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        
        ! Count star particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(tp(ipart).ne.0)then ! <-- will remove dark matter and sink particles 
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        
        ! Gather star particles
        if(npart2>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only star particles
              if(tp(ipart).ne.0)then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig  


                 call get_icell_from_pos (xp(ipart,:)/scale, ilevel, igrid_out, icell_out,ilevel2)
                 if(igrid.ne.igrid_out) then
                    error=.true.
                    print *, 'ERROR: ipart,ilevel,myid = ', ipart,ilevel,myid
                    print *, 'ERROR: igrid,igrid_out = ', igrid,igrid_out
                    print *, 'ERROR: xp = ', xp(ipart,:)/scale
                 endif
              endif
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do

  end do 
  ! End loop over cpus

  if(error)then
     call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
  endif

111 format('   Entering stlelar winds for level ',I2)

end subroutine check_igrid
!################################################################
!################################################################
!################################################################
!################################################################
