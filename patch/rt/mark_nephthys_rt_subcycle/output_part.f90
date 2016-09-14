! patch changes:
! - writing of mp0
! Joki added st_n_tp and st_n_sn variables (oct 2013)
! - writing of chp (Taysun, Jun 2016)
subroutine backup_part(filename)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(LEN=80)::filename

  integer::i,idim,ilun,ipart,ich
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::ii
  integer(i8b),allocatable,dimension(:)::ii8
  integer,allocatable,dimension(:)::ll
  logical,allocatable,dimension(:)::nb
  integer,parameter::tag=1122
  integer::dummy_io,info2
  
  if(verbose)write(*,*)'Entering backup_part'

  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif

  
  ilun=2*ncpu+myid+10

  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)
  open(unit=ilun,file=TRIM(fileloc),form='unformatted')
  rewind(ilun)
  ! Write header
  write(ilun)ncpu
  write(ilun)ndim
  write(ilun)npart
  write(ilun)localseed
  write(ilun)nstar_tot   
  write(ilun)mstar_tot   
  write(ilun)mstar_lost
  write(ilun)nsink
  ! Write position
  allocate(xdp(1:npart))
  do idim=1,ndim
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           xdp(ipart)=xp(i,idim)
         end if
     end do
     write(ilun)xdp
  end do
  ! Write velocity
  do idim=1,ndim
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           xdp(ipart)=vp(i,idim)
        end if
     end do
     write(ilun)xdp
  end do
  ! Write mass
  ipart=0
  do i=1,npartmax
     if(levelp(i)>0)then
        ipart=ipart+1
        xdp(ipart)=mp(i)
     end if
  end do
  write(ilun)xdp
  deallocate(xdp)
  ! Write identity
  allocate(ii8(1:npart))
  ipart=0
  do i=1,npartmax
     if(levelp(i)>0)then
        ipart=ipart+1
        ii8(ipart)=idp(i)
     end if
  end do
  write(ilun)ii8
  deallocate(ii8)
  ! Write level
  allocate(ll(1:npart))
  ipart=0
  do i=1,npartmax
     if(levelp(i)>0)then
        ipart=ipart+1
        ll(ipart)=levelp(i)
     end if
  end do
  write(ilun)ll
  deallocate(ll)

#ifdef OUTPUT_PARTICLE_POTENTIAL
  ! Write potential (added by AP)
  allocate(xdp(1:npart))
  ipart=0
  do i=1, npartmax
     if(levelp(i)>0) then
        ipart=ipart+1
        xdp(ipart)=ptcl_phi(i)
     end if
  end do
  write(ilun)xdp
  deallocate(xdp)
#endif

  ! Write birth epoch
  if(star.or.sink)then
     allocate(xdp(1:npart))
     ipart=0
     do i=1,npartmax
        if(levelp(i)>0)then
           ipart=ipart+1
           xdp(ipart)=tp(i)
        end if
     end do
     write(ilun)xdp
     if(write_stellar_densities) then
        ! Write gas density at birth
        ipart=0
        do i=1,npartmax
           if(levelp(i)>0)then
              ipart=ipart+1
              xdp(ipart)=st_n_tp(i)
           end if
        end do
        write(ilun)xdp
        ! Write gas density at SN
        ipart=0
        do i=1,npartmax
           if(levelp(i)>0)then
              ipart=ipart+1
              xdp(ipart)=st_n_SN(i)
           end if
        end do
        write(ilun)xdp
        ! Write SN energy injected
        ipart=0
        do i=1,npartmax
           if(levelp(i)>0)then
              ipart=ipart+1
              xdp(ipart)=st_e_SN(i)
           end if
        end do
        write(ilun)xdp
     endif ! if write_stellar_densities

     ! Write metallicity
     if(metal)then
        ipart=0
        do i=1,npartmax
           if(levelp(i)>0)then
              ipart=ipart+1
              xdp(ipart)=zp(i)
           end if
        end do
        write(ilun)xdp
     end if

     ! Write initial mass
     if(use_initial_mass)then
        ipart=0
        do i=1,npartmax
           if(levelp(i)>0)then
              ipart=ipart+1
              xdp(ipart)=mp0(i)
           end if
        end do
        write(ilun)xdp
     endif

     ! Write chemical abundance
     do ich=1,nchem
        ipart=0
        do i=1,npartmax
           if(levelp(i)>0)then
              ipart=ipart+1
              xdp(ipart)=chp(i,ich)
           end if
        end do
        write(ilun)xdp
     end do


     deallocate(xdp)
  end if
  close(ilun)

  ! Send the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
             & MPI_COMM_WORLD,info2)
     end if
  endif
#endif


end subroutine backup_part

