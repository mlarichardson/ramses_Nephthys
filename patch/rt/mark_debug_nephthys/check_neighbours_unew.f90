!################################################################
!################################################################
subroutine check_neighbours_unew_fine(ilevel)
  use amr_commons
  use pm_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine sends of a few grids and checks the value of unew for it's neighbours
  ! 
  ! This routine is called every fine time step.
  ! NOTICE:
  ! 1) will call on many grids so keep to small runs.
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,npart2,icpu
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  ! Gather star particles only.

#if NDIM==3
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
       if (son(igrid) == 0) then
         ig=ig+1
         ind_grid(ig)=igrid

      
         if(ig==min(6,nvector/5) .or. jgrid == numbl(icpu,ilevel))then
           ip = 5*ig
           call check_neighbours_unew(ind_grid,ig,ip,ilevel)
           ig=0
         end if

         if(ig==0)then
           exit
         endif
         igrid=next(igrid)   ! Go to next grid

       endif
    enddo

    if(ig > 0) then 
       ip = 5*ig
       call check_neighbours_unew(ind_grid,ig,ip,ilevel)
    end if

  end do 
  ! End loop over cpus

#endif

111 format('   Entering check neighbour fine for level ',I2)

end subroutine check_neighbours_unew_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine check_neighbours_unew(ind_grid,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine check neighbour unew_fine.
  ! It checks a max of 10 grids, each with 5 fake particles.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc, ich
  real(dp)::dx_min,vol_min
  real(dp)::dx,dx_loc,scale,birth_time,d
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mloss,mzloss,ethermal,ekinetic,dteff
  real(dp),dimension(1:nvector),save::vol_loc
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp)::msun2g=2d33
  ! stellar library
  real(dp)::mejecta,dfmloss,dfmzloss,log_deloss_erg
  real(dp)::mstar_ini,mstar_ini_msun,zstar
  real(dp)::unit_e_code
  real(dp)::dfmloss_spec(1:nchem)
  real(dp),dimension(1:nchem,1:nvector)::mloss_spec
  ! fractional abundances ; for ionisation fraction and ref, etc 
  real(dp),dimension(1:nvector,1:NVAR),save::fractions ! not compatible with delayed cooling
  integer::ivar,i_fractions,parent_grid

  logical::neighbour_missing=.false.

  ! starting index for passive variables except for imetal and chem
  i_fractions = imetal+nchem+1  

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc(1:nvector)=dx_loc**ndim
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

#if NDIM==3
  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale position at level ilevel - stellar winds has x on 0.5, 5.5
  ! kg = 1, 12, 14, 17, 24
  !   (1.2, 1.2, 1.2), (4.8,1.2,2.4), (2.4,3.6,3.6), (3.6,4.8,2.7), (4.4,3.3,4.8)
  x(1,1) = 1.2d0 ; x(1,2) = 1.2d0 ; x(1,3) = 1.2d0
  x(2,1) = 4.8d0 ; x(2,2) = 1.2d0 ; x(2,3) = 2.4d0
  x(3,1) = 2.4d0 ; x(3,2) = 3.6d0 ; x(3,3) = 3.6d0
  x(4,1) = 3.6d0 ; x(4,2) = 4.8d0 ; x(4,3) = 2.7d0
  x(5,1) = 4.4d0 ; x(5,2) = 3.3d0 ; x(5,3) = 4.8d0

  do idim=1,ndim
    do j=6,np
      x(j,idim) = x(j-5,idim) 
    enddo
  enddo

  do j=1,np
    ind_grid_part(j) = (j-1)/5 + 1
  enddo

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=int(x(j,idim))
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igd(j,idim)=id(j,idim)/2
     end do
  end do
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
  do j=1,np
     igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end if
     end do
  end do
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do

  ! Compute parent cell adresses
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     else
        indp(j) = nbors_father_cells(ind_grid_part(j),kg(j))
        vol_loc(j)=vol_loc(j)*2**ndim ! ilevel-1 cell volume
     end if
  end do

  do j=1,np
    parent_grid = nbors_father_cells(ind_grid_part(j),kg(j)) - (nbors_father_cells(ind_grid_part(j),kg(j))/ngridmax)*ngridmax - 1

    write(6,'(a,3i)') "CheckNGBR A: ", myid, j, ind_grid(ind_grid_part(j))
    write(6,'(a,i,3es14.6,i)') "CheckNGBR B: ", myid, x(j,1), x(j,2), x(j,3), kg(j)
    write(6,'(a,4i)') "CheckNGBR C: ", myid, nbors_father_cells(ind_grid_part(j),kg(j)), parent_grid, cpu_map(parent_grid)
    write(6,'(a,4i,es14.6)') "CheckNGBR D: ", myid, igrid(j), icell(j), indp(j), unew(indp(j),1)

    if ( unew(indp(j),1) .eq. 0. ) neighbour_missing = .true.
  enddo

  if (neighbour_missing) stop
 
#endif
  
end subroutine check_neighbours_unew
!################################################################
!################################################################
