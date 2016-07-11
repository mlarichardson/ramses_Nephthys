! Added  st_n_tp and st_n_sn variables (oct 2013)
! Also T threshold for star formation
!################################################################
!################################################################
!################################################################
!################################################################
subroutine star_formation(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use cooling_module, ONLY: XH=>X,rhoc,mH,twopi,T2_min_fix,mu_mol,kB 
  use random
#ifdef ATON
  use radiation_commons, ONLY: Erad
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Description: This subroutine spawns star-particle of constant mass
  ! using a Poissdiverging ity law if some gas condition are fulfilled. 
  ! It modifies hydrodynamic variables according to mass conservation 
  ! and assumes an isothermal transformation... 
  ! On exit, the gas velocity and sound speed are unchanged.
  ! New star particles are synchronized with other collisionless particles.
  ! Array flag2 is used as temporary work space.
  ! Yann Rasera  10/2002-01/2003
  !----------------------------------------------------------------------
  ! local constants
  real(dp) :: t0,d0,e0,mgas,mcell
  real(dp) :: scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,kms
  real(dp),dimension(1:twotondim,1:3) :: xc
  ! other variables
  integer  :: ncache,nnew,ivar,ngrid,icpu,index_star,ilun,irad
  integer  :: igrid,ix,iy,iz,ind,i,j,n,iskip,istar,inew,nx_loc
  integer  :: ntot,ntot_all,info,nstar_corrected,iche,ncell
  logical  :: ok_free,ok_all
  real(dp) :: d,x,y,z,u,v,w,e,zg,vdisp,dgas,uavg,vavg,wavg,dtot,tcell
  real(dp) :: mstar,dstar,tstar,nISM,nCOM
  real(dp) :: velc,uc,vc,wc,mass_load,ul,vl,wl,ur,vr,wr,divv,curlv,alpha
  real(dp) :: vxgauss,vygauss,vzgauss,birth_epoch,factG
  real(kind=8) :: mlost,mtot,mlost_all,mtot_all
  real(kind=8) :: RandNum,GaussNum,PoissMean,RandVel(1:3) 
  real(dp),dimension(1:3)::skip_loc
  real(dp) :: dx,dx_loc,scale,vol_loc,dx_min,vol_min,d1,d2,d3,d4,d5,d6
  real(dp) :: bx1,bx2,by1,by2,bz1,bz2
  real(dp) :: T2_EOS,T2,T2_old,polytropic_constant,boost,Zsolar,t_cool
  real(dp) :: t_ff,t_dyn,alpha0,trgv,temp,scrit,phi_t,e_cts
  real(dp) :: d_gmc,d_dc,sigs,c_s2,theta,lamjt
  real(dp) :: onepi,mcell1,mjeans1,sound_speed
  real(dp) :: mstar_nsn,scale_msun,nH,T2min
  real(dp),dimension(1:nvector)::effArr,vturb ! for Padoan method
  integer ,dimension(1:ncpu,1:IRandNumSize)    :: allseed
  integer ,dimension(1:nvector),save           :: ind_grid,ind_cell,nstar
  integer ,dimension(1:nvector),save           :: ind_grid_new,ind_cell_new,ind_part
  logical ,dimension(1:nvector),save           :: ok,ok_new=.true.,ok_true=.true.
  integer ,dimension(1:ncpu)                   :: ntot_star_cpu,ntot_star_all
  integer ,dimension(1:nvector,0:twondim)      :: ind_nbor
  real(dp), dimension(1:nvector)               :: sfr_ff
#ifdef NCHEM
  real(dp),dimension(1:nchem) :: chem1=0.0
#endif
  logical::ok_exist

  if (numbtot(1,ilevel)==0) return
  if (.not. hydro) return
  if (ndim.ne.3) return


  if (verbose) write(*,*)' Entering star_formation'
  

  ilun = 3*ncpu + myid + 10 
 
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  kms = scale_v/1d5
  scale_msun = scale_l**3*scale_d/1.989d33

  ! Mesh spacing in that level
  dx       = 0.5D0**ilevel 
  nx_loc   = (icoarse_max-icoarse_min+1)
  skip_loc = (/0.0d0,0.0d0,0.0d0/)
  if (ndim>0) skip_loc(1)=dble(icoarse_min)
  if (ndim>1) skip_loc(2)=dble(jcoarse_min)
  if (ndim>2) skip_loc(3)=dble(kcoarse_min)
  scale    = boxlen/dble(nx_loc)
  dx_loc   = dx*scale
  vol_loc  = dx_loc**ndim
  dx_min   = (0.5D0**nlevelmax)*scale
  vol_min  = dx_min**ndim

  ! Star formation time scale from Gyr to code units 
  t0    = t_star*(1d9*365.*24.*3600.)/scale_t
  ! ISM density threshold from H/cc to code units    
  nISM  = n_star
  nCOM  = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*XH/mH
  nISM  = MAX(nCOM,nISM)
  d0    = nISM/scale_nH

  !------------------------------------------------------
  ! Set the star particle mass from the number of SN [TK]
  !------------------------------------------------------
  if(nsn2mass>0.and.fstar_min<0)then
     ! Mesh spacing
     mstar_nsn  = (nsn2mass*M_SNII)/eta_sn/scale_msun
     ! ISM density threshold from H/cc to code units    
     mstar      = n_star/(scale_nH*aexp**3)*vol_min
     fstar_min  = mstar_nsn/mstar
     if(myid==1) write(*,*) ">>>TKNOTE: Mstar,min=",mstar_nsn*scale_msun,fstar_min
  endif

  mstar = n_star/(scale_nH*aexp**3)*vol_min*fstar_min
  dstar = mstar/vol_loc
  factG = 1d0
  if (cosmo) factG = 3d0/4d0/twopi*omega_m*aexp
  ! formation time of the star particle
  if(use_proper_time)then
     birth_epoch = texp
  else
     birth_epoch = t
  endif
  d_gmc = MAX(nCOM,n_gmc)/scale_nH
  d_dc  = MAX(nCOM,n_dc)/scale_nH    ! dark cloud
  onepi = twopi/2d0

  ! Cells center position relative to grid center position
  do ind=1,twotondim  
     iz        = (ind-1)/4
     iy        = (ind-1-4*iz)/2
     ix        = (ind-1-2*iy-4*iz)
     xc(ind,1) = (dble(ix)-0.5D0)*dx
     xc(ind,2) = (dble(iy)-0.5D0)*dx
     xc(ind,3) = (dble(iz)-0.5D0)*dx
  end do

  ! If necessary, initialize random number generator
  if (localseed(1)==-1) then
     call rans(ncpu,iseed,allseed)
     localseed = allseed(myid,1:IRandNumSize)
  end if

#if NDIM==3
  !------------------------------------------------
  ! Convert hydro variables to primitive variables
  !------------------------------------------------
  ncache = active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid = MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i) = active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim  
        iskip = ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i) = iskip+ind_grid(i)
        end do
        do i=1,ngrid
           d   = uold(ind_cell(i),1)
           u   = uold(ind_cell(i),2)/d
           v   = uold(ind_cell(i),3)/d
           w   = uold(ind_cell(i),4)/d
           e   = uold(ind_cell(i),5)
#ifdef SOLVERmhd
           bx1 = uold(ind_cell(i),6)
           by1 = uold(ind_cell(i),7)
           bz1 = uold(ind_cell(i),8)
           bx2 = uold(ind_cell(i),nvar+1)
           by2 = uold(ind_cell(i),nvar+2)
           bz2 = uold(ind_cell(i),nvar+3)
           e   = e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
           e   = e-0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=1,nener
              e=e-uold(ind_cell(i),5+irad)
           end do
#endif
           uold(ind_cell(i),1) = d
           uold(ind_cell(i),2) = u
           uold(ind_cell(i),3) = v
           uold(ind_cell(i),4) = w
           uold(ind_cell(i),5) = e/d
        end do
        do ivar=imetal,nvar
           do i=1,ngrid
              d = uold(ind_cell(i),1)
              w = uold(ind_cell(i),ivar)/d
              uold(ind_cell(i),ivar) = w
           end do
        end do
     end do
  end do

! get values of uold for density and velocities in virtual boundaries
#ifndef WITHOUTMPI
  do ivar=1,4
     call make_virtual_fine_dp(uold(1,ivar),ilevel)
  end do
#endif

  !------------------------------------------------
  ! Compute number of new stars in each cell
  !------------------------------------------------
  ntot = 0
  ! Loop over grids
  ncache = active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid = MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i) = active(ilevel)%igrid(igrid+i-1)
     end do
     ! Star formation criterion ---> logical array ok(i)
     do ind=1,twotondim
        iskip = ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i) = iskip+ind_grid(i)
        end do
        ! Flag leaf cells
        do i=1,ngrid
           ok(i) = (son(ind_cell(i))==0)
        end do


        if (TRIM(star_maker)=='density') then
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              ! Density criterion
              if(d<=d0)ok(i)=.false.
              ! Temperature criterion
              if(ok(i))then
                 nH    = d*scale_nH
                 T2min = T2_star*(nH/nISM)**(g_star-1.0)
                 !notice that uold is not primitive variable at this point
                 T2    = uold(ind_cell(i),5) 
                 T2    = (gamma-1)*T2*scale_T2 
                 T2    = max(T2-T2min,T2_min_fix)
                 if(T2.gt.T2thres_SF) ok(i)=.false.
              endif
           end do          
            
        else if (TRIM(star_maker)=='denmax') then
           ! Density criterion
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              if(d<=d0)then
                 ok(i)=.false.
              else
                 ! we need local density to determine if the cell we are interested 
                 ! is the local density maxima
                 ncell = 1 ! we just want the neighbors of that cell
                 call getnbor(ind_cell(i),ind_nbor,ncell,ilevel)
                 ! First calculate velocities for divergence
                 d1    = uold(ind_nbor(1,1),1)     ; d4    = uold(ind_nbor(1,4),1) 
                 d2    = uold(ind_nbor(1,2),1)     ; d5    = uold(ind_nbor(1,5),1) 
                 d3    = uold(ind_nbor(1,3),1)     ; d6    = uold(ind_nbor(1,6),1)
                 if(d<d1)ok(i)=.false.
                 if(d<d2)ok(i)=.false.
                 if(d<d3)ok(i)=.false.
                 if(d<d4)ok(i)=.false.
                 if(d<d5)ok(i)=.false.
                 if(d<d6)ok(i)=.false.
              endif
           end do           
        else if(TRIM(star_maker)=='hopkins')then
           ! Enforce self-gravitating criterion a la Hopkins et al 2013
           do i=1,ngrid
              ! if cell is a leaf cell
              if (ok(i)) then 
                 d = uold(ind_cell(i),1)
                 ! density threshold to avoid computing the criterion too often but otherwise not needed
                 if (d <= d_gmc) then
                    ok(i) = .false.
                 else
                    ! we need velocity divergence and curl estimates in the cell, so construct values of velocity field 
                    ! on the 6 faces of the cell using simple linear interpolation from neighbouring cell values and differentiate.
                    ! get neighbor cells if they exist, otherwise use straight injection from local cell
                    ncell = 1 ! we just want the neighbors of that cell
                    call getnbor(ind_cell(i),ind_nbor,ncell,ilevel)
                    ! First calculate velocities for divergence
                    d1    = uold(ind_nbor(1,1),1)     ; d4    = uold(ind_nbor(1,4),1) 
                    d2    = uold(ind_nbor(1,2),1)     ; d5    = uold(ind_nbor(1,5),1) 
                    d3    = uold(ind_nbor(1,3),1)     ; d6    = uold(ind_nbor(1,6),1)  
                    ul    = (d2*uold(ind_nbor(1,2),2) + d*uold(ind_cell(i),2))/(d2+d)
                    vl    = (d4*uold(ind_nbor(1,4),3) + d*uold(ind_cell(i),3))/(d4+d)
                    wl    = (d6*uold(ind_nbor(1,6),4) + d*uold(ind_cell(i),4))/(d6+d)
                    ur    = (d1*uold(ind_nbor(1,1),2) + d*uold(ind_cell(i),2))/(d1+d)
                    vr    = (d3*uold(ind_nbor(1,3),3) + d*uold(ind_cell(i),3))/(d3+d)
                    wr    = (d5*uold(ind_nbor(1,5),4) + d*uold(ind_cell(i),4))/(d5+d)
                    divv  = ( ((ur-ul)+(vr-vl)+(wr-wl)) / dx_loc )**2 
                    ! Then calculate velocities for curl
                    vl    = (d6*uold(ind_nbor(1,6),3) + d*uold(ind_cell(i),3))/(d6+d)
                    wl    = (d4*uold(ind_nbor(1,4),4) + d*uold(ind_cell(i),4))/(d4+d)
                    vr    = (d5*uold(ind_nbor(1,5),3) + d*uold(ind_cell(i),3))/(d5+d)
                    wr    = (d3*uold(ind_nbor(1,3),4) + d*uold(ind_cell(i),4))/(d3+d)
                    curlv = ((wr-wl)-(vr-vl))**2         ! first term
                    ul    = (d6*uold(ind_nbor(1,6),2) + d*uold(ind_cell(i),2))/(d6+d)
                    wl    = (d2*uold(ind_nbor(1,2),4) + d*uold(ind_cell(i),4))/(d2+d)
                    ur    = (d5*uold(ind_nbor(1,5),2) + d*uold(ind_cell(i),2))/(d5+d)
                    wr    = (d1*uold(ind_nbor(1,1),4) + d*uold(ind_cell(i),4))/(d1+d)
                    curlv = curlv + ((ur-ul)-(wr-wl))**2 ! second term
                    ul    = (d4*uold(ind_nbor(1,4),2) + d*uold(ind_cell(i),2))/(d4+d)
                    vl    = (d2*uold(ind_nbor(1,2),3) + d*uold(ind_cell(i),3))/(d2+d)
                    ur    = (d3*uold(ind_nbor(1,3),2) + d*uold(ind_cell(i),2))/(d3+d)
                    vr    = (d1*uold(ind_nbor(1,1),3) + d*uold(ind_cell(i),3))/(d1+d)
                    curlv = curlv + ((vr-vl)-(ur-ul))**2 ! third term
                    curlv = curlv / dx_loc**2 
                    ! Check if gas in cell is self-gravitating (alpha < 1)
                    alpha = 0.5d0*(divv+curlv)/(factG*d)
                    if (alpha > 1d0) ok(i) = .false. 
                 endif
              endif
           end do

        else if(TRIM(star_maker)=='federrath')then
           ! Enforce turbulence criterion + efficiency following Federrath & Klessen 2012
           do i=1,ngrid
              ! if cell is a leaf cell
              if (ok(i)) then 
                 d = uold(ind_cell(i),1)
                 ! set density threshold to avoid to compute tensor trace too often but otherwise not needed
                 if (d <= d_gmc) then
                    ok(i) = .false.
                 else
                    ! We need to estimate the norm of the gradient of the velocity field in the cell (tensor of 2nd rank)
                    ! i.e. || A ||^2 = trace( A A^T) where A = grad vec(v) is the tensor. 
                    ! So construct values of velocity field on the 6 faces of the cell using simple linear interpolation 
                    ! from neighbouring cell values and differentiate. 
                    ! Get neighbor cells if they exist, otherwise use straight injection from local cell
                    ncell = 1 ! we just want the neighbors of that cell
                    call getnbor(ind_cell(i),ind_nbor,ncell,ilevel)
                    d1    = uold(ind_nbor(1,1),1)     ; d4    = uold(ind_nbor(1,4),1) 
                    d2    = uold(ind_nbor(1,2),1)     ; d5    = uold(ind_nbor(1,5),1) 
                    d3    = uold(ind_nbor(1,3),1)     ; d6    = uold(ind_nbor(1,6),1)  
                    ul    = (d2*uold(ind_nbor(1,2),2) + d*uold(ind_cell(i),2))/(d2+d)
                    ur    = (d1*uold(ind_nbor(1,1),2) + d*uold(ind_cell(i),2))/(d1+d)
                    trgv  = (ur-ul)**2
                    ul    = (d4*uold(ind_nbor(1,4),3) + d*uold(ind_cell(i),3))/(d4+d)
                    ur    = (d3*uold(ind_nbor(1,3),3) + d*uold(ind_cell(i),3))/(d3+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d6*uold(ind_nbor(1,6),4) + d*uold(ind_cell(i),4))/(d6+d)
                    ur    = (d5*uold(ind_nbor(1,5),4) + d*uold(ind_cell(i),4))/(d5+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d6*uold(ind_nbor(1,6),3) + d*uold(ind_cell(i),3))/(d6+d)
                    ur    = (d5*uold(ind_nbor(1,5),3) + d*uold(ind_cell(i),3))/(d5+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d4*uold(ind_nbor(1,4),4) + d*uold(ind_cell(i),4))/(d4+d)
                    ur    = (d3*uold(ind_nbor(1,3),4) + d*uold(ind_cell(i),4))/(d3+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d6*uold(ind_nbor(1,6),2) + d*uold(ind_cell(i),2))/(d6+d)
                    ur    = (d5*uold(ind_nbor(1,5),2) + d*uold(ind_cell(i),2))/(d5+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d2*uold(ind_nbor(1,2),4) + d*uold(ind_cell(i),4))/(d2+d)
                    ur    = (d1*uold(ind_nbor(1,1),4) + d*uold(ind_cell(i),4))/(d1+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d4*uold(ind_nbor(1,4),2) + d*uold(ind_cell(i),2))/(d4+d)
                    ur    = (d3*uold(ind_nbor(1,3),2) + d*uold(ind_cell(i),2))/(d3+d)
                    trgv  = trgv + (ur-ul)**2
                    ul    = (d2*uold(ind_nbor(1,2),3) + d*uold(ind_cell(i),3))/(d2+d)
                    ur    = (d1*uold(ind_nbor(1,1),3) + d*uold(ind_cell(i),3))/(d1+d)
                    trgv  = trgv + (ur-ul)**2
                    trgv  = trgv 
                    ! now compute sound speed squared = cs^2 
                    temp  = uold(ind_cell(i),5)*(gamma -1.0)
                    ! prevent numerical crash due to negative temperature
                    temp  = max(temp,smallc**2)
                    c_s2  = temp    
                    ! Calculate "turbulent" Jeans length in cell units, lamjt 
                    ! (see e.g. Chandrasekhar 51, Bonazzola et al 87, Federrath & Klessen 2012 eq 36)
                    lamjt = (onepi*trgv + sqrt(onepi*onepi*trgv*trgv + 36.0*onepi*c_s2*factG*d*dx_loc**2))/(6.0*factG*d*dx_loc**2)
                    if (lamjt > 1d0) then ! Jeans length resolved: gas is stable 
                       ok(i) = .false. 
                    else ! Jeans length not resolved --> form stars to lower density and stabilise gas
                       ! corresponding virial parameter for homogeneous sphere <= 1.5 in turbulent dominated 
                       ! limit and <= 0.5 in pressure dominated limit (in good agreement with observations,
                       ! at least for massive (>= 10^5 M_sun) clouds see Kauffmann, Pillai, Goldsmith 2013, Fig.1)  
                       alpha0 = 5d0/(onepi*factG*d)*(trgv + c_s2)/dx_loc**2
                       ! compute star formation efficiency per free-fall time (Federrath & Klessen 2012 eq 41) 
                       ! e_cts is the unresolved (for us) proto-stellar feedback: i.e. the dense gas core-to-star efficiency
                       e_cts = 0.5  ! would be 1.0 without feedback (Federrath & Klessen 2012)
                       phi_t = 0.57 ; theta = 0.33 ! best fit values of the Padoan & Nordlund multi-scale sf model to GMC simulation data 
                       sigs  = log(1.0+0.16*trgv/c_s2) ! factor 0.16 is b^2 where b=0.4 for a mixture of turbulence forcing modes
                       scrit = log(0.067/theta**2*alpha0*trgv/c_s2) ! best fit from Padoan & Nordlund MS model again 
                       sfr_ff(i) = e_cts/2.0*phi_t*exp(3.0/8.0*sigs)*(2.0-erfc((sigs-scrit)/sqrt(2.0*sigs)))
                    endif
                 endif
              endif
           end do

        else if(TRIM(star_maker)=='dencon')then
           ! Cen & Ostriker scheme (converging flow)
           do i=1,ngrid
              ! if cell is a leaf cell
              if (ok(i)) then 
                 d  = uold(ind_cell(i),1)
                 if (d <= d0) then  ! density threshold
                    ok(i) = .false.
                 endif

                 if(ok(i))then ! converging flow check
                    ncell = 1 ! we just want the neighbors of that cell
                    call getnbor(ind_cell(i),ind_nbor,ncell,ilevel)
                    ! First calculate velocities for divergence
                    d1    = uold(ind_nbor(1,1),1)     ; d4    = uold(ind_nbor(1,4),1) 
                    d2    = uold(ind_nbor(1,2),1)     ; d5    = uold(ind_nbor(1,5),1) 
                    d3    = uold(ind_nbor(1,3),1)     ; d6    = uold(ind_nbor(1,6),1)  
                    ul    = d1*uold(ind_nbor(1,1),2)
                    ur    = d2*uold(ind_nbor(1,2),2)
                    vl    = d3*uold(ind_nbor(1,3),3)
                    vr    = d4*uold(ind_nbor(1,4),3)
                    wl    = d5*uold(ind_nbor(1,5),4)
                    wr    = d6*uold(ind_nbor(1,6),4)
                    divv  = (ur-ul)+(vr-vl)+(wr-wl)
                    if(divv>0) ok(i) = .false.  ! diverging flow
                 endif

              endif
           enddo !i=1,ngrid

        else if(TRIM(star_maker)=='padoan')then
           ! variable star formation efficiency based on the local kinematics (Padoan+12)
           do i=1,ngrid
              ! if cell is a leaf cell
              if (ok(i)) then 
                 d = uold(ind_cell(i),1)
                 ! density threshold to avoid computing the criterion too often but otherwise not needed
                 if (d <= d_gmc) then
                    ok(i) = .false.
                 else
                    u = uold(ind_cell(i),2)
                    v = uold(ind_cell(i),3)
                    w = uold(ind_cell(i),4)
                    ncell = 1 ! we just want the neighbors of that cell
                    call getnbor(ind_cell(i),ind_nbor,ncell,ilevel)

                    ! density at 6 neighboring cells (if no cells, straight injection)
                    d1    = uold(ind_nbor(1,1),1)     ; d4    = uold(ind_nbor(1,4),1) 
                    d2    = uold(ind_nbor(1,2),1)     ; d5    = uold(ind_nbor(1,5),1) 
                    d3    = uold(ind_nbor(1,3),1)     ; d6    = uold(ind_nbor(1,6),1) 
       
                    ! local density maxima
                    if(ok(i)) then
                      if(d1>d.or.d2>d.or.d3>d.or.d4>d.or.d5>d.or.d6>d) ok(i)=.false.
                    endif

                    ! converging flow check
                    ul    = uold(ind_nbor(1,1),2)
                    ur    = uold(ind_nbor(1,2),2)
                    vl    = uold(ind_nbor(1,3),3)
                    vr    = uold(ind_nbor(1,4),3)
                    wl    = uold(ind_nbor(1,5),4)
                    wr    = uold(ind_nbor(1,6),4)
                    divv  = (ur*d2-ul*d1)+(vr*d4-vl*d3)+(wr*d6-wl*d5)
                    if(divv>0) ok(i) = .false.  ! diverging flow
 
                    if(ok(i))then
                       vturb(i)=0d0

                       !average velocity
                       dtot  = d+d1+d2+d3+d4+d5+d6
                       uavg  = d1*uold(ind_nbor(1,1),2) + d2*uold(ind_nbor(1,2),2) &
                            &+ d3*uold(ind_nbor(1,3),2) + d4*uold(ind_nbor(1,4),2) &
                            &+ d5*uold(ind_nbor(1,5),2) + d6*uold(ind_nbor(1,6),2) &
                            &+ d *uold(ind_cell(i  ),2)
                       vavg  = d1*uold(ind_nbor(1,1),3) + d2*uold(ind_nbor(1,2),3) &
                            &+ d3*uold(ind_nbor(1,3),3) + d4*uold(ind_nbor(1,4),3) &
                            &+ d5*uold(ind_nbor(1,5),3) + d6*uold(ind_nbor(1,6),3) &
                            &+ d *uold(ind_cell(i  ),3)
                       wavg  = d1*uold(ind_nbor(1,1),4) + d2*uold(ind_nbor(1,2),4) &
                            &+ d3*uold(ind_nbor(1,3),4) + d4*uold(ind_nbor(1,4),4) &
                            &+ d5*uold(ind_nbor(1,5),4) + d6*uold(ind_nbor(1,6),4) &
                            &+ d *uold(ind_cell(i  ),4)
                       uavg  = uavg/dtot
                       vavg  = vavg/dtot
                       wavg  = wavg/dtot

                       !left face
                       ul    = uold(ind_nbor(1,1),2) - uavg
                       vl    = uold(ind_nbor(1,1),3) - vavg
                       wl    = uold(ind_nbor(1,1),4) - wavg
                       if(ul>0)ul=0d0  ! take out infalling velocity
                       vturb(i) = vturb(i) + (ul*ul+vl*vl+wl*wl)
                       !right face
                       ur    = uold(ind_nbor(1,2),2) - uavg
                       vr    = uold(ind_nbor(1,2),3) - vavg
                       wr    = uold(ind_nbor(1,2),4) - wavg
                       if(ur<0)ur=0d0
                       vturb(i) = vturb(i) + (ur*ur+vr*vr+wr*wr)
                       !left face
                       ul    = uold(ind_nbor(1,3),2) - uavg
                       vl    = uold(ind_nbor(1,3),3) - vavg
                       wl    = uold(ind_nbor(1,3),4) - wavg
                       if(vl>0)vl=0d0 
                       vturb(i) = vturb(i) + (ul*ul+vl*vl+wl*wl)
                       !right face
                       ur    = uold(ind_nbor(1,4),2) - uavg
                       vr    = uold(ind_nbor(1,4),3) - vavg
                       wr    = uold(ind_nbor(1,4),4) - wavg
                       if(vr<0)vr=0d0
                       vturb(i) = vturb(i) + (ur*ur+vr*vr+wr*wr)
                       !left face
                       ul    = uold(ind_nbor(1,5),2) - uavg
                       vl    = uold(ind_nbor(1,5),3) - vavg
                       wl    = uold(ind_nbor(1,5),4) - wavg
                       if(wl>0)wl=0d0 
                       vturb(i) = vturb(i) + (ul*ul+vl*vl+wl*wl)
                       !right face
                       ur    = uold(ind_nbor(1,6),2) - uavg
                       vr    = uold(ind_nbor(1,6),3) - vavg
                       wr    = uold(ind_nbor(1,6),4) - wavg
                       if(wr<0)wr=0d0 
                       vturb(i) = vturb(i) + (ur*ur+vr*vr+wr*wr)
                       
                       vturb(i) = sqrt(vturb(i)/6.) 

                       vturb(i) = max(vturb(i),1d-20)
                       t_ff  = dsqrt(3d0*onepi/(32d0*6.67d-8*d*scale_d)) ! [s]
                       t_dyn = 3d0*dx_loc*scale_l/(2.*vturb(i)*scale_v)
                       effArr(i)=0.5*exp(-1.6*t_ff/t_dyn)
                    endif  ! if converging

                 endif ! if d>d_gmc
              endif ! if leaf cell
           end do

        endif ! star maker


        ! Calculate number of new stars in each cell using Poisson statistics
        do i=1,ngrid
           nstar(i) = 0
           if (ok(i)) then
              ! Compute mean number of events
              d         = uold(ind_cell(i),1)
              mcell     = d*vol_loc
              if(TRIM(star_maker)=='hopkins')then
                 ! this is the free fall time of an homogeneous sphere 
                 sfr_ff(i) = 1.0 
                 tstar     = 0.5427*sqrt(1.0/(factG*d))
              else if(TRIM(star_maker)=='padoan')then
                 ! this is for Padoan prescription
                 sfr_ff(i) = 1.0 
                 tstar     = t0*sqrt(d0/d)/(effArr(i)/eps_star)
              else if(TRIM(star_maker)=='federrath')then
                 ! this is the free fall time of an homogeneous sphere 
                 tstar     = 0.5427*sqrt(1.0/(factG*d))
              else
                 sfr_ff(i) = 1.0 
                 tstar     = t0*sqrt(d0/d) 
              endif
              PoissMean = dtnew(ilevel)*sfr_ff(i)/tstar*mcell/mstar
              ! If catastrophic star formation (massive star cluster) wants to occur, we need to limit the 
              ! maximal mass of the star particle we want to create in a cell.
              PoissMean = min(PoissMean,10.0)
              ! Compute Poisson realisation
              call poissdev(localseed,PoissMean,nstar(i))
              ! Compute depleted gas mass
              mgas      = nstar(i)*mstar

              ! Security to prevent more than 90% of gas depletion
              if (mgas > 0.9*mcell) then
                 nstar_corrected = int(0.9*mcell/mstar)
                 mstar_lost      = mstar_lost+(nstar(i)-nstar_corrected)*mstar
                 nstar(i)        = nstar_corrected
              endif

              ! Compute new stars local statistics
              if (nstar(i)>0) then
                 ntot = ntot+1
              else
                 nstar(i)=0
              endif
              mstar_tot = mstar_tot+nstar(i)*mstar
           endif
        enddo
        ! Store nstar in array flag2
        do i=1,ngrid
           flag2(ind_cell(i)) = nstar(i)
        end do
     end do
  end do

  !---------------------------------
  ! Check for free particle memory
  !---------------------------------
  ok_free = ((numbp_free-ntot)>=0)
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbp_free_tot = numbp_free
#endif
  if(.not. ok_free)then
     write(*,*)'No more free memory for particles'
     write(*,*)'Increase npartmax'
#ifndef WITHOUTMPI
    call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
    stop
#endif
  end if

  !---------------------------------
  ! Compute global stars statistics
  !---------------------------------
#ifndef WITHOUTMPI
  mlost = mstar_lost; mtot = mstar_tot
  call MPI_ALLREDUCE(ntot,ntot_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mtot,mtot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mlost,mlost_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
  ntot_all  = ntot
  mtot_all  = mstar_tot
  mlost_all = mstar_lost
#endif
  ntot_star_cpu = 0; ntot_star_all = 0
  ntot_star_cpu(myid) = ntot
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntot_star_cpu,ntot_star_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntot_star_cpu(1) = ntot_star_all(1)
#endif
  do icpu=2,ncpu
     ntot_star_cpu(icpu) = ntot_star_cpu(icpu-1)+ntot_star_all(icpu)
  end do
  nstar_tot = nstar_tot+ntot_all
  if (myid==1) then
     if (ntot_all.gt.0) then
        !print *,'##################################################################################'
        write(*,'(" Level=",I6," New star=",I6," Tot=",I10," Mass=",1PE10.3," Lost=",0PF4.1,"%")')&
             & ilevel,ntot_all,nstar_tot,mtot_all,mlost_all/mtot_all*100.
        !print *,'##################################################################################'
     endif
  end if

  !------------------------------
  ! Create new star particles
  !------------------------------
  ! Starting identity number
  if (myid==1) then
     index_star = nstar_tot-ntot_all
  else
     index_star = nstar_tot-ntot_all+ntot_star_cpu(myid-1)
  end if

  ! Loop over grids
  ncache = active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid = MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i) = active(ilevel)%igrid(igrid+i-1)
     end do

     ! Loop over cells
     do ind=1,twotondim
        iskip = ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i) = iskip+ind_grid(i)
        end do

        ! Flag cells with at least one new star
        do i=1,ngrid
           ok(i) = flag2(ind_cell(i))>0
        end do

        ! Gather new star arrays
        nnew = 0
        do i=1,ngrid
           if (ok(i)) then
              nnew               = nnew+1
              ind_grid_new(nnew) = ind_grid(i)
              ind_cell_new(nnew) = ind_cell(i)
           end if
        end do

        ! Update linked list for stars
        call remove_free(ind_part,nnew)
        call add_list(ind_part,ind_grid_new,ok_new,nnew)

        ! Calculate new star particle and modify gas density
        do i=1,nnew
           index_star = index_star+1

           ! Get gas variables
           n = flag2(ind_cell_new(i))
           d = uold(ind_cell_new(i),1)
           u = uold(ind_cell_new(i),2)
           v = uold(ind_cell_new(i),3)
           w = uold(ind_cell_new(i),4)
           x = (xg(ind_grid_new(i),1)+xc(ind,1)-skip_loc(1))*scale
           y = (xg(ind_grid_new(i),2)+xc(ind,2)-skip_loc(2))*scale
           z = (xg(ind_grid_new(i),3)+xc(ind,3)-skip_loc(3))*scale
           if (metal) then
              zg = uold(ind_cell_new(i),imetal)
#ifdef NCHEM
              do iche=1,nchem
                 chem1(iche) = uold(ind_cell_new(i),ichem+iche-1)
              enddo
#endif
           endif

           ! Set star particle variables
           tp(ind_part(i))     = birth_epoch  ! Birth epoch
           mp(ind_part(i))  = n*mstar      ! Mass
           if(use_initial_mass) then
              mp0(ind_part(i)) = mp(ind_part(i)) ! Initial Mass
           endif
           levelp(ind_part(i)) = ilevel       ! Level
           idp(ind_part(i))    = -index_star  ! Star identity
           st_n_tp(ind_part(i))=d       ! Cell density                   !SD
           st_n_SN(ind_part(i))=0d0                                      !SD
           st_e_SN(ind_part(i))=0d0                                      !SD
           xp(ind_part(i),1)   = x
           xp(ind_part(i),2)   = y
           xp(ind_part(i),3)   = z
           vp(ind_part(i),1)   = u
           vp(ind_part(i),2)   = v
           vp(ind_part(i),3)   = w
           if (metal) then
              zp(ind_part(i)) = zg  ! Initial star metallicity
#ifdef NCHEM
              do iche=1,nchem
                 chp(ind_part(i),iche) = chem1(iche)  ! Initial chemical abudance
              enddo
#endif
           endif
        end do
        ! End loop over new star particles

        ! Modify gas density according to mass depletion
        do i=1,ngrid
           if (flag2(ind_cell(i))>0) then
              n = flag2(ind_cell(i))
              d = uold(ind_cell(i),1)
              uold(ind_cell(i),1) = d-n*dstar
           endif
        end do

     end do
     ! End loop over cells
  end do
  ! End loop over grids
  
  !---------------------------------------------------------
  ! Convert hydro variables back to conservative variables
  !---------------------------------------------------------
  ncache = active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid = MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i) = active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim  
        iskip = ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i) = iskip+ind_grid(i)
        end do
        do i=1,ngrid
           d   = uold(ind_cell(i),1)
           u   = uold(ind_cell(i),2)
           v   = uold(ind_cell(i),3)
           w   = uold(ind_cell(i),4)
           e   = uold(ind_cell(i),5)*d
#ifdef SOLVERmhd
           bx1 = uold(ind_cell(i),6)
           by1 = uold(ind_cell(i),7)
           bz1 = uold(ind_cell(i),8)
           bx2 = uold(ind_cell(i),nvar+1)
           by2 = uold(ind_cell(i),nvar+2)
           bz2 = uold(ind_cell(i),nvar+3)
           e   = e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
           e   = e+0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=1,nener
              e=e+uold(ind_cell(i),5+irad)
           end do
#endif
           uold(ind_cell(i),1) = d
           uold(ind_cell(i),2) = d*u
           uold(ind_cell(i),3) = d*v
           uold(ind_cell(i),4) = d*w
           uold(ind_cell(i),5) = e
        end do
        do ivar=imetal,nvar
           do i=1,ngrid
              d = uold(ind_cell(i),1)
              w = uold(ind_cell(i),ivar)
              uold(ind_cell(i),ivar) = d*w
           end do
        end do
     end do
  end do

#endif


end subroutine star_formation 
!################################################################
!################################################################
!################################################################
!################################################################
subroutine getnbor(ind_cell,ind_father,ncell,ilevel)
  use amr_commons
  implicit none
  integer::ncell,ilevel
  integer,dimension(1:nvector)::ind_cell
  integer,dimension(1:nvector,0:twondim)::ind_father
  !-----------------------------------------------------------------
  ! This subroutine determines the 2*ndim neighboring cells
  ! cells of the input cell (ind_cell). 
  ! If for some reasons they don't exist, the routine returns 
  ! the input cell.
  !-----------------------------------------------------------------
  integer::nxny,i,idim,j,iok,ind
  integer,dimension(1:3)::ibound,iskip1,iskip2
  integer,dimension(1:nvector,1:3),save::ix
  integer,dimension(1:nvector),save::ind_grid_father,pos
  integer,dimension(1:nvector,0:twondim),save::igridn,igridn_ok
  integer,dimension(1:nvector,1:twondim),save::icelln_ok


  if(ilevel==1)then 
     write(*,*) 'Warning: attempting to form stars on level 1 --> this is not allowed ...'
     return
  endif

  ! Get father cell
  do i=1,ncell
     ind_father(i,0)=ind_cell(i)
  end do
  
  ! Get father cell position in the grid
  do i=1,ncell
     pos(i)=(ind_father(i,0)-ncoarse-1)/ngridmax+1
  end do
  
  ! Get father grid
  do i=1,ncell
     ind_grid_father(i)=ind_father(i,0)-ncoarse-(pos(i)-1)*ngridmax
  end do
  
  ! Get neighboring father grids
  call getnborgrids(ind_grid_father,igridn,ncell)
  
  ! Loop over position
  do ind=1,twotondim
     
     ! Select father cells that sit at position ind
     do j=0,twondim
        iok=0
        do i=1,ncell
           if(pos(i)==ind)then
              iok=iok+1
              igridn_ok(iok,j)=igridn(i,j)
           end if
        end do
     end do
     
     ! Get neighboring cells for selected cells
     if(iok>0)call getnborcells(igridn_ok,ind,icelln_ok,iok)
     
     ! Update neighboring father cells for selected cells
     do j=1,twondim
        iok=0
        do i=1,ncell
           if(pos(i)==ind)then
              iok=iok+1
              if(icelln_ok(iok,j)>0)then
                 ind_father(i,j)=icelln_ok(iok,j)
                 !write(*,*) 'index first if',ind_father(i,j) 
              else
                 ind_father(i,j)=ind_cell(i)
                 !write(*,*) 'index second if',ind_father(i,j) 
              end if
           end if
        end do
     end do
     
  end do
     
    
end subroutine getnbor
!##############################################################
!##############################################################
!##############################################################
!##############################################################
function erfc(x)

! complementary error function
  use amr_commons, ONLY: dp 
  implicit none
  real(dp) erfc
  real(dp) x, y
  real(kind=8) pv, ph 
  real(kind=8) q0, q1, q2, q3, q4, q5, q6, q7
  real(kind=8) p0, p1, p2, p3, p4, p5, p6, p7
  parameter(pv= 1.26974899965115684d+01, ph= 6.10399733098688199d+00) 
  parameter(p0= 2.96316885199227378d-01, p1= 1.81581125134637070d-01) 
  parameter(p2= 6.81866451424939493d-02, p3= 1.56907543161966709d-02) 
  parameter(p4= 2.21290116681517573d-03, p5= 1.91395813098742864d-04) 
  parameter(p6= 9.71013284010551623d-06, p7= 1.66642447174307753d-07)
  parameter(q0= 6.12158644495538758d-02, q1= 5.50942780056002085d-01) 
  parameter(q2= 1.53039662058770397d+00, q3= 2.99957952311300634d+00) 
  parameter(q4= 4.95867777128246701d+00, q5= 7.41471251099335407d+00) 
  parameter(q6= 1.04765104356545238d+01, q7= 1.48455557345597957d+01)
  
  y = x*x
  y = exp(-y)*x*(p7/(y+q7)+p6/(y+q6) + p5/(y+q5)+p4/(y+q4)+p3/(y+q3) &   
       &       + p2/(y+q2)+p1/(y+q1)+p0/(y+q0))
  if (x < ph) y = y+2d0/(exp(pv*x)+1.0)
  erfc = y
  
  return
  
end function erfc
!##############################################################
!##############################################################
!##############################################################
subroutine MaxwellianCDF(x,ywant,fn,df)
   implicit none
   real(kind=8)::x, fn, df, ywant
   real(kind=8)::pi=3.14159265358979323
   fn = derf(x) - 2d0/dsqrt(pi)*x*dexp(-x*x) - ywant
   df = dsqrt(8d0/pi)*x*x*dexp(-x*x)
end subroutine MaxwellianCDF
!##############################################################
!##############################################################
!##############################################################
subroutine newton_raphson_safe(x1,x2,ywant,xeps,rtsafe)
   IMPLICIT NONE
   REAL(KIND=8)::x1,x2,ywant,xeps,rtsafe
   !Using a combination of Newton-Raphson and bisection, 
   !and the root of a function bracketed between x1 and x2. 
   !The root, returned as the function value rtsafe, will be returned 
   !until its accuracy is known within xeps. funcd is a user-supplied 
   !subroutine which returns both the function value and the 
   !first derivative of the function.
   INTEGER::j,MAXIT=100
   REAL(KIND=8)::df,dx,dxold,f,fh,fl,temp,xh,xl

   call MaxwellianCDF(x1,ywant,fl,df)
   call MaxwellianCDF(x2,ywant,fh,df)
   if((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.)) &
      write(*,*) 'root must be bracketed in rtsafe'
   if(fl.eq.0.)then
      rtsafe=x1
      return
   else if(fh.eq.0.)then
      rtsafe=x2
      return
   else if(fl.lt.0.)then !Orient the search so that f(xl) < 0.
      xl=x1
      xh=x2
   else
      xh=x1
      xl=x2
   endif
   rtsafe=.5*(x1+x2) !Initialize the guess for root,
   dxold=dabs(x2-x1) !the \stepsize before last,"
   dx=dxold !and the last step.
   call MaxwellianCDF(rtsafe,ywant,f,df)
   do j=1,MAXIT !Loop over allowed iterations.
      if((((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).gt.0)& !Bisect if Newton out of range,
         .or.(dabs(2d0*f).gt.dabs(dxold*df)) ) then !or not decreasing fast enough.
         dxold=dx
         dx=0.5d0*(xh-xl)
         rtsafe=xl+dx
         if(xl.eq.rtsafe)return !Change in root is negligible.
      else !Newton step acceptable. Take it.
         dxold=dx
         dx=f/df
         temp=rtsafe
         rtsafe=rtsafe-dx
         if(temp.eq.rtsafe)return
      endif
      if(dabs(dx).lt.xeps) return !Convergence criterion.
      call MaxwellianCDF(rtsafe,ywant,f,df) !The one new function evaluation per iteration.
      if(f.lt.0.) then !Maintain the bracket on the root.
         xl=rtsafe
      else
         xh=rtsafe
      endif
    enddo
    write(*,*) 'rtsafe exceeding maximum iterations'
end subroutine newton_raphson_safe
!##############################################################
!##############################################################
!##############################################################
subroutine open_stat_files
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  character(len=5)::nchar
  character(len=256)::filename,mkdircmd
  character(len=7)::SFdir='Useful/'
  integer::ilun

  call title(myid,nchar)
  filename=SFdir//'/stat_'//nchar//'.dat'
  ilun = 3*ncpu+myid+10
  open(ilun,file=TRIM(filename),access='APPEND')

end subroutine open_stat_files
