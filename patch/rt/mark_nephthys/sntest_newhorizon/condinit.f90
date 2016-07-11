!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_commons
  use pm_commons
  use poisson_parameters, ONLY: gravity_params
  use random
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar,i
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp)::xx,yy,zz,r,rr,v,xc,yc,zc,M,rho,dzz,zint,HH,rdisk,dpdr,dmax
  real(dp)::rgal,sum,sum2,dmin,zmin,zmax,c,fgas,pi,tol
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,M_b,az,eps
  real(dp)::rmin,rmax,Tiso,sigmaNT,scale_s,dxmin2
  real(dp)::RandNum,phi,Rrand,SS,CC,UU
  integer::ix,iy,iz
  real(dp)::d_background,d_surround,d_blob
  ! Add here, if you wish, some user-defined initial conditions
  dxmin2 = boxlen/2d0**levelmin/2d0
  xc  =boxlen/2.0 + dxmin2
  yc  =boxlen/2.0 + dxmin2
  zc  =boxlen/2.0 + dxmin2

  Tiso = gravity_params(2)/1.2195

  d_background = gravity_params(1)

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  pi=acos(-1.)
  do i=1,nn
     xx=x(i,1)-xc
     yy=x(i,2)-yc
     zz=x(i,3)-zc

     q(i,1)=d_background/scale_nH       
     q(i,2:4)=0d0 ! static
     q(i,ndim+2)=Tiso/scale_T2*q(i,1)

     if(metal)q(i,imetal)=z_ave*0.02
     do ivar=imetal+1,nvar
        q(i,ivar)=0d0
     end do

     ! Default value for hydrogen fraction (chemistry)
     do ivar=ichem,ichem+nchem-1
        if(chem_list(ivar-ichem+1)=='H ') q(i,ivar)=0.76-z_ave*0.02
     enddo


  enddo

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif
  ! kinetic energy
  u(1:nn,ndim+2)=0.0d0
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,4)**2
#endif
  ! pressure -> total fluid energy
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)
  ! passive scalars
  do ivar=imetal,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do

end subroutine condinit
