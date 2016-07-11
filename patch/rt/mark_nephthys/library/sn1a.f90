module amr_commons
   integer, parameter::dp=8
   real(dp)::A_snIa=0.0013
end module amr_commons

program main
   implicit none
   integer::id_star
   real(8)::mass0,nsnIa
  
   mass0=1d4

   id_star=-1
   call get_number_of_snIa(id_star,mass0,nsnIa)


end program main

subroutine get_number_of_snIa ( id_star, mass0, nsnIa )
   use amr_commons, ONLY:dp, A_snIa
   implicit none
   integer ::id_star
   real(dp)::mass0 ! initial mass of a star particle in Msun
   real(dp)::nsnIa ! number of snIa
!-------------------------------------------------------------
!	Use the inverse method to generate random numbers for snIa
!-------------------------------------------------------------
   real(dp)::A_DTD,t_ini,t_fin,xdum,ydum
   integer ::localseed,i,nsnIa_tot
   real(dp),external::ran1


!  DTD = A_DTD* t^-1
!  1 = int A_DTD / t dt
!    = A_DTD*(log_e (t_f) - log_e (t_i))
!  A_DTD = 1d0/(alog(t_fin) - alog(t_ini))

!  n_sn(<t) = A_DTD*(log (t) - log(t_ini)) ; 0 < n_sn < 1
!  log(t) = n_sn/A_DTD + log(t_ini)
!  t = exp(n_sn/A_DTD + log(t_ini))


   t_ini = 50d6
   t_fin = 1.37d10
   A_DTD = 1d0 / (log(t_fin) - log(t_ini))

   nsnIa_tot = NINT(mass0 * A_snIa)
   localseed = -ABS(id_star)
 
   do i=1,nsnIa_tot
      xdum = ran1(localseed)
      ydum = exp(xdum / A_DTD + log(t_ini))
      print *, ydum
   end do 

end subroutine get_number_of_snIa


function ran1(idum)
   implicit none 
   integer:: idum,IA,IM,IQ,IR,NTAB,NDIV
   real(kind=8):: ran1,AM,EPS,RNMX
   parameter(IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,&
            &NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
   integer::j,k,iv(NTAB),iy
   save iv,iy
   data iv /NTAB*0/, iy /0/  
   if (idum.le.0.or.iy.eq.0) then! initialize
      idum=max(-idum,1)
      do j=NTAB+8,1,-1
         k=idum/IQ
         idum=IA*(idum-k*IQ)-IR*k
         if (idum.lt.0) idum=idum+IM
         if (j.le.NTAB) iv(j)=idum
      end do
      iy=iv(1)
   end if
   k=idum/IQ
   idum=IA*(idum-k*IQ)-IR*k
   if (idum.lt.0) idum=idum+IM
   j=1+iy/NDIV
   iy=iv(j)
   iv(j)=idum
   ran1=min(AM*iy,RNMX)
   return
end
