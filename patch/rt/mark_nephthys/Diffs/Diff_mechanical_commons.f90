5,6c5,6
<    use amr_commons, ONLY:dp,ndim,ncpu,nchem
<    integer(kind=4),parameter::nvarSN=13+nchem
---
>    use amr_commons, ONLY:dp,ndim,ncpu
>    integer(kind=4),parameter::nvarSN=11
15c15
<    real(dp)::f_LOAD,f_LEFT,f_ESN,f_PCAN
---
>    real(dp)::f_LOAD,f_LEFT,f_ESN 
19c19
<    integer ::ncomm_SN   ! the number of cells to be communicated (specific to myid)
---
>    integer ::nSN_comm   ! the number of cells to be communicated (specific to myid)
30,31d29
<    real(dp),parameter:: expN_SN_boost=-0.15
<    real(dp),parameter:: expE_SN_boost=16d0/17d0
34d31
<   ! common lists for SNe across different cpus
38c35
<    real(dp),dimension(1:ncomm_max)::nSN_comm          ! number of SNe
---
>    !integer ,dimension(1:ncomm_max)::lSN_comm  ! level
49,51c46
<    integer,dimension(:)  ,allocatable::ncomm_SN_cpu,ncomm_SN_mpi
<    real(dp),dimension(1:ncomm_max)::rSt_comm          ! Stromgren radius in pc
<    real(dp),dimension(1:ncomm_max,1:nchem)::mchloadSN_comm     ! chemical species ejected + entrained
---
>    integer,dimension(:)  ,allocatable::nSN_comm_cpu,nSN_comm_mpi
71,77d65
<    ! For chemical abundance due to SN II
<    real(dp)::Zejecta_chem_II(1:nchem)
< 
<    ! For chemical abundance due to SN Ia
<    real(dp)::mejecta_Ia
<    real(dp)::Zejecta_chem_Ia(1:nchem)
< 
85d72
<    use hydro_parameters, ONLY:ichem
91,92c78
<    integer::ncshell,iring,i2,j2,k2,nrad,irad,ich
<    character(len=2)::element_name
---
>    integer::ncshell,iring,i2,j2,k2,nrad,irad
99c85
<    if(.not.metal.and.mechanical_feedback)then
---
>    if(.not.metal.and.mechanical_feedback>0)then
106,107c92,93
<    allocate(ncomm_SN_cpu(1:ncpu))
<    allocate(ncomm_SN_mpi(1:ncpu))
---
>    allocate(nSN_comm_cpu(1:ncpu))
>    allocate(nSN_comm_mpi(1:ncpu))
114,116d99
<    f_PCAN = 0.9387  ! correction due to direct momentum cancellation 
<                     ! due to the assumption of 48 neighboring cells 
<                     ! even in the uniform case where there are 18 immediate neighbors
144c127
<          !indall(i+(j-1)*4+(k-1)*4*4) = ind      
---
>          !indall(i+(j-1)*4+(k-1)*4*4) = ind       
206,240d188
<    if(.not.variable_yield_SNII)then ! assume solar case
<       do ich=1,nchem
<          element_name=chem_list(ich)
<          select case (element_name)
<             case ('H ')
<                Zejecta_chem_II(ich) = 10.**(-0.30967822)
<             case ('He')
<                Zejecta_chem_II(ich) = 10.**(-0.40330181)
<             case ('C ')
<                Zejecta_chem_II(ich) = 10.**(-1.9626259)
<             case ('N ')
<                Zejecta_chem_II(ich) = 10.**(-2.4260355)
<             case ('O ')
<                Zejecta_chem_II(ich) = 10.**(-1.1213435)
<             case ('Mg')
<                Zejecta_chem_II(ich) = 10.**(-2.3706062)
<             case ('Si')
<                Zejecta_chem_II(ich) = 10.**(-2.0431845)
<             case ('S ')
<                Zejecta_chem_II(ich) = 10.**(-2.2964020)
<             case ('Fe')
<                Zejecta_chem_II(ich) = 10.**(-2.2126987)
<             case default
<                Zejecta_chem_II(ich)=0
<          end select
< 
<       end do
<    endif
< 
<   if (snIa) call init_snIa_yield
< 
<   ! This is not to double count SN feedback
<   if (mechanical_feedback) f_w = -1
< 
< 
245,301d192
< subroutine init_snIa_yield
<    use amr_commons
<    use mechanical_commons, ONLY: mejecta_Ia, Zejecta_chem_Ia
<    implicit none
<    integer::ich
<    real(kind=8)::yield_snIa(1:66)
<    character(len=2)::element_name
< !----------------------------------------------------------------------------
< !  Nomoto et al. (1997,1984) W7 (carbon-deflagration model)
< !            (better fit with Tycho observation)
< !----------------------------------------------------------------------------
< !       12C ,     13C ,   14N ,    15N ,    16O ,    17O ,    18O ,    19F ,
< !       20Ne,     21Ne,   22Ne,    23Na,    24Mg,    25Mg,    26Mg,    27Al,
< !       28Si,     29Si,   30Si,    31P ,    32S ,    33S ,    34S ,    36S ,
< !       35Cl,     37Cl,   36Ar,    38Ar,    40Ar,    39K ,    41K ,    40Ca,
< !       42Ca,     43Ca,   44Ca,    46Ca,    48Ca,    45Sc,    46Ti,    47Ti,
< !       48Ti,     49Ti,   50Ti,    50V ,    51V ,    50Cr,    52Cr,    53Cr,
< !       54Cr,     55Mn,   54Fe,    56Fe,    57Fe,    58Fe,    59Co,    58Ni,
< !       60Ni,     61Ni,   62Ni,    64Ni,    63Cu,    65Cu,    64Zn,    66Zn,
< !       67Zn,     68Zn                                                   
< !----------------------------------------------------------------------------
< 
<    yield_snIa = (/ &  ! Msun per SN
<     &4.83E-02,1.40E-06,1.16E-06,1.32E-09,1.43E-01,3.54E-08,8.25E-10,5.67E-10,&
<     &2.02E-03,8.46E-06,2.49E-03,6.32E-05,8.50E-03,4.05E-05,3.18E-05,9.86E-04,&
<     &1.50E-01,8.61E-04,1.74E-03,4.18E-04,8.41E-02,4.50E-04,1.90E-03,3.15E-07,&
<     &1.34E-04,3.98E-05,1.49E-02,1.06E-03,1.26E-08,8.52E-05,7.44E-06,1.23E-02,&
<     &3.52E-05,1.03E-07,8.86E-06,1.99E-09,7.10E-12,2.47E-07,1.71E-05,6.04E-07,&
<     &2.03E-04,1.69E-05,1.26E-05,8.28E-09,5.15E-05,2.71E-04,5.15E-03,7.85E-04,&
<     &1.90E-04,8.23E-03,1.04E-01,6.13E-01,2.55E-02,9.63E-04,1.02E-03,1.28E-01,&
<     &1.05E-02,2.51E-04,2.66E-03,1.31E-06,1.79E-06,6.83E-07,1.22E-05,2.12E-05,&
<     &1.34E-08,1.02E-08/)
< 
<    mejecta_Ia=sum(yield_snIa)
<    do ich=1,nchem
<       element_name=chem_list(ich)
<       select case (element_name)
<          case ('C ')
<             Zejecta_chem_Ia(ich)=sum(yield_snIa(1:2))/mejecta_Ia
<          case ('N ')
<             Zejecta_chem_Ia(ich)=sum(yield_snIa(3:4))/mejecta_Ia
<          case ('O ')
<             Zejecta_chem_Ia(ich)=sum(yield_snIa(5:7))/mejecta_Ia
<          case ('Mg')
<             Zejecta_chem_Ia(ich)=sum(yield_snIa(13:15))/mejecta_Ia
<          case ('Si')
<             Zejecta_chem_Ia(ich)=sum(yield_snIa(17:19))/mejecta_Ia
<          case ('S ')
<             Zejecta_chem_Ia(ich)=sum(yield_snIa(21:24))/mejecta_Ia
<          case ('Fe')
<             Zejecta_chem_Ia(ich)=sum(yield_snIa(51:54))/mejecta_Ia
<          case default
<             Zejecta_chem_Ia(ich)=0d0
<       end select     
<    enddo
< 
< end subroutine init_snIa_yield
