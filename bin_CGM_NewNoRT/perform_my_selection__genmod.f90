        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul 19 00:13:59 2016
        MODULE PERFORM_MY_SELECTION__genmod
          INTERFACE 
            SUBROUTINE PERFORM_MY_SELECTION(JUSTCOUNT,Z1,Z2,OM0IN,OMLIN,&
     &HUBIN,LBOX,OBSERVER,THETAY,THETAZ,THETA,PHI,POS,VEL,NPART,POSOUT, &
     &VELOUT,ZOUT,NPARTOUT,VERBOSE)
              INTEGER(KIND=4) :: NPARTOUT
              LOGICAL(KIND=4) :: JUSTCOUNT
              REAL(KIND=8) :: Z1
              REAL(KIND=8) :: Z2
              REAL(KIND=8) :: OM0IN
              REAL(KIND=8) :: OMLIN
              REAL(KIND=8) :: HUBIN
              REAL(KIND=8) :: LBOX
              REAL(KIND=8) :: OBSERVER(3)
              REAL(KIND=8) :: THETAY
              REAL(KIND=8) :: THETAZ
              REAL(KIND=8) :: THETA
              REAL(KIND=8) :: PHI
              REAL(KIND=8) :: POS(1:3,1:32)
              REAL(KIND=8) :: VEL(1:3,1:32)
              INTEGER(KIND=4) :: NPART
              REAL(KIND=8) :: POSOUT(3,NPARTOUT)
              REAL(KIND=8) :: VELOUT(3,NPARTOUT)
              REAL(KIND=8) :: ZOUT(NPARTOUT)
              LOGICAL(KIND=4) :: VERBOSE
            END SUBROUTINE PERFORM_MY_SELECTION
          END INTERFACE 
        END MODULE PERFORM_MY_SELECTION__genmod
