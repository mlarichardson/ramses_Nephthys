        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul 19 00:14:03 2016
        MODULE CONSUP__genmod
          INTERFACE 
            SUBROUTINE CONSUP(UIN,FLUX,DIV,DT,NGRID)
              REAL(KIND=8) :: UIN(1:32,-1:4,-1:4,-1:4,1:15)
              REAL(KIND=8) :: FLUX(1:32,1:3,1:3,1:3,1:15,1:3)
              REAL(KIND=8) :: DIV(1:32,1:3,1:3,1:3)
              REAL(KIND=8) :: DT
              INTEGER(KIND=4) :: NGRID
            END SUBROUTINE CONSUP
          END INTERFACE 
        END MODULE CONSUP__genmod
