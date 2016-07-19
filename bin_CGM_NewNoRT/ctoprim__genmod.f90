        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul 19 00:14:04 2016
        MODULE CTOPRIM__genmod
          INTERFACE 
            SUBROUTINE CTOPRIM(UIN,Q,C,GRAVIN,DT,NGRID)
              REAL(KIND=8) :: UIN(1:32,-1:4,-1:4,-1:4,1:15)
              REAL(KIND=8) :: Q(1:32,-1:4,-1:4,-1:4,1:15)
              REAL(KIND=8) :: C(1:32,-1:4,-1:4,-1:4)
              REAL(KIND=8) :: GRAVIN(1:32,-1:4,-1:4,-1:4,1:3)
              REAL(KIND=8) :: DT
              INTEGER(KIND=4) :: NGRID
            END SUBROUTINE CTOPRIM
          END INTERFACE 
        END MODULE CTOPRIM__genmod
