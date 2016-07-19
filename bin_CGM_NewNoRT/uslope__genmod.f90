        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul 19 00:14:04 2016
        MODULE USLOPE__genmod
          INTERFACE 
            SUBROUTINE USLOPE(Q,DQ,DX,DT,NGRID)
              REAL(KIND=8) :: Q(1:32,-1:4,-1:4,-1:4,1:15)
              REAL(KIND=8) :: DQ(1:32,-1:4,-1:4,-1:4,1:15,1:3)
              REAL(KIND=8) :: DX
              REAL(KIND=8) :: DT
              INTEGER(KIND=4) :: NGRID
            END SUBROUTINE USLOPE
          END INTERFACE 
        END MODULE USLOPE__genmod
