        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul 19 00:14:03 2016
        MODULE CMPDIVU__genmod
          INTERFACE 
            SUBROUTINE CMPDIVU(Q,DIV,DX,DY,DZ,NGRID)
              REAL(KIND=8) :: Q(1:32,-1:4,-1:4,-1:4,1:15)
              REAL(KIND=8) :: DIV(1:32,1:3,1:3,1:3)
              REAL(KIND=8) :: DX
              REAL(KIND=8) :: DY
              REAL(KIND=8) :: DZ
              INTEGER(KIND=4) :: NGRID
            END SUBROUTINE CMPDIVU
          END INTERFACE 
        END MODULE CMPDIVU__genmod
