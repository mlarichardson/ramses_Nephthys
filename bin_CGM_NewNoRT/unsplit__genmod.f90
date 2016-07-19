        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul 19 00:14:04 2016
        MODULE UNSPLIT__genmod
          INTERFACE 
            SUBROUTINE UNSPLIT(UIN,GRAVIN,FLUX,TMP,DX,DY,DZ,DT,NGRID)
              REAL(KIND=8) :: UIN(1:32,-1:4,-1:4,-1:4,1:15)
              REAL(KIND=8) :: GRAVIN(1:32,-1:4,-1:4,-1:4,1:3)
              REAL(KIND=8) :: FLUX(1:32,1:3,1:3,1:3,1:15,1:3)
              REAL(KIND=8) :: TMP(1:32,1:3,1:3,1:3,1:2,1:3)
              REAL(KIND=8) :: DX
              REAL(KIND=8) :: DY
              REAL(KIND=8) :: DZ
              REAL(KIND=8) :: DT
              INTEGER(KIND=4) :: NGRID
            END SUBROUTINE UNSPLIT
          END INTERFACE 
        END MODULE UNSPLIT__genmod
