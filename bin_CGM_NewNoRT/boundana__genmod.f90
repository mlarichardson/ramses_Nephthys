        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul 19 00:14:08 2016
        MODULE BOUNDANA__genmod
          INTERFACE 
            SUBROUTINE BOUNDANA(X,U,DX,IBOUND,NCELL)
              REAL(KIND=8) :: X(1:32,1:3)
              REAL(KIND=8) :: U(1:32,1:15)
              REAL(KIND=8) :: DX
              INTEGER(KIND=4) :: IBOUND
              INTEGER(KIND=4) :: NCELL
            END SUBROUTINE BOUNDANA
          END INTERFACE 
        END MODULE BOUNDANA__genmod
