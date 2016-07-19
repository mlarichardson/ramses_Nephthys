        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul 19 00:14:11 2016
        MODULE CIC_AMR__genmod
          INTERFACE 
            SUBROUTINE CIC_AMR(IND_CELL,IND_PART,IND_GRID_PART,X0,NG,NP,&
     &ILEVEL)
              INTEGER(KIND=4) :: IND_CELL(1:32)
              INTEGER(KIND=4) :: IND_PART(1:32)
              INTEGER(KIND=4) :: IND_GRID_PART(1:32)
              REAL(KIND=8) :: X0(1:32,1:3)
              INTEGER(KIND=4) :: NG
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: ILEVEL
            END SUBROUTINE CIC_AMR
          END INTERFACE 
        END MODULE CIC_AMR__genmod
