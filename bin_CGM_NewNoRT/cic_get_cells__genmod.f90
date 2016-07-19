        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul 19 00:14:20 2016
        MODULE CIC_GET_CELLS__genmod
          INTERFACE 
            SUBROUTINE CIC_GET_CELLS(INDP,XX,VOL,OK,IND_GRID,XPART,     &
     &IND_GRID_PART,NG,NP,ILEVEL)
              INTEGER(KIND=4) :: INDP(1:32,1:8)
              REAL(KIND=8) :: XX(1:32,1:3,8)
              REAL(KIND=8) :: VOL(1:32,1:8)
              LOGICAL(KIND=4) :: OK(1:32,1:8)
              INTEGER(KIND=4) :: IND_GRID(1:32)
              REAL(KIND=8) :: XPART(1:32,1:3)
              INTEGER(KIND=4) :: IND_GRID_PART(1:32)
              INTEGER(KIND=4) :: NG
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: ILEVEL
            END SUBROUTINE CIC_GET_CELLS
          END INTERFACE 
        END MODULE CIC_GET_CELLS__genmod
