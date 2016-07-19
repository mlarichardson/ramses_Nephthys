        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul 19 00:13:49 2016
        MODULE FRIEDMAN__genmod
          INTERFACE 
            SUBROUTINE FRIEDMAN(O_MAT_0,O_VAC_0,O_K_0,ALPHA,AXP_MIN,    &
     &AXP_OUT,HEXP_OUT,TAU_OUT,T_OUT,NTABLE)
              INTEGER(KIND=4) :: NTABLE
              REAL(KIND=8) :: O_MAT_0
              REAL(KIND=8) :: O_VAC_0
              REAL(KIND=8) :: O_K_0
              REAL(KIND=8) :: ALPHA
              REAL(KIND=8) :: AXP_MIN
              REAL(KIND=8) :: AXP_OUT(0:NTABLE)
              REAL(KIND=8) :: HEXP_OUT(0:NTABLE)
              REAL(KIND=8) :: TAU_OUT(0:NTABLE)
              REAL(KIND=8) :: T_OUT(0:NTABLE)
            END SUBROUTINE FRIEDMAN
          END INTERFACE 
        END MODULE FRIEDMAN__genmod
