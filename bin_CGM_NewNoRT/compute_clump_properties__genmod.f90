        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul 19 00:14:25 2016
        MODULE COMPUTE_CLUMP_PROPERTIES__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_CLUMP_PROPERTIES(XX)
              USE AMR_COMMONS
              USE POISSON_COMMONS, ONLY :                               &
     &          PHI,                                                    &
     &          F
              REAL(KIND=8) :: XX(1:NCOARSE+NGRIDMAX*8)
            END SUBROUTINE COMPUTE_CLUMP_PROPERTIES
          END INTERFACE 
        END MODULE COMPUTE_CLUMP_PROPERTIES__genmod
