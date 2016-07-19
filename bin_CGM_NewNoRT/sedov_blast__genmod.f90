        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul 19 00:14:22 2016
        MODULE SEDOV_BLAST__genmod
          INTERFACE 
            SUBROUTINE SEDOV_BLAST(XSN,MSN,INDSN,VOL_GAS,DQ,EKBLAST,NSN,&
     &MLOADSN,ZLOADSN,VLOADSN)
              INTEGER(KIND=4) :: NSN
              REAL(KIND=8) :: XSN(1:NSN,1:3)
              REAL(KIND=8) :: MSN(1:NSN)
              INTEGER(KIND=4) :: INDSN(1:NSN)
              REAL(KIND=8) :: VOL_GAS(1:NSN)
              REAL(KIND=8) :: DQ(1:NSN,1:3)
              REAL(KIND=8) :: EKBLAST(1:NSN)
              REAL(KIND=8) :: MLOADSN(1:NSN)
              REAL(KIND=8) :: ZLOADSN(1:NSN)
              REAL(KIND=8) :: VLOADSN(1:NSN,1:3)
            END SUBROUTINE SEDOV_BLAST
          END INTERFACE 
        END MODULE SEDOV_BLAST__genmod
