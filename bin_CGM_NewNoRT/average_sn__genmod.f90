        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul 19 00:14:22 2016
        MODULE AVERAGE_SN__genmod
          INTERFACE 
            SUBROUTINE AVERAGE_SN(XSN,VSN,VOL_GAS,DQ,EKBLAST,IND_BLAST, &
     &NSN,NSN_TOT,ISN_MYID,MSN,MLOADSN,ZSN,ZLOADSN,VLOADSN,IND_PART)
              INTEGER(KIND=4) :: NSN_TOT
              INTEGER(KIND=4) :: NSN
              REAL(KIND=8) :: XSN(1:NSN,1:3)
              REAL(KIND=8) :: VSN(1:NSN,1:3)
              REAL(KIND=8) :: VOL_GAS(1:NSN)
              REAL(KIND=8) :: DQ(1:NSN,1:3)
              REAL(KIND=8) :: EKBLAST(1:NSN)
              INTEGER(KIND=4) :: IND_BLAST(1:NSN)
              INTEGER(KIND=4) :: ISN_MYID(1:NSN)
              REAL(KIND=8) :: MSN(1:NSN)
              REAL(KIND=8) :: MLOADSN(1:NSN)
              REAL(KIND=8) :: ZSN(1:NSN)
              REAL(KIND=8) :: ZLOADSN(1:NSN)
              REAL(KIND=8) :: VLOADSN(1:NSN,1:3)
              INTEGER(KIND=4) :: IND_PART(1:NSN)
            END SUBROUTINE AVERAGE_SN
          END INTERFACE 
        END MODULE AVERAGE_SN__genmod
