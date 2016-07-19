        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul 19 00:14:30 2016
        MODULE MECH_FINE_SNIA__genmod
          INTERFACE 
            SUBROUTINE MECH_FINE_SNIA(IND_GRID,IND_POS_CELL,NP,ILEVEL,  &
     &DTEFF,NSN,MSN,PSN,MZSN,NPHSN,MCHSN)
              INTEGER(KIND=4) :: IND_GRID(1:32)
              INTEGER(KIND=4) :: IND_POS_CELL(1:32)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: ILEVEL
              REAL(KIND=8) :: DTEFF
              REAL(KIND=8) :: NSN(1:32)
              REAL(KIND=8) :: MSN(1:32)
              REAL(KIND=8) :: PSN(1:32,1:3)
              REAL(KIND=8) :: MZSN(1:32)
              REAL(KIND=8) :: NPHSN(1:32)
              REAL(KIND=8) :: MCHSN(1:32,1:8)
            END SUBROUTINE MECH_FINE_SNIA
          END INTERFACE 
        END MODULE MECH_FINE_SNIA__genmod
