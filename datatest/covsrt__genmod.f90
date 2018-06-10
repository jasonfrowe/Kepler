        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  9 18:24:41 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COVSRT__genmod
          INTERFACE 
            SUBROUTINE COVSRT(COVAR,NPC,MA,IA,MFIT)
              INTEGER(KIND=4) :: MA
              INTEGER(KIND=4) :: NPC
              REAL(KIND=8) :: COVAR(NPC,NPC)
              INTEGER(KIND=4) :: IA(MA)
              INTEGER(KIND=4) :: MFIT
            END SUBROUTINE COVSRT
          END INTERFACE 
        END MODULE COVSRT__genmod
