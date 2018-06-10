        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  9 18:17:14 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LFIT__genmod
          INTERFACE 
            SUBROUTINE LFIT(X,Y,SIG,NDAT,A,IA,MA,COVAR,NPC,CHISQ)
              INTEGER(KIND=4) :: NPC
              INTEGER(KIND=4) :: MA
              INTEGER(KIND=4) :: NDAT
              REAL(KIND=8) :: X(NDAT)
              REAL(KIND=8) :: Y(NDAT)
              REAL(KIND=8) :: SIG(NDAT)
              REAL(KIND=8) :: A(MA)
              INTEGER(KIND=4) :: IA(MA)
              REAL(KIND=8) :: COVAR(NPC,NPC)
              REAL(KIND=8) :: CHISQ
            END SUBROUTINE LFIT
          END INTERFACE 
        END MODULE LFIT__genmod
