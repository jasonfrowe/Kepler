        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  9 18:17:15 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MRQCOF__genmod
          INTERFACE 
            SUBROUTINE MRQCOF(X,Y,SIG,DTYPE,NDATA,A,IA,MA,ALPHA,BETA,   &
     &NALP,CHISQ)
              INTEGER(KIND=4) :: NALP
              INTEGER(KIND=4) :: MA
              INTEGER(KIND=4) :: NDATA
              REAL(KIND=8) :: X(NDATA)
              REAL(KIND=8) :: Y(NDATA)
              REAL(KIND=8) :: SIG(NDATA)
              INTEGER(KIND=4) :: DTYPE(NDATA)
              REAL(KIND=8) :: A(MA)
              INTEGER(KIND=4) :: IA(MA)
              REAL(KIND=8) :: ALPHA(NALP,NALP)
              REAL(KIND=8) :: BETA(MA)
              REAL(KIND=8) :: CHISQ
            END SUBROUTINE MRQCOF
          END INTERFACE 
        END MODULE MRQCOF__genmod
