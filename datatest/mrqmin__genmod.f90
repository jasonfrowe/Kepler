        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  9 18:17:15 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MRQMIN__genmod
          INTERFACE 
            SUBROUTINE MRQMIN(X,Y,SIG,DTYPE,NDATA,A,IA,MA,COVAR,ALPHA,  &
     &NCA,CHISQ,ALAMDA)
              INTEGER(KIND=4) :: NCA
              INTEGER(KIND=4) :: MA
              INTEGER(KIND=4) :: NDATA
              REAL(KIND=8) :: X(NDATA)
              REAL(KIND=8) :: Y(NDATA)
              REAL(KIND=8) :: SIG(NDATA)
              INTEGER(KIND=4) :: DTYPE(NDATA)
              REAL(KIND=8) :: A(MA)
              INTEGER(KIND=4) :: IA(MA)
              REAL(KIND=8) :: COVAR(NCA,NCA)
              REAL(KIND=8) :: ALPHA(NCA,NCA)
              REAL(KIND=8) :: CHISQ
              REAL(KIND=8) :: ALAMDA
            END SUBROUTINE MRQMIN
          END INTERFACE 
        END MODULE MRQMIN__genmod
