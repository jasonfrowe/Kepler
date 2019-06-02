        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:08 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LMPAR__genmod
          INTERFACE 
            SUBROUTINE LMPAR(N,R,LDR,IPVT,DIAG,QTB,DELTA,PAR,X,SDIAG,WA1&
     &,WA2)
              INTEGER(KIND=4) :: LDR
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: R(LDR,N)
              INTEGER(KIND=4) :: IPVT(N)
              REAL(KIND=8) :: DIAG(N)
              REAL(KIND=8) :: QTB(N)
              REAL(KIND=8) :: DELTA
              REAL(KIND=8) :: PAR
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: SDIAG(N)
              REAL(KIND=8) :: WA1(N)
              REAL(KIND=8) :: WA2(N)
            END SUBROUTINE LMPAR
          END INTERFACE 
        END MODULE LMPAR__genmod
