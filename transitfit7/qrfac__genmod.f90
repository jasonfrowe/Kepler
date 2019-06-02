        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:08 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE QRFAC__genmod
          INTERFACE 
            SUBROUTINE QRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,RDIAG,ACNORM,WA)
              INTEGER(KIND=4) :: LIPVT
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: A(LDA,N)
              LOGICAL(KIND=4) :: PIVOT
              INTEGER(KIND=4) :: IPVT(LIPVT)
              REAL(KIND=8) :: RDIAG(N)
              REAL(KIND=8) :: ACNORM(N)
              REAL(KIND=8) :: WA(N)
            END SUBROUTINE QRFAC
          END INTERFACE 
        END MODULE QRFAC__genmod
