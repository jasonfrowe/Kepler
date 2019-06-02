        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:08 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LMDIF__genmod
          INTERFACE 
            SUBROUTINE LMDIF(FCN,M,N,X,FVEC,FTOL,XTOL,GTOL,MAXFEV,EPSFCN&
     &,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,FJAC,LDFJAC,IPVT,QTF,WA1,WA2,  &
     &WA3,WA4)
              INTEGER(KIND=4) :: LDFJAC
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: M
              EXTERNAL FCN
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: FVEC(M)
              REAL(KIND=8) :: FTOL
              REAL(KIND=8) :: XTOL
              REAL(KIND=8) :: GTOL
              INTEGER(KIND=4) :: MAXFEV
              REAL(KIND=8) :: EPSFCN
              REAL(KIND=8) :: DIAG(N)
              INTEGER(KIND=4) :: MODE
              REAL(KIND=8) :: FACTOR
              INTEGER(KIND=4) :: NPRINT
              INTEGER(KIND=4) :: INFO
              INTEGER(KIND=4) :: NFEV
              REAL(KIND=8) :: FJAC(LDFJAC,N)
              INTEGER(KIND=4) :: IPVT(N)
              REAL(KIND=8) :: QTF(N)
              REAL(KIND=8) :: WA1(N)
              REAL(KIND=8) :: WA2(N)
              REAL(KIND=8) :: WA3(N)
              REAL(KIND=8) :: WA4(M)
            END SUBROUTINE LMDIF
          END INTERFACE 
        END MODULE LMDIF__genmod
