        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:08 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LMDIF1__genmod
          INTERFACE 
            SUBROUTINE LMDIF1(FCN,M,N,X,FVEC,TOL,INFO,IWA,WA,LWA)
              INTEGER(KIND=4) :: LWA
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: M
              EXTERNAL FCN
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: FVEC(M)
              REAL(KIND=8) :: TOL
              INTEGER(KIND=4) :: INFO
              INTEGER(KIND=4) :: IWA(N)
              REAL(KIND=8) :: WA(LWA)
            END SUBROUTINE LMDIF1
          END INTERFACE 
        END MODULE LMDIF1__genmod
