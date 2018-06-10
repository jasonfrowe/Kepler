        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  9 18:24:41 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GAUSSJ__genmod
          INTERFACE 
            SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
              INTEGER(KIND=4) :: MP
              INTEGER(KIND=4) :: NP
              REAL(KIND=8) :: A(NP,NP)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: B(NP,MP)
              INTEGER(KIND=4) :: M
            END SUBROUTINE GAUSSJ
          END INTERFACE 
        END MODULE GAUSSJ__genmod
