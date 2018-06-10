        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  9 18:17:15 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FUNCS__genmod
          INTERFACE 
            SUBROUTINE FUNCS(X,A,Y,DYDA,NA,N,DTYPE)
              INTEGER(KIND=4) :: NA
              REAL(KIND=8) :: X
              REAL(KIND=8) :: A(NA)
              REAL(KIND=8) :: Y
              REAL(KIND=8) :: DYDA(NA)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: DTYPE
            END SUBROUTINE FUNCS
          END INTERFACE 
        END MODULE FUNCS__genmod
