        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:10 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FCNB__genmod
          INTERFACE 
            SUBROUTINE FCNB(NPT,N,X,FVEC,IFLAG)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NPT
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: FVEC(NPT)
              INTEGER(KIND=4) :: IFLAG
            END SUBROUTINE FCNB
          END INTERFACE 
        END MODULE FCNB__genmod
