        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:10 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FCN__genmod
          INTERFACE 
            SUBROUTINE FCN(NPT,N,X,FVEC,IFLAG)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NPT
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: FVEC(NPT)
              INTEGER(KIND=4) :: IFLAG
            END SUBROUTINE FCN
          END INTERFACE 
        END MODULE FCN__genmod
