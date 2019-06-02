        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:05 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READINPUTSOL__genmod
          INTERFACE 
            SUBROUTINE READINPUTSOL(INPUTSOL,SOL,SERR)
              CHARACTER(LEN=80) :: INPUTSOL
              REAL(KIND=8) :: SOL(:)
              REAL(KIND=8) :: SERR(:,:)
            END SUBROUTINE READINPUTSOL
          END INTERFACE 
        END MODULE READINPUTSOL__genmod
