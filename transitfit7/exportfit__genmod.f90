        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:08 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EXPORTFIT__genmod
          INTERFACE 
            SUBROUTINE EXPORTFIT(NBODIES,SOL,SERR)
              INTEGER(KIND=4) :: NBODIES
              REAL(KIND=8) :: SOL(:)
              REAL(KIND=8) :: SERR(:,:)
            END SUBROUTINE EXPORTFIT
          END INTERFACE 
        END MODULE EXPORTFIT__genmod
