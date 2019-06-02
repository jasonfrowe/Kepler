        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:10 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CALCIMPACT__genmod
          INTERFACE 
            SUBROUTINE CALCIMPACT(NBODIES,Y,SOL,B_CUR)
              INTEGER(KIND=4) :: NBODIES
              REAL(KIND=8) :: Y(:)
              REAL(KIND=8) :: SOL(:)
              REAL(KIND=8) :: B_CUR(:)
            END SUBROUTINE CALCIMPACT
          END INTERFACE 
        END MODULE CALCIMPACT__genmod
