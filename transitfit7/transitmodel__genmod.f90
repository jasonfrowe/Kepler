        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:07 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TRANSITMODEL__genmod
          INTERFACE 
            SUBROUTINE TRANSITMODEL(NBODIES,NINTG,XPOS,YPOS,ZPOS,SOL,   &
     &TMODEL)
              INTEGER(KIND=4) :: NBODIES
              INTEGER(KIND=4) :: NINTG
              REAL(KIND=8) :: XPOS(:,:)
              REAL(KIND=8) :: YPOS(:,:)
              REAL(KIND=8) :: ZPOS(:,:)
              REAL(KIND=8) :: SOL(:)
              REAL(KIND=8) :: TMODEL
            END SUBROUTINE TRANSITMODEL
          END INTERFACE 
        END MODULE TRANSITMODEL__genmod
