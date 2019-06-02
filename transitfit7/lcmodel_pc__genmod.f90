        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:10 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LCMODEL_PC__genmod
          INTERFACE 
            SUBROUTINE LCMODEL_PC(NBODIES,NPT,TOL,SOL,TIME,NTMID,TMID,  &
     &PERCOR,COLFLAG,ITPRINT)
              INTEGER(KIND=4) :: NBODIES
              INTEGER(KIND=4) :: NPT
              REAL(KIND=8) :: TOL
              REAL(KIND=8) :: SOL(:)
              REAL(KIND=8) :: TIME(:)
              INTEGER(KIND=4) :: NTMID(:)
              REAL(KIND=8) :: TMID(:,:)
              REAL(KIND=8) :: PERCOR(:)
              INTEGER(KIND=4) :: COLFLAG
              INTEGER(KIND=4) :: ITPRINT
            END SUBROUTINE LCMODEL_PC
          END INTERFACE 
        END MODULE LCMODEL_PC__genmod
