        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:06 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LCMODEL__genmod
          INTERFACE 
            SUBROUTINE LCMODEL(NBODIES,NPT,TOL,SOL,TIME,ITIME,NTMID,TMID&
     &,PERCOR,ANS,COLFLAG,ITPRINT,ITMODEL)
              INTEGER(KIND=4) :: NBODIES
              INTEGER(KIND=4) :: NPT
              REAL(KIND=8) :: TOL
              REAL(KIND=8) :: SOL(:)
              REAL(KIND=8) :: TIME(:)
              REAL(KIND=8) :: ITIME(:)
              INTEGER(KIND=4) :: NTMID(:)
              REAL(KIND=8) :: TMID(:,:)
              REAL(KIND=8) :: PERCOR(:)
              REAL(KIND=8) :: ANS(:)
              INTEGER(KIND=4) :: COLFLAG
              INTEGER(KIND=4) :: ITPRINT
              INTEGER(KIND=4) :: ITMODEL
            END SUBROUTINE LCMODEL
          END INTERFACE 
        END MODULE LCMODEL__genmod
