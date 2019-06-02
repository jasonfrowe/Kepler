        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:10 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FITTRANSITMODEL__genmod
          INTERFACE 
            SUBROUTINE FITTRANSITMODEL(NBODIES,NPT,TOL,SOL,SERR,TIME,   &
     &FLUX,FERR,ITIME,NTMIDMAX)
              INTEGER(KIND=4) ,TARGET :: NBODIES
              INTEGER(KIND=4) :: NPT
              REAL(KIND=8) ,TARGET :: TOL
              REAL(KIND=8) ,TARGET :: SOL(:)
              REAL(KIND=8) ,TARGET :: SERR(:,:)
              REAL(KIND=8) ,TARGET :: TIME(:)
              REAL(KIND=8) ,TARGET :: FLUX(:)
              REAL(KIND=8) ,TARGET :: FERR(:)
              REAL(KIND=8) ,TARGET :: ITIME(:)
              INTEGER(KIND=4) ,TARGET :: NTMIDMAX
            END SUBROUTINE FITTRANSITMODEL
          END INTERFACE 
        END MODULE FITTRANSITMODEL__genmod
