        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:09 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PERCORCALC__genmod
          INTERFACE 
            SUBROUTINE PERCORCALC(NBODIES,SOL,NTMIDMAX,NTMID,TMID,PERCOR&
     &)
              INTEGER(KIND=4) :: NBODIES
              REAL(KIND=8) :: SOL(:)
              INTEGER(KIND=4) :: NTMIDMAX
              INTEGER(KIND=4) :: NTMID(:)
              REAL(KIND=8) :: TMID(:,:)
              REAL(KIND=8) :: PERCOR(:)
            END SUBROUTINE PERCORCALC
          END INTERFACE 
        END MODULE PERCORCALC__genmod
