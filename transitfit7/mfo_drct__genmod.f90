        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:09 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MFO_DRCT__genmod
          INTERFACE 
            SUBROUTINE MFO_DRCT(I0,NBOD,NBIG,M,X,RCRIT,A,STAT)
              INTEGER(KIND=4) :: NBOD
              INTEGER(KIND=4) :: I0
              INTEGER(KIND=4) :: NBIG
              REAL(KIND=8) :: M(NBOD)
              REAL(KIND=8) :: X(3,NBOD)
              REAL(KIND=8) :: RCRIT(NBOD)
              REAL(KIND=8) :: A(3,NBOD)
              INTEGER(KIND=4) :: STAT(NBOD)
            END SUBROUTINE MFO_DRCT
          END INTERFACE 
        END MODULE MFO_DRCT__genmod
