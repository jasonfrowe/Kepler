        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:09 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MFO_OBL__genmod
          INTERFACE 
            SUBROUTINE MFO_OBL(JCEN,NBOD,M,X,A,ACEN)
              INTEGER(KIND=4) :: NBOD
              REAL(KIND=8) :: JCEN(3)
              REAL(KIND=8) :: M(NBOD)
              REAL(KIND=8) :: X(3,NBOD)
              REAL(KIND=8) :: A(3,NBOD)
              REAL(KIND=8) :: ACEN(3)
            END SUBROUTINE MFO_OBL
          END INTERFACE 
        END MODULE MFO_OBL__genmod
