        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:09 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MCE_MERG__genmod
          INTERFACE 
            SUBROUTINE MCE_MERG(JCEN,I,J,NBOD,NBIG,M,XH,VH,S,STAT,ELOST)
              INTEGER(KIND=4) :: NBOD
              REAL(KIND=8) :: JCEN(3)
              INTEGER(KIND=4) :: I
              INTEGER(KIND=4) :: J
              INTEGER(KIND=4) :: NBIG
              REAL(KIND=8) :: M(NBOD)
              REAL(KIND=8) :: XH(3,NBOD)
              REAL(KIND=8) :: VH(3,NBOD)
              REAL(KIND=8) :: S(3,NBOD)
              INTEGER(KIND=4) :: STAT(NBOD)
              REAL(KIND=8) :: ELOST
            END SUBROUTINE MCE_MERG
          END INTERFACE 
        END MODULE MCE_MERG__genmod
