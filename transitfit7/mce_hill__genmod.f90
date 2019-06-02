        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:09 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MCE_HILL__genmod
          INTERFACE 
            SUBROUTINE MCE_HILL(NBOD,M,X,V,HILL,A)
              INTEGER(KIND=4) :: NBOD
              REAL(KIND=8) :: M(NBOD)
              REAL(KIND=8) :: X(3,NBOD)
              REAL(KIND=8) :: V(3,NBOD)
              REAL(KIND=8) :: HILL(NBOD)
              REAL(KIND=8) :: A(NBOD)
            END SUBROUTINE MCE_HILL
          END INTERFACE 
        END MODULE MCE_HILL__genmod
