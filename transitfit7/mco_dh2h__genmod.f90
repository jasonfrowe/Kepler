        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:09 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MCO_DH2H__genmod
          INTERFACE 
            SUBROUTINE MCO_DH2H(TIME,JCEN,NBOD,NBIG,H,M,X,V,XH,VH)
              INTEGER(KIND=4) :: NBOD
              REAL(KIND=8) :: TIME
              REAL(KIND=8) :: JCEN(3)
              INTEGER(KIND=4) :: NBIG
              REAL(KIND=8) :: H
              REAL(KIND=8) :: M(NBOD)
              REAL(KIND=8) :: X(3,NBOD)
              REAL(KIND=8) :: V(3,NBOD)
              REAL(KIND=8) :: XH(3,NBOD)
              REAL(KIND=8) :: VH(3,NBOD)
            END SUBROUTINE MCO_DH2H
          END INTERFACE 
        END MODULE MCO_DH2H__genmod
