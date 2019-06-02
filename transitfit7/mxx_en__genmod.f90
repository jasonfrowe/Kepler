        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:09 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MXX_EN__genmod
          INTERFACE 
            SUBROUTINE MXX_EN(JCEN,NBOD,NBIG,M,XH,VH,S,E,L2)
              INTEGER(KIND=4) :: NBOD
              REAL(KIND=8) :: JCEN(3)
              INTEGER(KIND=4) :: NBIG
              REAL(KIND=8) :: M(NBOD)
              REAL(KIND=8) :: XH(3,NBOD)
              REAL(KIND=8) :: VH(3,NBOD)
              REAL(KIND=8) :: S(3,NBOD)
              REAL(KIND=8) :: E
              REAL(KIND=8) :: L2
            END SUBROUTINE MXX_EN
          END INTERFACE 
        END MODULE MXX_EN__genmod
