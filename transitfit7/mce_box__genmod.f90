        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:09 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MCE_BOX__genmod
          INTERFACE 
            SUBROUTINE MCE_BOX(NBOD,H,X0,V0,X1,V1,XMIN,XMAX,YMIN,YMAX)
              INTEGER(KIND=4) :: NBOD
              REAL(KIND=8) :: H
              REAL(KIND=8) :: X0(3,NBOD)
              REAL(KIND=8) :: V0(3,NBOD)
              REAL(KIND=8) :: X1(3,NBOD)
              REAL(KIND=8) :: V1(3,NBOD)
              REAL(KIND=8) :: XMIN(NBOD)
              REAL(KIND=8) :: XMAX(NBOD)
              REAL(KIND=8) :: YMIN(NBOD)
              REAL(KIND=8) :: YMAX(NBOD)
            END SUBROUTINE MCE_BOX
          END INTERFACE 
        END MODULE MCE_BOX__genmod
