        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:09 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MCE_SNIF__genmod
          INTERFACE 
            SUBROUTINE MCE_SNIF(H,I0,NBOD,NBIG,X0,V0,X1,V1,RCRIT,CE,NCE,&
     &ICE,JCE)
              INTEGER(KIND=4) :: NBOD
              REAL(KIND=8) :: H
              INTEGER(KIND=4) :: I0
              INTEGER(KIND=4) :: NBIG
              REAL(KIND=8) :: X0(3,NBOD)
              REAL(KIND=8) :: V0(3,NBOD)
              REAL(KIND=8) :: X1(3,NBOD)
              REAL(KIND=8) :: V1(3,NBOD)
              REAL(KIND=8) :: RCRIT(NBOD)
              INTEGER(KIND=4) :: CE(NBOD)
              INTEGER(KIND=4) :: NCE
              INTEGER(KIND=4) :: ICE(2000)
              INTEGER(KIND=4) :: JCE(2000)
            END SUBROUTINE MCE_SNIF
          END INTERFACE 
        END MODULE MCE_SNIF__genmod
