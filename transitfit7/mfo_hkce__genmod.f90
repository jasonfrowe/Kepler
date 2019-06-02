        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:09 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MFO_HKCE__genmod
          INTERFACE 
            SUBROUTINE MFO_HKCE(TIME,JCEN,NBOD,NBIG,M,X,V,SPIN,RCRIT,A, &
     &STAT,NGF,NGFLAG,OPT,NCE,ICE,JCE)
              INTEGER(KIND=4) :: NCE
              INTEGER(KIND=4) :: NBOD
              REAL(KIND=8) :: TIME
              REAL(KIND=8) :: JCEN(3)
              INTEGER(KIND=4) :: NBIG
              REAL(KIND=8) :: M(NBOD)
              REAL(KIND=8) :: X(3,NBOD)
              REAL(KIND=8) :: V(3,NBOD)
              REAL(KIND=8) :: SPIN(3,NBOD)
              REAL(KIND=8) :: RCRIT(NBOD)
              REAL(KIND=8) :: A(3,NBOD)
              INTEGER(KIND=4) :: STAT(NBOD)
              REAL(KIND=8) :: NGF(4,NBOD)
              INTEGER(KIND=4) :: NGFLAG
              INTEGER(KIND=4) :: OPT(8)
              INTEGER(KIND=4) :: ICE(NCE)
              INTEGER(KIND=4) :: JCE(NCE)
            END SUBROUTINE MFO_HKCE
          END INTERFACE 
        END MODULE MFO_HKCE__genmod
