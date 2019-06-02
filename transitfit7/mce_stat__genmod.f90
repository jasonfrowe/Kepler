        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:09 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MCE_STAT__genmod
          INTERFACE 
            SUBROUTINE MCE_STAT(TIME,H,RCEN,NBOD,NBIG,M,X0,V0,X1,V1,RCE,&
     &RPHYS,NCLO,ICLO,JCLO,DCLO,TCLO,IXVCLO,JXVCLO,NHIT,IHIT,JHIT,CHIT, &
     &DHIT,THIT,THIT1,NOWFLAG,STAT)
              INTEGER(KIND=4) :: NBOD
              REAL(KIND=8) :: TIME
              REAL(KIND=8) :: H
              REAL(KIND=8) :: RCEN
              INTEGER(KIND=4) :: NBIG
              REAL(KIND=8) :: M(NBOD)
              REAL(KIND=8) :: X0(3,NBOD)
              REAL(KIND=8) :: V0(3,NBOD)
              REAL(KIND=8) :: X1(3,NBOD)
              REAL(KIND=8) :: V1(3,NBOD)
              REAL(KIND=8) :: RCE(NBOD)
              REAL(KIND=8) :: RPHYS(NBOD)
              INTEGER(KIND=4) :: NCLO
              INTEGER(KIND=4) :: ICLO(50)
              INTEGER(KIND=4) :: JCLO(50)
              REAL(KIND=8) :: DCLO(50)
              REAL(KIND=8) :: TCLO(50)
              REAL(KIND=8) :: IXVCLO(6,50)
              REAL(KIND=8) :: JXVCLO(6,50)
              INTEGER(KIND=4) :: NHIT
              INTEGER(KIND=4) :: IHIT(50)
              INTEGER(KIND=4) :: JHIT(50)
              INTEGER(KIND=4) :: CHIT(50)
              REAL(KIND=8) :: DHIT(50)
              REAL(KIND=8) :: THIT(50)
              REAL(KIND=8) :: THIT1
              INTEGER(KIND=4) :: NOWFLAG
              INTEGER(KIND=4) :: STAT(NBOD)
            END SUBROUTINE MCE_STAT
          END INTERFACE 
        END MODULE MCE_STAT__genmod
