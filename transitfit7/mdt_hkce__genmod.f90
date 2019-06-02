        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:09 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MDT_HKCE__genmod
          INTERFACE 
            SUBROUTINE MDT_HKCE(TIME,TSTART,H0,HREC,TOL,RMAX,ELOST,JCEN,&
     &RCEN,NBOD,NBIG,M,X,V,S,RPHY,RCRIT,RCE,STAT,ALGOR,OPT,NGFLAG,      &
     &COLFLAG,CE,NCE,ICE,JCE,NCLO,ICLO,JCLO,DCLO,TCLO,IXVCLO,JXVCLO)
              INTEGER(KIND=4) :: NCE
              INTEGER(KIND=4) :: NBOD
              REAL(KIND=8) :: TIME
              REAL(KIND=8) :: TSTART
              REAL(KIND=8) :: H0
              REAL(KIND=8) :: HREC
              REAL(KIND=8) :: TOL
              REAL(KIND=8) :: RMAX
              REAL(KIND=8) :: ELOST
              REAL(KIND=8) :: JCEN(3)
              REAL(KIND=8) :: RCEN
              INTEGER(KIND=4) :: NBIG
              REAL(KIND=8) :: M(NBOD)
              REAL(KIND=8) :: X(3,NBOD)
              REAL(KIND=8) :: V(3,NBOD)
              REAL(KIND=8) :: S(3,NBOD)
              REAL(KIND=8) :: RPHY(NBOD)
              REAL(KIND=8) :: RCRIT(NBOD)
              REAL(KIND=8) :: RCE(NBOD)
              INTEGER(KIND=4) :: STAT(NBOD)
              INTEGER(KIND=4) :: ALGOR
              INTEGER(KIND=4) :: OPT(8)
              INTEGER(KIND=4) :: NGFLAG
              INTEGER(KIND=4) :: COLFLAG
              INTEGER(KIND=4) :: CE(NBOD)
              INTEGER(KIND=4) :: ICE(NCE)
              INTEGER(KIND=4) :: JCE(NCE)
              INTEGER(KIND=4) :: NCLO
              INTEGER(KIND=4) :: ICLO(50)
              INTEGER(KIND=4) :: JCLO(50)
              REAL(KIND=8) :: DCLO(50)
              REAL(KIND=8) :: TCLO(50)
              REAL(KIND=8) :: IXVCLO(6,50)
              REAL(KIND=8) :: JXVCLO(6,50)
            END SUBROUTINE MDT_HKCE
          END INTERFACE 
        END MODULE MDT_HKCE__genmod
