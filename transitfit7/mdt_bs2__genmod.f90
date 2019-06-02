        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:09 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MDT_BS2__genmod
          INTERFACE 
            SUBROUTINE MDT_BS2(TIME,H0,HDID,TOL,JCEN,NBOD,NBIG,MASS,X0, &
     &V0,S,RPHYS,RCRIT,NGF,STAT,DTFLAG,NGFLAG,OPT,NCE,ICE,JCE)
              INTEGER(KIND=4) :: NCE
              INTEGER(KIND=4) :: NBOD
              REAL(KIND=8) :: TIME
              REAL(KIND=8) :: H0
              REAL(KIND=8) :: HDID
              REAL(KIND=8) :: TOL
              REAL(KIND=8) :: JCEN(3)
              INTEGER(KIND=4) :: NBIG
              REAL(KIND=8) :: MASS(NBOD)
              REAL(KIND=8) :: X0(3,NBOD)
              REAL(KIND=8) :: V0(3,NBOD)
              REAL(KIND=8) :: S(3,NBOD)
              REAL(KIND=8) :: RPHYS(NBOD)
              REAL(KIND=8) :: RCRIT(NBOD)
              REAL(KIND=8) :: NGF(4,NBOD)
              INTEGER(KIND=4) :: STAT(NBOD)
              INTEGER(KIND=4) :: DTFLAG
              INTEGER(KIND=4) :: NGFLAG
              INTEGER(KIND=4) :: OPT(8)
              INTEGER(KIND=4) :: ICE(NCE)
              INTEGER(KIND=4) :: JCE(NCE)
            END SUBROUTINE MDT_BS2
          END INTERFACE 
        END MODULE MDT_BS2__genmod
