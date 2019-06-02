        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:09 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MCE_INIT__genmod
          INTERFACE 
            SUBROUTINE MCE_INIT(TSTART,ALGOR,H,JCEN,RCEN,RMAX,CEFAC,NBOD&
     &,NBIG,M,X,V,S,RHO,RCEH,RPHYS,RCE,RCRIT,RCRITFLAG)
              INTEGER(KIND=4) :: NBOD
              REAL(KIND=8) :: TSTART
              INTEGER(KIND=4) :: ALGOR
              REAL(KIND=8) :: H
              REAL(KIND=8) :: JCEN(3)
              REAL(KIND=8) :: RCEN
              REAL(KIND=8) :: RMAX
              REAL(KIND=8) :: CEFAC
              INTEGER(KIND=4) :: NBIG
              REAL(KIND=8) :: M(NBOD)
              REAL(KIND=8) :: X(3,NBOD)
              REAL(KIND=8) :: V(3,NBOD)
              REAL(KIND=8) :: S(3,NBOD)
              REAL(KIND=8) :: RHO(NBOD)
              REAL(KIND=8) :: RCEH(NBOD)
              REAL(KIND=8) :: RPHYS(NBOD)
              REAL(KIND=8) :: RCE(NBOD)
              REAL(KIND=8) :: RCRIT(NBOD)
              INTEGER(KIND=4) :: RCRITFLAG
            END SUBROUTINE MCE_INIT
          END INTERFACE 
        END MODULE MCE_INIT__genmod
