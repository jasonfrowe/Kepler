        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:10 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CALCBMIN__genmod
          INTERFACE 
            SUBROUTINE CALCBMIN(NBC,NBODIES,T,SOL,TOL,NBOD,M,X,V,ALGOR, &
     &NBIG,NGFLAG,OPFLAG,COLFLAG,OPT,STAT,RCEN,RMAX,TSTART,JCEN,EN,AM,  &
     &RPHYS,RCE,RCRIT,S,A,HREC,ANGF,AUSR,TTRAN)
              INTEGER(KIND=4) ,TARGET :: NBC
              INTEGER(KIND=4) ,TARGET :: NBODIES
              REAL(KIND=8) ,TARGET :: T
              REAL(KIND=8) ,TARGET :: SOL(:)
              REAL(KIND=8) ,TARGET :: TOL
              INTEGER(KIND=4) ,TARGET :: NBOD
              REAL(KIND=8) ,TARGET :: M(:)
              REAL(KIND=8) ,TARGET :: X(:,:)
              REAL(KIND=8) ,TARGET :: V(:,:)
              INTEGER(KIND=4) ,TARGET :: ALGOR
              INTEGER(KIND=4) ,TARGET :: NBIG
              INTEGER(KIND=4) ,TARGET :: NGFLAG
              INTEGER(KIND=4) ,TARGET :: OPFLAG
              INTEGER(KIND=4) ,TARGET :: COLFLAG
              INTEGER(KIND=4) ,TARGET :: OPT(8)
              INTEGER(KIND=4) ,TARGET :: STAT(:)
              REAL(KIND=8) ,TARGET :: RCEN
              REAL(KIND=8) ,TARGET :: RMAX
              REAL(KIND=8) ,TARGET :: TSTART
              REAL(KIND=8) ,TARGET :: JCEN(3)
              REAL(KIND=8) ,TARGET :: EN(3)
              REAL(KIND=8) ,TARGET :: AM(3)
              REAL(KIND=8) ,TARGET :: RPHYS(:)
              REAL(KIND=8) ,TARGET :: RCE(:)
              REAL(KIND=8) ,TARGET :: RCRIT(:)
              REAL(KIND=8) ,TARGET :: S(:,:)
              REAL(KIND=8) ,TARGET :: A(3,2000)
              REAL(KIND=8) ,TARGET :: HREC
              REAL(KIND=8) ,TARGET :: ANGF(3,2000)
              REAL(KIND=8) ,TARGET :: AUSR(3,2000)
              REAL(KIND=8) :: TTRAN
            END SUBROUTINE CALCBMIN
          END INTERFACE 
        END MODULE CALCBMIN__genmod
