        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  9 18:17:14 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READTTFILE__genmod
          INTERFACE 
            SUBROUTINE READTTFILE(NUNIT,NPLANETMAX,NMAX,NPLANET,NTT,TOBS&
     &,OMC)
              INTEGER(KIND=4) :: NPLANET
              INTEGER(KIND=4) :: NMAX
              INTEGER(KIND=4) :: NPLANETMAX
              INTEGER(KIND=4) :: NUNIT
              INTEGER(KIND=4) :: NTT(NPLANET)
              REAL(KIND=8) :: TOBS(NPLANETMAX,NMAX)
              REAL(KIND=8) :: OMC(NPLANETMAX,NMAX)
            END SUBROUTINE READTTFILE
          END INTERFACE 
        END MODULE READTTFILE__genmod
