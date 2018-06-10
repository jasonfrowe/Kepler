        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  9 18:17:13 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GETFITPARS__genmod
          INTERFACE 
            SUBROUTINE GETFITPARS(NUNIT,NFIT,NPLANET,SOL,SERR,ERR)
              INTEGER(KIND=4) :: NFIT
              INTEGER(KIND=4) :: NUNIT
              INTEGER(KIND=4) :: NPLANET
              REAL(KIND=8) :: SOL(NFIT)
              REAL(KIND=8) :: SERR(NFIT,2)
              REAL(KIND=8) :: ERR(NFIT)
            END SUBROUTINE GETFITPARS
          END INTERFACE 
        END MODULE GETFITPARS__genmod
