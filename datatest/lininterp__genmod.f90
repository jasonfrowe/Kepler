        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  9 18:17:14 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LININTERP__genmod
          INTERFACE 
            SUBROUTINE LININTERP(X,Y,NPMAX,NMAX,NP,NPT,XIN,YOUT)
              INTEGER(KIND=4) :: NMAX
              INTEGER(KIND=4) :: NPMAX
              REAL(KIND=8) :: X(NPMAX,NMAX)
              REAL(KIND=8) :: Y(NPMAX,NMAX)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: NPT(NPMAX)
              REAL(KIND=8) :: XIN
              REAL(KIND=8) :: YOUT
            END SUBROUTINE LININTERP
          END INTERFACE 
        END MODULE LININTERP__genmod
