        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:07 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MANDELAGOL__genmod
          INTERFACE 
            FUNCTION MANDELAGOL(NINTG,R1,R2,X1,X2,Y1,Y2,C,B0,MU,MULIMB0,&
     &MULIMBF,DIST)
              INTEGER(KIND=4) :: NINTG
              REAL(KIND=8) :: R1
              REAL(KIND=8) :: R2
              REAL(KIND=8) :: X1
              REAL(KIND=8) :: X2(NINTG)
              REAL(KIND=8) :: Y1
              REAL(KIND=8) :: Y2(NINTG)
              REAL(KIND=8) :: C(4)
              REAL(KIND=8) :: B0(NINTG)
              REAL(KIND=8) :: MU(NINTG)
              REAL(KIND=8) :: MULIMB0(NINTG)
              REAL(KIND=8) :: MULIMBF(NINTG,5)
              REAL(KIND=8) :: DIST(NINTG)
              REAL(KIND=8) :: MANDELAGOL
            END FUNCTION MANDELAGOL
          END INTERFACE 
        END MODULE MANDELAGOL__genmod
