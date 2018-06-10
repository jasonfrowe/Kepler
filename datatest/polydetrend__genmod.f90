        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  9 18:17:14 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE POLYDETREND__genmod
          INTERFACE 
            SUBROUTINE POLYDETREND(NPT,TIME,MAG,MERR,NFIT,TZERO,OFF)
              INTEGER(KIND=4) :: NPT
              REAL(KIND=8) :: TIME(:)
              REAL(KIND=8) :: MAG(:)
              REAL(KIND=8) :: MERR(:)
              INTEGER(KIND=4) :: NFIT
              REAL(KIND=8) :: TZERO
              REAL(KIND=8) :: OFF
            END SUBROUTINE POLYDETREND
          END INTERFACE 
        END MODULE POLYDETREND__genmod
