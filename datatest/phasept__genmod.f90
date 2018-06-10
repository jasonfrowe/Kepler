        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  9 18:17:14 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PHASEPT__genmod
          INTERFACE 
            SUBROUTINE PHASEPT(NPT,TIME,PHASE,PERIOD,TOFF)
              INTEGER(KIND=4) :: NPT
              REAL(KIND=8) :: TIME(NPT)
              REAL(KIND=8) :: PHASE(NPT)
              REAL(KIND=8) :: PERIOD
              REAL(KIND=8) :: TOFF
            END SUBROUTINE PHASEPT
          END INTERFACE 
        END MODULE PHASEPT__genmod
