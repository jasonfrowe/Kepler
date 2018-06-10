        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  9 18:17:14 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TRANSITDUR__genmod
          INTERFACE 
            FUNCTION TRANSITDUR(NP,NFIT,SOL)
              INTEGER(KIND=4) :: NFIT
              INTEGER(KIND=4) :: NP
              REAL(KIND=8) :: SOL(NFIT)
              REAL(KIND=8) :: TRANSITDUR
            END FUNCTION TRANSITDUR
          END INTERFACE 
        END MODULE TRANSITDUR__genmod
