        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:08 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE OCCULTNL__genmod
          INTERFACE 
            SUBROUTINE OCCULTNL(RL,C1,C2,C3,C4,B0,MULIMB0,MULIMBF,NB)
              INTEGER(KIND=4) :: NB
              REAL(KIND=8) :: RL
              REAL(KIND=8) :: C1
              REAL(KIND=8) :: C2
              REAL(KIND=8) :: C3
              REAL(KIND=8) :: C4
              REAL(KIND=8) :: B0(NB)
              REAL(KIND=8) :: MULIMB0(NB)
              REAL(KIND=8) :: MULIMBF(NB,5)
            END SUBROUTINE OCCULTNL
          END INTERFACE 
        END MODULE OCCULTNL__genmod
