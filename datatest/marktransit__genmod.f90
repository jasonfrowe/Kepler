        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  9 18:17:14 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MARKTRANSIT__genmod
          INTERFACE 
            SUBROUTINE MARKTRANSIT(NP,NPLANET,NPT,TIME,TFLAG,NFIT,SOL,  &
     &NTT,TOBS,OMC)
              INTEGER(KIND=4) :: NP
              INTEGER(KIND=4) :: NPLANET
              INTEGER(KIND=4) :: NPT
              REAL(KIND=8) :: TIME(:)
              INTEGER(KIND=4) :: TFLAG(:)
              INTEGER(KIND=4) :: NFIT
              REAL(KIND=8) :: SOL(:)
              INTEGER(KIND=4) :: NTT(:)
              REAL(KIND=8) :: TOBS(:,:)
              REAL(KIND=8) :: OMC(:,:)
            END SUBROUTINE MARKTRANSIT
          END INTERFACE 
        END MODULE MARKTRANSIT__genmod
