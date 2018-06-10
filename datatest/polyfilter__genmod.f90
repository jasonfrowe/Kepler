        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  9 18:17:14 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE POLYFILTER__genmod
          INTERFACE 
            SUBROUTINE POLYFILTER(NPT,TIME,FLUX,FERR,TFLAG,BOXBIN,NFITP)
              INTEGER(KIND=4) :: NPT
              REAL(KIND=8) :: TIME(:)
              REAL(KIND=8) :: FLUX(:)
              REAL(KIND=8) :: FERR(:)
              INTEGER(KIND=4) :: TFLAG(:)
              REAL(KIND=8) :: BOXBIN
              INTEGER(KIND=4) :: NFITP
            END SUBROUTINE POLYFILTER
          END INTERFACE 
        END MODULE POLYFILTER__genmod
