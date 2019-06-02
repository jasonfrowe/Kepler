        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:05 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READKEPLC__genmod
          INTERFACE 
            SUBROUTINE READKEPLC(PHOTFILE,NPT,TIME,FLUX,FERR,ITIME)
              CHARACTER(LEN=80) :: PHOTFILE
              INTEGER(KIND=4) :: NPT
              REAL(KIND=8) :: TIME(:)
              REAL(KIND=8) :: FLUX(:)
              REAL(KIND=8) :: FERR(:)
              REAL(KIND=8) :: ITIME(:)
            END SUBROUTINE READKEPLC
          END INTERFACE 
        END MODULE READKEPLC__genmod
