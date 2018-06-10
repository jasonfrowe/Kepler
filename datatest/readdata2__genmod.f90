        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  9 18:17:13 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READDATA2__genmod
          INTERFACE 
            SUBROUTINE READDATA2(OBSFILE,NMAX,NPT,TIME,FLUX,FERR,ITIME, &
     &ZTIME)
              CHARACTER(LEN=80) :: OBSFILE
              INTEGER(KIND=4) :: NMAX
              INTEGER(KIND=4) :: NPT
              REAL(KIND=8) :: TIME(:)
              REAL(KIND=8) :: FLUX(:)
              REAL(KIND=8) :: FERR(:)
              REAL(KIND=8) :: ITIME(:)
              REAL(KIND=8) :: ZTIME
            END SUBROUTINE READDATA2
          END INTERFACE 
        END MODULE READDATA2__genmod
