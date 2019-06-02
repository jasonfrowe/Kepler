        !COMPILER-GENERATED INTERFACE MODULE: Sat Jun  1 20:15:08 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE OCTIMING__genmod
          INTERFACE 
            SUBROUTINE OCTIMING(NBODIES,NINTG,XPOS,YPOS,ZPOS,SOL,OPOS,TC&
     &,TCALC,TOLD,BOLD,NTMID,TMID)
              INTEGER(KIND=4) :: NBODIES
              INTEGER(KIND=4) :: NINTG
              REAL(KIND=8) :: XPOS(:,:)
              REAL(KIND=8) :: YPOS(:,:)
              REAL(KIND=8) :: ZPOS(:,:)
              REAL(KIND=8) :: SOL(:)
              REAL(KIND=8) :: OPOS(:)
              INTEGER(KIND=4) :: TC(:)
              REAL(KIND=8) :: TCALC(:)
              REAL(KIND=8) :: TOLD(:,:)
              REAL(KIND=8) :: BOLD(:,:)
              INTEGER(KIND=4) :: NTMID(:)
              REAL(KIND=8) :: TMID(:,:)
            END SUBROUTINE OCTIMING
          END INTERFACE 
        END MODULE OCTIMING__genmod
