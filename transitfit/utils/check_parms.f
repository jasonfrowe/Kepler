      !-----------------------------------------------------------------
      ! function  : check_parms (integer)
      !
      ! package   : F1
      !
      ! Language  : Fortran 90
      !
      ! author    : F. Colavecchia (flavioc@lanl.gov)
      !                            (flavioc@cab.cnea.gov.ar)
      !
      ! date      : 8/20/02      version: 1.0
      !
      ! purpose   : Check the parameters of the F1 function, 
      !             looking for pathological problems. According
      !             to the value of flag (that refers to the analytic
      !             continuation chosen) it tests the parameters
      !             looking for negative integers that 
      !             cause trobles in the Gammas functions of each
      !             analytic continuation.
      !
      ! input     :    a  ->  complex variable
      !                b1  -> complex variable
      !                b2  -> complex variable
      !                c   -> complex variable
      !                flag-> integer variable    
      !
      ! output    :  Integer representing what kind of trouble 
      !              we have:
      !              0 : normal out, no trouble.
      !              1 : 
      !
      ! notes     :  most of the code is f77, few things come from f90.
      !
      !-------------------------------------------------------------------
      integer function check_parms(a,b1,b2,c,flag)
      implicit none
    
      complex*16 a,b1,b2,c
      integer flag,isneg
    
      check_parms = 0
      !
      !   Eq. (23)
      !
      if(flag.eq.23) then
        !
        ! First Term: b2-a
        !
        !
        ! Second Term: c-b2
        !
        if(isneg(c-b2).le.0) then
            check_parms = 1
            return
        end if
      end if    
      !
      !   Eq. (24)
      !
      if(flag.eq.24) then
        !
        ! First Term: b1-a
        !
        !
        ! Second Term: c-b1
        !
        if(isneg(c-b1).le.0) then
            check_parms = 1
            return
        end if
      end if       
      return
      end    
