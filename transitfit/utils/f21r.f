!----------------------------------------------------------------------
!     function f21                                              3/2/98
!     ____________
!
!     by:               f. colavecchia
!     flavioc@cab.cnea.gov.ar
!     
!     wrapper for the hypergeometric function 
!     2f1
!
!     input:    
!     complex*16 a,b,c,z
!     integer n
!     
!     output: complex*16 2f1(a,b,c,z)
!
!
!---------------------------------------------------------------------
!

        complex*16 function f21(a,b,c,z)
        implicit none

        integer*4 isneg
        real*8 er,zero
        real*8 re,im
        complex*16 cf,a,b,c,z,cgamma,cgammar,hypgeo
        complex*16 ci,chyp

        ci = (0d0,1d0)
        f21 = (1d0,0d0)
        er = 1.0e-9
        zero = 1e-5*er
!
        if(cdabs(z).lt.zero) return
        if(cdabs(a-c).lt.zero) then
                f21 = (1-z)**(-b)
                return
        end if  
        if(cdabs(b-c).lt.zero) then
                f21 = (1-z)**(-a)
                return
        end if  
!       discard the case z=1
        if(cdabs(z-1d0).lt.zero) then
                f21 = cgamma(c)*cgamma(c-a-b)*cgammar(c-a)*cgammar(c-b)
                return
        end if  
        if(cdabs(a).le.zero.or.cdabs(b).le.zero) then
                f21 = (1d0,0d0)
                return
        else if(isneg(c).eq.-1) then
                stop 'Stop in f21, c is a negative integer'
        end if

!       Our F21
!
!        call cf21d(a,b,c,z,er,cf)
!
!       NR F21
       cf = hypgeo(a,b,c,z)
!       
!       Forrey's F21
!
!       cf = chyp(a,b,c,z)
!
        f21 = cf

        return
        end

        
