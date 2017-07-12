      !-------------------------------------------------------------------
      ! function  : f1bnl (complex*16)
      !
      ! package   : F1
      !
      ! Language  : Fortran 90
      !
      ! author    : F. Colavecchia (flavioc@lanl.gov)
      !                            (flavioc@cab.cnea.gov.ar)
      !
      ! date      : 3/26/97      version: 0.1
      ! revision  : 6/25/02      version: 1.0
      !
      ! purpose   :  Computes the F1 function in the convergence
      !              region by the series of Burchnall and Chaundy,
      !              Eq. (31) of paper CPC 138 (2001) 29.
      !
      ! input     :    a  -> complex parameter of Appell's F1
      !                b1 -> complex parameter of Appell's F1
      !                b2 -> complex parameter of Appell's F1
      !                c  -> complex parameter of Appell's F1
      !                x  -> real variable
      !                y  -> real variable
      !
      !
      ! output    :    f1bnl  -> F1 Appell's hypergeometric function
      !
      ! modules   :
      !
      !
      ! common    :
      !
      !
      ! notes     :  most of the code is f77, few things come from f90.
      !
      !-----------------------------------------------------------------
      complex*16 function f1bnl(ca,cb,cbp,cc,x,y)

      implicit none
      real*8 zero
      parameter(zero=1d-12)
      complex*16 ca,cb,cbp,cc,cx,cy
      complex*16 c1,csum,cnum,cden,cterm,c2f1x,c2f1y,ctmp,cz,f21
      complex*16 cgamma,cgammar
      real*8 abx,aby,one,test,eps,x,y
      integer*4 n,nmax,isneg

      eps = 1d-18
      nmax = 300
      c1 = (1d0,0d0)
      one = 1d0
      !   
      !   Check whether x*y are less than one
      !
      cx = c1*x
      cy = c1*y
      abx = dabs(x)
      aby = dabs(y)
      if(dabs(x*y).gt.one) then
            write(*,*) 'Error in function f1bnl',abx,aby
      endif
      !
      !   Discard problematic cases
      !
            if (x.eq.1 .and. y.eq.1) then 
                    f1bnl =cgamma(cc)*cgamma(cc-ca-cb-cbp)*
     .                  (cgammar(cc-ca)*cgammar(cc-cb-cbp))
                    return
            else if(abs(x-y).lt.zero) then
                    f1bnl = f21(ca,cb+cbp,cc,cx)
                    return
            else if(cdabs(ca-cc).lt.zero) then
                    f1bnl = (1-x)**(-cb)*(1-y)**(-cbp)
                    return
            else if(abs(y-1).lt.zero) then
                    f1bnl =cgamma(cc)*cgamma(cc-ca-cbp)*
     .             f21(ca,cb,cc-cbp,cx)*(cgammar(cc-ca)*cgammar(cc-cbp))
                    return
            else if(abs(x-1).lt.zero) then
                    f1bnl =cgamma(cc)*cgamma(cc-ca-cb)*
     .              f21(ca,cbp,cc-cb,cy)*(cgammar(cc-ca)*cgammar(cc-cb))
                    return
            else if(abx.lt.zero) then
                    f1bnl = f21(ca,cbp,cc,cy)
                    return
            else if(aby.lt.zero) then
                    f1bnl = f21(ca,cb,cc,cx)
                    return
            else if (cdabs(cb).lt.zero) then
                    f1bnl = f21(ca,cbp,cc,cy)
                    return
            else if (cdabs(cbp).lt.zero) then
                    f1bnl = f21(ca,cb,cc,cx)
                    return
            endif  
      if(isneg(cc).lt.0) then
            Stop 'Stop in f1bnl: c is a negative integer'
      endif

      cz = cx*cy
      ctmp = c1
      !
      !   Zero order of the series
      !
      c2f1x = f21(ca,cb,cc,cx)
      c2f1y = f21(ca,cbp,cc,cy)
      csum = c2f1x*c2f1y
      !
      !   First order of the series
      !
      cnum = ca*cb*cbp*(cc-ca)
      cden = cc*cc*(cc+1)
      c2f1x = f21(ca+1,cb+1,cc+2,cx)
      c2f1y = f21(ca+1,cbp+1,cc+2,cy)
      ctmp = cnum*cz/cden
      cterm= ctmp*c2f1x*c2f1y
      csum = csum+cterm
      test = cdabs(cterm/csum)
      if(test.lt.eps)then
                    f1bnl = csum
      !           write(*,*) "test", test
                    return
      endif
      !
      !   Recurrence relation for n>2
      !
      do n=2,nmax
            cnum = (ca+n-1)*(cb+n-1)*(cbp+n-1)*(cc-ca+n-1)
            cden =(cc+2*n-2)*(cc+2*n-3)*(cc+2*n-1)*(cc+2*n-2)*n/(cc+n-2)
            c2f1x = f21(ca+n,cb+n,cc+2*n,cx)
            c2f1y = f21(ca+n,cbp+n,cc+2*n,cy)
            ctmp = ctmp*cnum*cz/cden
            cterm= ctmp*c2f1x*c2f1y
            csum = csum+cterm
            test = cdabs(cterm/csum)
            if(test.lt.eps)then
                    f1bnl = csum
      !           write(*,*) n,test
                    return
            endif
      end do
      !
      !   Series do not converge
      !
      write(*,*) 'Series did not converge'
      call writef1(6,ca,cb,cbp,cc,cx,cy,csum)
      return

      end
      !-----------------------------------------------------------------
      ! function  : cf1bnl (complex*16)
      !
      ! package   : F1
      !
      ! Language  : Fortran 90
      !
      ! author    : F. Colavecchia (flavioc@lanl.gov)
      !                            (flavioc@cab.cnea.gov.ar)
      !
      ! date      : 3/26/97      version: 0.1
      ! revision  : 6/25/02      version: 1.0
      !
      ! purpose   :  Computes the F1 function in the convergence
      !              region by the series of Burchnall and Chaundy,
      !              Eq. (31) of paper CPC 138 (2001) 29. SEE NOTES.
      !
      ! input     :    a  -> complex parameter of Appell's F1
      !                b1 -> complex parameter of Appell's F1
      !                b2 -> complex parameter of Appell's F1
      !                c  -> complex parameter of Appell's F1
      !                x  -> complex variable
      !                y  -> complex variable
      !
      !
      ! output    :    cf1bnl  -> F1 Appell's hypergeometric function
      !
      ! modules   :
      !
      !
      ! common    :
      !
      !
      ! notes     :  most of the code is f77, few things come from f90.
      !                            Actually, this only masks the f1bnl for
      !                            complex
      !                            values, BUT since f1bnl manages REAL
      !                            variables,
      !                           DO NOT use with full complex ones.    
      !
      !-----------------------------------------------------------------





      !-----------------------------------------------------------------
      !   function cf1bnl
      !   2/2/97
      !   ______________
      !
      !   by:             f. colavecchia
      !           
      !                   calculates the hypergeometric function 
      !                   f1(a,b,b',c,x,y) with formula OF Burchnall &
      !                   Chaundy for |x| and |y| <1
      !
      !   input:  
      !                   complex*16 ca,cb,cbp,cc,cx,cy
      !                   
      !   output: complex*16 cf1bnl
      !
      !   limitations: 

      !
      !
      !
      !-----------------------------------------------------------------
      !
      complex*16 function cf1bnl(ca,cb,cbp,cc,cx,cy)

      implicit none

      complex*16 ca,cb,cbp,cc,cx,cy
      complex*16 c1,csum,cnum,cden,cterm,c2f1x,c2f1y,ctmp,cz,f21,f2
      complex*16 cgamma,cgammar,f1bnl
      real*8 abx,aby,one,test,eps,half
      integer*4 n,nmax

      eps = 1d-12
      nmax = 100
      c1 = (1d0,0d0)
      one = 1d0
      half = 0.5

      cf1bnl = f1bnl(ca,cb,cbp,cc,dreal(cx),dreal(cy))

      return
      end

