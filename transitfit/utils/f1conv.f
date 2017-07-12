      !-----------------------------------------------------------------
      ! function  : f1 (complex*16)
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
      !              region, following the ideas in
      !              the paper CPC 138 (2001) 29, section 3.1.3.
      ! 
      ! input     :    a  -> complex parameter of Appell's F1
      !                b1 -> complex parameter of Appell's F1
      !                b2 -> complex parameter of Appell's F1
      !                c  -> complex parameter of Appell's F1
      !                x  -> real variable
      !                y  -> real variable
      !
      !
      ! output    :    f1  -> F1 Appell's hypergeometric function
      !
      !-----------------------------------------------------------------
          complex*16 function f1conv(ca,cb,cbp,cc,x,y)
          implicit none
          real*8 zero
          parameter(zero=1d-10)
          complex*16 ca,cb,cbp,cc,cx,cy
          complex*16 c1,csum,cnum,cden,cterm,c2f1x,c2f1y,ctmp,cz,f21
          complex*16 cgamma,cgammar,f1bnl,hypgeof1
          real*8 abx,aby,one,test,eps,x,y
          integer*4 n,nmax,isneg

          eps = 1d-18
          nmax = 100
          c1 = (1d0,0d0)
          one = 1d0
          cx = c1*x
          cy = c1*y
          abx = dabs(x)
          aby = dabs(y)
          !
          !     Discard problematic cases
          !
          if (dabs(x-1).lt.zero .and. dabs(y-1).lt.zero) then 
                  f1conv =
     .  cgamma(cc)*cgamma(cc-ca-cb-cbp)*(cgammar(cc-ca)*
     .  cgammar(cc-cb-cbp))
                  return
          else if(dabs(x-y).lt.zero) then
                  f1conv = f21(ca,cb+cbp,cc,cx)
                  return
          else if(cdabs(ca-cc).lt.zero) then
                  f1conv = (1-x)**(-cb)*(1-y)**(-cbp)
                  return
          else if(dabs(y-1).lt.zero) then
                  f1conv =
     .  cgamma(cc)*cgamma(cc-ca-cbp)*f21(ca,cb,cc-cbp,cx)*
     .  (cgammar(cc-ca)*cgammar(cc-cbp))
                  return
          else if(dabs(x-1).lt.zero) then
                  f1conv =
     .  cgamma(cc)*cgamma(cc-ca-cb)*f21(ca,cbp,cc-cb,cy)*
     .  (cgammar(cc-ca)*cgammar(cc-cb))
                  return
          else if(abx.lt.zero) then
                  f1conv = f21(ca,cbp,cc,cy)
                  return
          else if(aby.lt.zero) then
                  f1conv = f21(ca,cb,cc,cx)
                  return
          else if (cdabs(cb).lt.zero) then
                  f1conv = f21(ca,cbp,cc,cy)
                  return
          else if (cdabs(cbp).lt.zero) then
                  f1conv = f21(ca,cb,cc,cx)
                  return
          endif  
          if(isneg(cc).lt.0) then
                  Stop 'Stop in f1conv: c is a negative integer'
          endif
    !
    !   Two choices
    !
    !   Ode integration
    !
          if(dabs(x).gt.0.5d0.and.dabs(y).gt.0.5) then
                  f1conv = hypgeof1(ca,cb,cbp,cc,cx,cy)
                  return
          endif
    !
    !   Series Expansion
    !
          if(dabs(x*y).lt.one) then
                  f1conv = f1bnl(ca,cb,cbp,cc,x,y)
                  return
          endif
    !
    !   Sorry
    !           
          write(*,*) 'f1conv  : Not possible',x,y
          return
          end
