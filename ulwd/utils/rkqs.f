CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n,NMAX
      REAL*8 eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
C     EXTERNAL derivs
      PARAMETER (NMAX=50) !Maximum number of equations.
C     USES derivs,rkck
C     Fifth-order Runge-Kutta step with monitoring of local truncation 
C     error to ensure accuracy and adjust stepsize. Input are the 
C     dependent variable vector y(1:n) and its derivative dydx(1:n) at 
C     the starting value of the independent variable x. Also input are 
C     the stepsize to be attempted htry, the required accuracy eps, and 
C     the vector yscal(1:n) against which the error is scaled. On 
C     output, y and x are replaced by their new values, hdid is the
C     stepsize that was actually accomplished, and hnext is the 
C     estimated next stepsize. derivs is the user-supplied subroutine 
C     that computes the right-hand side derivatives.
      INTEGER i
      REAL*8 errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,
     .  PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89d-4)
C     The value ERRCON equals (5/SAFETY)**(1/PGROW), see use below.
      h=htry !Set stepsize to the initial trial value.
 1    call rkck(y,dydx,n,x,h,ytemp,yerr)
      errmax=0. !Evaluate accuracy.
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
 11   continue
      errmax=errmax/eps !Scale relative to required tolerance.
      if(errmax.gt.1.)then !Truncation error too large, reduce stepsize.
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(abs(htemp),0.1*abs(h)),h)!No more than a factor of 10
        xnew=x+h
        if(xnew.eq.x)pause 'stepsize underflow in rkqs'
        goto 1! For another try.
      else !Step succeeded. Compute size of next step.
        if(errmax.gt.ERRCON)then
            hnext=SAFETY*h*(errmax**PGROW)
        else !No more than a factor of 5 increase.
            hnext=5.*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
            y(i)=ytemp(i)
 12     continue
        return
      endif
      END 
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n,NMAX
      REAL*8 h,x,dydx(n),y(n),yerr(n),yout(n)
C     EXTERNAL derivs
      PARAMETER (NMAX=2) !Set to the maximum number of functions.
C     USES derivs
C     Given values for n variables y and their derivatives dydx known at
C     x, use the fifth-order Cash-Karp Runge-Kutta method to advance the
C     solution over an interval h and return the incremented variables 
C     as yout. Also return an estimate of the local truncation error in 
C     yout using the embedded fourth-order method. The user supplies the 
C     subroutine derivs(x,y,dydx), which returns derivatives dydx at x.
      INTEGER i
      REAL*8 ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),
     *  ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,
     *  B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,
     *  DC4,DC5,DC6
      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,
     *  B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,
     *  B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,
     *  B63=575./13824.,B64=44275./110592.,B65=253./4096.,
     *  C1=37./378.,C3=250./621.,C4=125./594.,C6=512./1771.,
     *  DC1=C1-2825./27648.,DC3=C3-18575./48384.,
     *  DC4=C4-13525./55296.,DC5=-277./14336.,DC6=C6-.25)
      do 11 i=1,n !First step.
        ytemp(i)=y(i)+B21*h*dydx(i)
 11   continue
      call derivs(x+A2*h,ytemp,ak2)! Second step.
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
 12   continue
      call derivs(x+A3*h,ytemp,ak3) !Third step.
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
 13   continue
      call derivs(x+A4*h,ytemp,ak4) !Fourth step.
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+
     *      B54*ak4(i))
 14   continue
      call derivs(x+A5*h,ytemp,ak5) !Fifth step.
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+
     *      B64*ak4(i)+B65*ak5(i))
 15   continue
      call derivs(x+A6*h,ytemp,ak6) !Sixth step.
      do 16 i=1,n ! Accumulate increments with proper weights.
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+
     *      C6*ak6(i))
 16   continue
      do 17 i=1,n
C     Estimate error as difference between fourth and fifth order 
C     methods.
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)
     *      +DC6*ak6(i))
 17   continue
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine derivs(x,y,dydx)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 x,y(2),dydx(2),onedzc
      common /fitting/ onedzc
      
c      write(0,*) "y",y(1),y(2)
      dydx(2)=y(1)
      dydx(1)=-2.0d0/x*y(1)-(y(2)*y(2)-onedzc*onedzc)**(3.0d0/2.0d0)
c      write(6,*) "dydx",dydx(1),dydx(2)
 
      return
      end