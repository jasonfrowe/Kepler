CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER nv,NMAX,KMAXX,IMAX
      REAL*8 eps,hdid,hnext,htry,x,dydx(nv),y(nv),yscal(nv),SAFE1,
     * SAFE2,REDMAX,REDMIN,TINY,SCALMX
      PARAMETER (NMAX=6000,KMAXX=8,IMAX=KMAXX+1,SAFE1=.25,SAFE2=.7,
     * REDMAX=1.d-5,REDMIN=.7,TINY=1.d-30,SCALMX=.1)
C     USES derivs,mmid,pzextr
c     Bulirsch-Stoer step with monitoring of local truncation error to
C     ensure accuracy and adjust stepsize. Input are the dependent
C     variable vector y(1:nv) and its derivative dydx(1:nv) at the
C     starting value of the independent variable x. Also input are the
C     stepsize to be attempted htry, the required accuracy eps, and the
C     vector yscal(1:nv) against which the error is scaled. On output,
C     y and x are replaced by their new values, hdid is the stepsize
C     that was actually accomplished, and hnext is the estimated next
C     stepsize. derivs is the user-supplied subroutine that computes the
C     right-hand side derivatives. Be sure to set htry on successive
C     steps to the value of hnext returned from the previous step, as is
C     the case if the routine is called by odeint.
C
c     Parameters: NMAX is the maximum value of nv; KMAXX is the maximum
C     row number used in the extrapolation; IMAX is the next row number;
C     SAFE1 and SAFE2 are safety factors; REDMAX is the maximum factor
C     used when a stepsize is reduced, REDMIN the minimum; TINY prevents
C     division by zero; 1/SCALMX is the maximum factor by which a
C     stepsize can be increased.

      INTEGER i,iq,k,kk,km,kmax,kopt,nseq(IMAX)
      REAL*8 eps1,epsold,errmax,fact,h,red,scale,work,wrkmin,xest,
     * xnew,a(IMAX),alf(KMAXX,KMAXX),err(KMAXX),yerr(NMAX),
     * ysav(NMAX),yseq(NMAX)
      LOGICAL first,reduct
      SAVE a,alf,epsold,first,kmax,kopt,nseq,xnew
      EXTERNAL derivs
      DATA first/.true./,epsold/-1./
      DATA nseq /2,4,6,8,10,12,14,16,18/
      if(eps.ne.epsold)then !A new tolerance, so reinitialize.
         hnext=-1.d29 !“Impossible” values.
         xnew=-1.d29
         eps1=SAFE1*eps
         a(1)=nseq(1)+1 !Compute work coefficients Ak.
         do 11 k=1,KMAXX
            a(k+1)=a(k)+nseq(k+1)
 11      continue
         do 13 iq=2,KMAXX !Compute α(k, q).
            do 12 k=1,iq-1
               alf(k,iq)=eps1**((a(k+1)-a(iq+1))/
     *            ((a(iq+1)-a(1)+1.)*(2*k+1)))
 12         continue
 13      continue
         epsold=eps
         do 14 kopt=2,KMAXX-1 !Determine optimal row number for convergence.
            if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1 
 14      continue
  1      kmax=kopt
      endif
      h=htry
      do 15 i=1,nv !Save the starting values.
         ysav(i)=y(i)
 15   continue
      if(h.ne.hnext.or.x.ne.xnew)then !A new stepsize or a new integration: re-establish
         first=.true. !the order window.
         kopt=kmax
      endif
      reduct=.false.
 2    do 17 k=1,kmax !Evaluate the sequence of modified midpoint
         xnew=x+h !integrations.
         if(xnew.eq.x)then
            write(6,*) "step size underflow in bsstep",h,x
            read(5,*)
         endif
         call mmid(ysav,dydx,nv,x,h,nseq(k),yseq,derivs)
         xest=(h/nseq(k))**2 !Squared, since error series is even.
         call pzextr(k,xest,yseq,y,yerr,nv) !Perform extrapolation.
         if(k.ne.1)then !Compute normalized error estimate (k).
            errmax=TINY
            do 16 i=1,nv
               errmax=max(errmax,abs(yerr(i)/yscal(i)))
 16         continue
            errmax=errmax/eps !Scale error relative to tolerance.
            km=k-1
            err(km)=(errmax/SAFE1)**(1./(2*km+1))
         endif
         if(k.ne.1.and.(k.ge.kopt-1.or.first))then !In order window.
            if(errmax.lt.1.)goto 4 !Converged.
            if(k.eq.kmax.or.k.eq.kopt+1)then !Check for possible stepsize reduction.
               red=SAFE2/err(km)
               goto 3
            else if(k.eq.kopt)then
               if(alf(kopt-1,kopt).lt.err(km))then
                  red=1./err(km)
                  goto 3
               endif
            else if(kopt.eq.kmax)then
               if(alf(km,kmax-1).lt.err(km))then
                  red=alf(km,kmax-1)*
     *               SAFE2/err(km)
                  goto 3
               endif
            else if(alf(km,kopt).lt.err(km))then
               red=alf(km,kopt-1)/err(km)
               goto 3
            endif
         endif
 17   continue
 3    red=min(red,REDMIN) !Reduce stepsize by at least REDMIN and at
      red=max(red,REDMAX) !most REDMAX.
      h=h*red
      reduct=.true.
      goto 2 !Try again.
 4    x=xnew !Successful step taken.
      hdid=h
      first=.false.
      wrkmin=1.d35 !Compute optimal row for convergence and
      do 18 kk=1,km !corresponding stepsize.
         fact=max(err(kk),SCALMX)
         work=fact*a(kk+1)
         if(work.lt.wrkmin)then
            scale=fact
            wrkmin=work
            kopt=kk+1
         endif
 18   continue
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then !Check for possible order increase,
C     but not if stepsize was just reduced.
         fact=max(scale/alf(kopt-1,kopt),SCALMX)
         if(a(kopt+1)*fact.le.wrkmin)then
            hnext=h/fact
            kopt=kopt+1
         endif
      endif
      return
      END
     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      SUBROUTINE pzextr(iest,xest,yest,yz,dy,nv)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER iest,nv,IMAX,NMAX
      REAL*8 xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13,NMAX=6000)
C     Use polynomial extrapolation to evaluate nv functions at x = 0 by
C     fitting a polynomial to a sequence of estimates with progressively
C     smaller values x = xest, and corresponding function vectors
C     yest(1:nv). This call is number iest in the sequence of calls.
C     Extrapolated function values are output as yz(1:nv), and their
C     estimated error is output as dy(1:nv). Parameters: Maximum
C     expected value of iest is IMAX; of nv is NMAX.
      INTEGER j,k1
      REAL*8 delta,f1,f2,q,d(NMAX),qcol(NMAX,IMAX),x(IMAX)
      SAVE qcol,x
      x(iest)=xest !Save current independent variable.
      do 11 j=1,nv
         dy(j)=yest(j)
         yz(j)=yest(j)
 11   continue
      if(iest.eq.1) then !Store first estimate in first column.
         do 12 j=1,nv
            qcol(j,1)=yest(j)
 12      continue
      else
         do 13 j=1,nv
            d(j)=yest(j)
 13      continue
         do 15 k1=1,iest-1
            delta=1./(x(iest-k1)-xest)
            f1=xest*delta
            f2=x(iest-k1)*delta
            do 14 j=1,nv !Propagate tableau 1 diagonal more.
               q=qcol(j,k1)
               qcol(j,k1)=dy(j)
               delta=d(j)-q
               dy(j)=f1*delta
               d(j)=f2*delta
               yz(j)=yz(j)+dy(j)
 14         continue
 15      continue
         do 16 j=1,nv
            qcol(j,iest)=dy(j)
 16      continue
      endif
      return
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER nstep,nvar,NMAX
      REAL*8 htot,xs,dydx(nvar),y(nvar),yout(nvar)
      EXTERNAL derivs
      PARAMETER (NMAX=6000)
C     Modified midpoint step. Dependent variable vector y(1:nvar) and
C     its derivative vector dydx(1:nvar) are input at xs. Also input is
C     htot, the total step to be made, and nstep, the number of substeps
C     to be used. The output is returned as yout(1:nvar), which need
C     not be a distinct array from y; if it is distinct, however, then
C     y and dydx are returned undamaged.
      INTEGER i,n
      REAL*8 h,h2,swap,x,ym(NMAX),yn(NMAX)
      h=htot/nstep !Stepsize this trip.
      do 11 i=1,nvar
         ym(i)=y(i)
         yn(i)=y(i)+h*dydx(i) !First step.
 11   continue
      x=xs+h
C     Will use yout for temporary storage of derivatives.
      call derivs(x,yn,yout)
      h2=2.*h
      do 13 n=2,nstep !General step.
         do 12 i=1,nvar
            swap=ym(i)+h2*yout(i)
            ym(i)=yn(i)
            yn(i)=swap
 12      continue
         x=x+h
         call derivs(x,yn,yout)
 13   continue
      do 14 i=1,nvar !Last step.
         yout(i)=0.5*(ym(i)+yn(i)+h*yout(i))
 14   continue
      return
      END

