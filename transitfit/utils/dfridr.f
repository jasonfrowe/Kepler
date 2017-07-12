CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION dfridr(x,h,err,n,sol,na,time,itime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER NTAB
      REAL*8 dfridr,err,h,x,dfuncs,CON,CON2,BIG,SAFE
      PARAMETER (CON=1.4D0,CON2=CON*CON,BIG=1.D30,NTAB=10,SAFE=2.0D0)
      integer n,na
      real*8 sol(na),time,itime
      EXTERNAL dfunc
C     Returns the derivative of a function func at a point x by Ridders' 
C     method of polynomial extrapolation. The value h is input as an 
C     estimated initial stepsize; it need not be small, but rather 
C     should be an increment in x over which func changes substantially. 
C     An estimate of the error in the derivative is returned as err.
C     Parameters: Stepsize is decreased by CON at each iteration. 
C     Max size of tableau is set by NTAB. 
C     Return when error is SAFE worse than the best so far.  
      INTEGER i,j
      REAL*8 errt,fac,hh,a(NTAB,NTAB)
      if(h.eq.0.) pause 'h must be nonzero in dfridr'
      hh=h
c      a(1,1)=(func(x+hh)-func(x-hh))/(2.0*hh)
      a(1,1)=(dfuncs(x+hh,n,sol,na,time,itime)-
     .  dfuncs(x-hh,n,sol,na,time,itime))/(2.0*hh)
c      write(0,*) "a(1,1):",a(1,1)
      err=BIG
C     Successive columns in the Neville tableau will go to smaller 
C     stepsizes and higher orders of extrapolation.

      do 12 i=2,NTAB 
        hh=hh/CON 
C       Try new, smaller stepsize.
c        a(1,i)=(func(x+hh)-func(x-hh))/(2.0*hh)
         a(1,i)=(dfuncs(x+hh,n,sol,na,time,itime)-
     .       dfuncs(x-hh,n,sol,na,time,itime))/(2.0*hh)
        fac=CON2
C       Compute extrapolations of various orders, requiring no new
C       function evaluations.
        do 11 j=2,i 
            a(j,i)=(a(j-1,i)*fac-a(j-1,i-1))/(fac-1.) 
            fac=CON2*fac
            errt=max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
C           The error strategy is to compare each new extrapolation to 
C           one order lower, both at the present stepsize and the 
C           previous one.

C           If error is decreased, save the improved answer.
            if (errt.le.err) then 
                err=errt
                dfridr=a(j,i)
            endif
 11     continue
        if(abs(a(i,i)-a(i-1,i-1)).ge.SAFE*err)return
C       If higher order is worse by a significant factor SAFE, then quit 
C       early.
 12   continue
      return
      END   
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function dfuncs(x,n,sol,na,time,itime) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER na,n
      REAL*8 x,sol(na),y
      integer nfit,i
      parameter(nfit=16)
      real*8 sol2(nfit),time,itime
      
C     x is the parameter to get the derivative for (pars(n))
C     a(na) are all parameters for the solution
C     time,itime are needed by the model

      do 10 i=1,na
        sol2(i)=sol(i)
 10   continue
      sol2(n)=x
      
C     model is returned in 'y'
      call transitmodel(1,time,itime,y,na,sol2)
 
      dfuncs=y !return this value

      return 
      END   