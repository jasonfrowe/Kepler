CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program asemidradstar
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
C     Calculate the ratio of the semi-major axis to the stellar radius
      double precision Rs,Rp,ideg,R1,R2,incl
C     Physical Constants
      double precision Pi,G,Msun,Mearth,Mjup,Rsun,Rjup,a,b,c,x1,x2,
     .  rtsafe,adrs
      common /fsolve/ a,c
C     Useful constants
      Pi=acos(-1.d0)   !Pi
      G=6.674d-11 !N m^2 kg^-2  Gravitation constamt
      Msun=1.9891d30 !kg  mass of Sun
      Mearth=5.974d24 !kg mass of Earth
      Mjup=317.833d0*Mearth !kg  mass of Jupiter
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      Rjup=142980.0d0*1000.0d0/2.0d0 !m  radius of Jupiter
      
      Rs=1.84d0
      Rp=1.38d0
      ideg=84.65
      
      R1=Rs*Rsun
      R2=Rp*Rjup
      incl=Pi*ideg/180.0d0
      
      a=(1.0d0+R2/R1)*(1.0d0+R2/R1)
      c=(cos(incl))**2.0
      x1=0.0
      x2=10.0
c      adrs=rtsafe(x1,x2,1.0d-6)
      adrs=sqrt(a*(1-c)/(c*c))
      
      write(6,*) "ADRS: ",adrs

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine funcd(b,y,dy)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision a,b,c,y,dy
      common /fsolve/ a,c
      
      y=b**3*c**2-b**2*c*(1.0d0-c**2)-b**2*c**2-a*b+a*c
      dy=3.0d0*b**2*c**2-2.0d0*b*c**2-2.0d0*b*(1.0d0-c**2)*c-a

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION rtsafe(x1,x2,xacc)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER MAXIT
      REAL*8 rtsafe,x1,x2,xacc
      PARAMETER (MAXIT=100) !Maximum allowed number of iterations.
      INTEGER j
      REAL*8 df,dx,dxold,f,fh,fl,temp,xh,xl
      call funcd(x1,fl,df)
      call funcd(x2,fh,df)
      if((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.))
     .  pause 'root must be bracketed in rtsafe'
      if(fl.eq.0.)then
        rtsafe=x1
        return
      elseif(fh.eq.0.)then
        rtsafe=x2
        return
      elseif(fl.lt.0.)then !Orient the search so that f(xl) < 0.
        xl=x1
        xh=x2
      else
        xh=x1
        xl=x2
      endif
      rtsafe=.5*(x1+x2) !Initialize the guess for root,
      dxold=abs(x2-x1) !the “stepsize before last,”
      dx=dxold !and the last step.
      call funcd(rtsafe,f,df)
      do 11 j=1,MAXIT !Loop over allowed iterations.
C       Bisect if Newton out of range, or not decreasing fast enough.
        if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).gt.0. 
     .    .or. abs(2.*f).gt.abs(dxold*df) ) then 
            dxold=dx
            dx=0.5*(xh-xl)
            rtsafe=xl+dx
            if(xl.eq.rtsafe)return !Change in root is negligible.
        else !Newton step acceptable. Take it.
            dxold=dx
            dx=f/df
            temp=rtsafe
            rtsafe=rtsafe-dx
            if(temp.eq.rtsafe)return
        endif
        if(abs(dx).lt.xacc) return !Convergence criterion.
C       The one new function evaluation per iteration.
        call funcd(rtsafe,f,df) 
        if(f.lt.0.) then !Maintain the bracket on the root.
            xl=rtsafe
        else
            xh=rtsafe
        endif
 11   continue
      pause 'rtsafe exceeding maximum iterations'
      return
      END
     