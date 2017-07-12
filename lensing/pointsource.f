      program pointsource
C     This program tests Eqn. 4 for Sahu & Gilliland 2003
      implicit none
      integer i,nmax,nc
      parameter(nmax=4000) !number of points in our curves.
      real xp(nmax),yp(nmax) !use real for ploting only
      double precision G,c,Msun,Rsun,einsrad,x,y,z,Pi,l,amin,amax,da,
     .   phi,dphi,AU,a,Ap,Re,Ml,D,u,Rs,phi0,incl,yoffset


C     Constants
      G=6.674d-11 !N m^2 kg^-2  Gravitation constamt
      c=2.99792458d8 !m/s Speed of light
      Msun=1.9891d30 !kg  mass of Sun
      Rsun=6.9599d8 !m radius of Sun
      Pi=acos(-1.d0)   !Pi
      AU=1.4959787069d11 !m  astronomical unit
      phi0=1.0d-4
      dphi=2.0d0*Phi0/dble(nmax) !delta change in phase

C     We are going to work with a circular orbit as a test.
C     phi=true anomaly [-Pi:Pi]
C     x = distance to left and right of source (AU)
C     y = distance up or down relative to source (AU)
C     z = distance to or away from the observer (AU)

C     Orbital Parameters
      amin = 0.1d0  !minimum semi-major axis to test (AU)
      amax = 1.0d0  !maximum semi-major axis to test (AU)
      da = 0.1  !step size in semi-major axis (AU)
      Ml=0.6*Msun !mass of lens
      yoffset=0.5e6 !how close to crossing the point source (m)

      a=amin !initial orbital separation (AU)

C     All the 'pg..' calls are for plotting purposes only.
      call pgopen('?')
      call PGPAP ( 6.0 ,1.0)
      call pgvport(0.15,0.85,0.2,0.9)
      call pgpage()
      call pgwindow(real(-phi0),real(phi0),0.9,50.0)
      call pgslw(3)
      call pgsch(1.5)
      CALL PGBOX('BCNTS2',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Phi (radians)","Amplification","")
      nc=2


      do 11 while(a.le.amax) !loop over all separations
C        inclination of system (Pi/2 is edge on)
         incl=Pi/2.0d0-tan(yoffset/(a*AU))

C        We will define the center of transits at phi=0.
         phi=-phi0 !initize phase
         do 10 i=1,nmax !lets loop over 1 phase and calculate Ap (Eqn 4)
            x=a*sin(phi)
            y=a*cos(phi)*cos(incl)
            z=a*cos(phi)

C           only calculate magnitifcation when lens is in front
            if(z.gt.0.0d0)then
C              Calculate Einstein Radius
               D=z*AU !projected lense-source distance (Fig 1)
               Re=einsrad(Ml,D,G,c) !Equation 1
               Rs=AU*sqrt(x*x+y*y) !projected distance of lens-source
               u=Rs/Re !see Equation 3
               Ap=(u*u+2.0d0)/(u*sqrt(u*u+4.0d0)) !Equation 4
            else
               Ap=0.0d0
            endif
            xp(i)=real(phi) !values for plotting
            yp(i)=real(Ap)

            phi=phi+dphi !update phi
 10      continue
         call pgsci(nc)
         call pgline(nmax,xp,yp)
         nc=nc+1

      a=a+da !next orbital separation
 11   continue


      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function einsrad(Ml,D,G,c)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This function calculates the Einstein Radius.
C     Assumes that D=Dls*Dl/Ds~Dls
      implicit none
      double precision Ml,D,G,c

      einsrad=4.0d0*G*Ml*D/(c*c) !Equation 1
      einsrad=sqrt(einsrad)

      return
      end



