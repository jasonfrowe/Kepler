      program microlensing
C     This program calculating microlensing for sphere of mass
C     transiting a star
      implicit none
      integer nmax,i,nrad,ntheta
      parameter(nmax=4000)
c      real xp(nmax),yp(nmax)
      double precision G,c,Msun,Rsun,Pi,AU,asemi,Rstar,Ml,Rl,phi0,phi,
     .   b,incl,day2sec,per,Mstar,dt,time,Psec,x,y,z,nl(4),lensing,
     .   Flux,F0,TFlux,PFlux

C     Constants
      G=6.674d-11 !N m^2 kg^-2  Gravitation constamt
      c=2.99792458d8 !m/s Speed of light
      Msun=1.9891d30 !kg  mass of Sun
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      Pi=acos(-1.d0)   !Pi
      AU=1.4959787069d11 !m  astronomical unit
      day2sec=24.0d0*60.0d0*60.0d0 !sec/Day
      nrad=100   !number of rings to integrate star/planet surface
      ntheta=100 !maximum number of angular divisions

C     We are going to work with a circular orbit as a test.
C     phi=true anomaly [-Pi:Pi] - transit when phi=0.
C     x = distance to left and right of source (AU)
C     y = distance up or down relative to source (AU)
C     z = distance to or away from the observer (AU)

C     Orbital Parameters (e.g. model parameters)
      Per=365.25d0 !orbital period (days)
      Mstar=1.0d0 !mass of star (Msun)
      Rstar=1.0d0 !radius of star (Rsun)
      Ml=0.6d0 !mass of lens (Msun)
      Rl=0.01d0 !radius of lens (Rsun)
      b=0.5d0 !minimum impact parameter
      phi0=-pi/2 !phase offset (when the transit occurs)

C     Limb-darkening co-efficients
      nl(1)= 0.5118
      nl(2)= 0.0525
      nl(3)= 0.4590
      nl(4)=-0.2727

C     derive semi-major axis (AU)
      Psec=Per*day2sec !period in seconds
      asemi=Psec*Psec*G*(Mstar+Ml)*Msun/(4.0d0*pi*pi)
      asemi=asemi**(1.0d0/3.0d0)/AU
c      write(0,*) asemi

C     orbital inclination - central transit is incl=pi/2
      incl=Pi/2.0d0-tan(b*Rstar*Rsun/(asemi*AU)) !(radians)
c      write(6,*) incl/(Pi/2.0d0)*90.0d0

C     Loop over 1 orbital period and sample nmax times.
      dt=Psec/dble(nmax) !time step (s)
      time=0.0d0 !initialize time (s)

C     Get unlensed/untransited flux for normalization
      F0=Tflux(nrad,ntheta,nl,Pi)
c      write(6,*) "F0:",F0

      time=7.86e6
      dt=0.0001e6
      do 10 i=1,600
c
c       time=7.89d6
         phi=2.0*Pi*(time/Psec-int(time/Psec))+phi0 !true anomaly
         x=asemi*sin(phi)
         y=asemi*cos(phi)*cos(incl)
         z=asemi*cos(phi)

         if(z.gt.0.0)then !only lensing when transiting
            Flux=PFlux(nrad,ntheta,x,y,z,Rstar,Rl,Ml,asemi,nl,Pi,G,c,
     .         Msun,Rsun,AU)
            write(6,*) time/day2sec,Flux/F0
         endif


         time=time+dt
 10   continue

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function PFlux(nrad,ntheta,x,y,z,Rstar,Rl,Ml,
     .   asemi,nl,Pi,G,c,Msun,Rsun,AU)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nrad,ntheta,ntheta2
      double precision nl(4),Pi,Ml,x,y,z,Rstar,Rl,G,c,Msun,Rsun,AU,Re,
     .   einsrad,thres,MaxAmpm1,asemi,dphi,dR,R,dtheta,dntheta,
     .   ringarea,aplus,ap,dist


      thres=1.0e-7 !tells us when amplification is negliable

      Re=einsrad(Ml,z,G,c,Msun,AU) !determine Einstein Radius (m)

C     Start my inspecting closest element..
      dist=AU*sqrt(x*x+y*y)-Rstar*Rsun
      ap=aplus(dist,Re)
      write(0,*) "ap:",ap

      dphi=1.0e-9 !smallest circle to consider for integration (radian)
      MaxAmpm1=0.1d0 !we check max Amplitude minus 1 vs threshold

      R=0.0 !initial integral radius (Rsun)
      dR=AU*asemi*sin(dphi)/Rsun !delta radius (Rsun)
      Pflux=0.0d0
c      do 10 while(MaxAmpm1.gt.thres)
         R=R+dR !increment integral radius
C        Calculate area (Rsun^2)
         ringarea=pi*((R+dR)*(R+dR)-(R-dR)*(R-dR)) !piR2^2 - piR1^2
C        split up ring
         ntheta2=ntheta!max(int(dble(ntheta)*R),10)
         dtheta=2.0d0*Pi/dble(ntheta2) !(radians)
         dntheta=dble(ntheta2) !precompute number of rings

cC        Now we sum up each segment of the ring
c         do 11 j=1,ntheta2
c            theta=twopi*dble(j)/dntheta !theta (radians)
c            zonearea=ringarea*dtheta/twopi !dtheta/2pi (Rsun^2)
c
cC           calculate the (x,y) position of the surface element relative
cC           to the center of the source projected on the lens plane
c            xp=Rsun*R*cos(theta)+x  !(m)
c            yp=Rsun*R*sin(theta)+y  !(m)
c
c 11      continue

c 10   continue

c
c
c
cC           This line is the same as Equation 10 (y)
c            dist=sqrt((x*AU-xp*rstar*Rsun)**2.0d0+
c     .                (y*AU-yp*rstar*Rsun)**2.0d0)
c
c            if(dist.lt.Rl*Rsun)then
c               Tr=0.0 !transiting
c            else
c               Tr=1.0 !non-transiting
c            endif
c
c            Amp=0.0 !initialize Amplificiation to 1.0
c
cC           Now we check to see if + image is occulted
c            yplus=0.5*(dist+sqrt(dist*dist+4.0d0*Re*Re))
c            ydre=dist/Re
cc            write(6,*) "yplus",yplus,rl*Rsun
c            if(yplus.ge.Rl*Rsun)then
c              Aplus=(ydre*ydre+2.0d0)/(2*ydre*sqrt(ydre*ydre+4.0d0))+0.5
c            else
c               Aplus=0.0d0
c            endif
cC           Now we check to see if + image is occulted
c            ymin=0.5*(dist-sqrt(dist*dist+4.0d0*Re*Re))
c            if(ymin.ge.Rl*Rsun)then
c              Amin=(ydre*ydre+2.0d0)/(2*ydre*sqrt(ydre*ydre+4.0d0))-0.5
c            else
c               Amin=0.0d0
c            endif
c            Amp=Aplus+Amin
cc            write(6,*) Amp,Aplus,Amin
cc            if(Amp.lt.1.0d0) Amp=1.0d0
c
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function aplus(y,Re)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision y,Re,ydRe

      ydRe=y/Re
      aplus=ydRe*ydRe+2.0d0
      aplus=aplus/(2.0d0*ydRe*sqrt(ydRe*ydRe+4.0d0))+0.5

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function Tflux(nrad,ntheta,nl,Pi)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nrad,ntheta,i,j,ntheta2,ii
      double precision x,y,z,Rstar,Rl,Pi,twoPi,dR,twon,nl(4),R,
     .   dtheta,dntheta,zonearea,ringarea,xp,yp,theta,In,u,limb

      twoPi=2.0d0*Pi !2Pi
      dR=1.0/dble(nrad)/2.0d0 !number of divisions in radius
      twon=dble(2*nrad) !preconvert to dble

c      write(6,*) "Re:",Re

      TFlux=0.0d0 !Initialize sum for Total Flux.
      do 10 i=1,nrad
        R=1.0d0*dble(2*i-1)/twon !R
        ringarea=pi*((R+dR)*(R+dR)-(R-dR)*(R-dR)) !piR2^2 - piR1^2
        ntheta2=max(int(dble(ntheta)*R),10)
        dtheta=2.0d0*Pi/dble(ntheta2)
        dntheta=dble(ntheta2)
        do 11 j=1,ntheta2
            theta=twopi*dble(j)/dntheta !theta
            zonearea=ringarea*dtheta/twopi !dtheta/2pi

C           calculate the x,y position of the surface element
            xp=R*cos(theta)
            yp=R*sin(theta)

            In=1.0d0 !intensity of element

            u=sqrt(1-(xp*xp+yp*yp))  !add effects of limb-darkening
            limb=0.0d0
            do 14 ii=1,4
                limb=limb+nl(ii)*(1.0d0-u**dble(ii/2.0d0))
 14         continue
            limb=1.0d0-limb !set limb=1.0 to disable limb-darkening

            TFlux=TFlux+zonearea*In*limb

 11     continue
 10   continue

c      write(6,*) TFlux

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function einsrad(Ml,D,G,c,Msun,AU)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This function calculates the Einstein Radius.
C     Assumes that D=Dls*Dl/Ds~Dls
      implicit none
      double precision Ml,D,G,c,Msun,AU

      einsrad=4.0d0*G*Ml*Msun*D*AU/(c*c) !Equation 1
      einsrad=sqrt(einsrad)

      return
      end
