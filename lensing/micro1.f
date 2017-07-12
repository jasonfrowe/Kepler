      program microlensing
C     This program calculating microlensing for sphere of mass
C     transiting a star
      implicit none
      integer nmax,i,nrad,ntheta
      parameter(nmax=8000)
c      real xp(nmax),yp(nmax)
      double precision G,c,Msun,Rsun,Pi,AU,asemi,Rstar,Ml,Rl,phi0,phi,
     .   b,incl,day2sec,per,Mstar,dt,time,Psec,x,y,z,nl(4),
     .   Flux,F0,TFlux,PFlux,thres

C     Constants
      G=6.674d-11 !N m^2 kg^-2  Gravitation constamt
      c=2.99792458d8 !m/s Speed of light
      Msun=1.9891d30 !kg  mass of Sun
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      Pi=acos(-1.d0)   !Pi
      AU=1.4959787069d11 !m  astronomical unit
      day2sec=24.0d0*60.0d0*60.0d0 !sec/Day
      nrad=15706   !number of rings to integrate star/planet surface
      ntheta=15706 !maximum number of angular divisions
C     thres is used to skip part of the phase when lensing is not
C     important (saves a lot of time!)
      thres=1.0e-7 !tells us when amplification is negliable

C     We are going to work with a circular orbit as a test.
C     phi=true anomaly [-Pi:Pi] - transit when phi=0.
C     x = distance to left and right of source (AU)
C     y = distance up or down relative to source (AU)
C     z = distance to or away from the observer (AU)

C     Orbital Parameters (e.g. model parameters)
      Per=288.0!365.25d0 !orbital period (days)
      Mstar=1.0d0 !mass of star (Msun)
      Rstar=1.0d0 !radius of star (Rsun)
      Ml=0.6d0 !mass of lens (Msun)
      Rl=0.01d0 !radius of lens (Rsun)
      b=0.5d0 !minimum impact parameter
      phi0=-pi/2 !phase offset (when the transit occurs)

c      Per=23.8d0 !orbital period (days)
c      Mstar=2.7d0 !mass of star (Msun)
c      Rstar=2.9d0 !radius of star (Rsun)
c      Ml=0.2d0 !mass of lens (Msun)
c      Rl=0.1d0 !radius of lens (Rsun)
c      b=0.5d0 !minimum impact parameter
c      phi0=-pi/2 !phase offset (when the transit occurs)
c
c      Per=8.6d0 !orbital period (days)
c      Mstar=0.8d0 !mass of star (Msun)
c      Rstar=0.00955d0 !radius of star (Rsun)
c      Ml=1.97d0 !mass of lens (Msun)
c      Rl=0.0000143781 !radius of lens (Rsun)
c      b=37.5d0 !minimum impact parameter
c      phi0=-pi/2 !phase offset (when the transit occurs)
c
c      Per=0.625d0 !orbital period (days)
c      Mstar=0.460d0 !mass of star (Msun)
c      Rstar=1.00d0 !radius of star (Rsun)
c      Ml=1.5d0 !mass of lens (Msun)
c      Rl=0.01!0014362d0 !radius of lens (Rsun)
c      b=0.542d0 !minimum impact parameter
c      phi0=-pi/2 !phase offset (when the transit occurs)

C     Limb-darkening co-efficients
      nl(1)= 0.5118  !these values are for the Sun with Kepler bandpass
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
      write(0,*) incl/(Pi/2.0d0)*90.0d0

C     Loop over 1 orbital period and sample nmax times.
      dt=Psec/dble(nmax) !time step (s)
      time=0.0d0 !initialize time (s)

C     Get unlensed/untransited flux for normalization
      F0=Tflux(nrad,ntheta,nl,Pi)
      write(0,*) "F0:",F0

c      time=7.86e6
c      dt=0.0001e6
c      do 10 i=1,600


      do 10 i=1,nmax
c       time=7.89d6
         phi=2.0*Pi*(time/Psec-int(time/Psec))+phi0 !true anomaly
         x=asemi*sin(phi)
         y=asemi*cos(phi)*cos(incl)
         z=asemi*cos(phi)

         if(z.gt.0.0)then !only lensing when transiting
            Flux=PFlux(nrad,ntheta,x,y,z,Rstar,Rl,Ml,nl,Pi,G,c,
     .         Msun,Rsun,AU,F0,thres)
            write(6,*) time/day2sec,Flux/F0
         else
C           Add occultation stuff here.
            Flux=1.0d0
         endif


         time=time+dt
 10   continue

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function PFlux(nrad,ntheta,x,y,z,Rstar,Rl,Ml,
     .   nl,Pi,G,c,Msun,Rsun,AU,F0,thres)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nrad,ntheta,ntheta2,i,j,ii
      double precision nl(4),Pi,Ml,x,y,z,Rstar,Rl,G,c,Msun,Rsun,AU,Re,
     .   einsrad,thres,dphi,dR,R,dtheta,dntheta,
     .   ringarea,aplus,ap,dist,amin,am,twoPi,twon,theta,zonearea,xp,yp,
     .   In,u,limb,Tr,Amp,yplus,ymin,F0

      twoPi=2.0d0*Pi !2Pi
      dR=1.0/dble(nrad)/2.0d0 !number of divisions in radius
      twon=dble(2*nrad) !preconvert to dble

      Re=einsrad(Ml,z,G,c,Msun,AU) !determine Einstein Radius (m)
!      write(0,*) "Re: ",Re

C     Start by inspecting closest element..
      dist=AU*sqrt(x*x+y*y)-Rstar*Rsun
      ap=aplus(dist,Re) !this is the largest A+ can be.
C     If the max-amplification is too small, then return unamped flux
C     and exit
      if((ap-1.0d0.lt.thres).and.(dist.gt.0.0))then
         Pflux=F0
         return
      endif


      PFlux=0.0d0 !Initialize sum for Total Flux.
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

            In=1.0d0 !intensity of an element

C     The next 6 lines compute limb-darkening.
            u=sqrt(1-(xp*xp+yp*yp))  !add effects of limb-darkening
            limb=0.0d0
            do 14 ii=1,4
                limb=limb+nl(ii)*(1.0d0-u**dble(ii/2.0d0))
 14         continue
            limb=1.0d0-limb !set limb=1.0 to disable limb-darkening

C           This line is the same as Equation 10 (y)
            dist=sqrt((x*AU-xp*rstar*Rsun)**2.0d0+
     .                (y*AU-yp*rstar*Rsun)**2.0d0)


C     These few lines are for transits, but I've disabled it for now
C     Not sure if this is right.
c            if(dist.lt.Rl*Rsun)then
c               Tr=0.0 !transiting
c            else
               Tr=1.0 !non-transiting
c            endif

C           Now we check to see if + image is occulted
            yplus=0.5*(dist+sqrt(dist*dist+4.0d0*Re*Re))
            if(yplus.gt.Rl*Rsun)then
               Ap=Aplus(dist,Re)
            else
               Ap=0.0d0
            endif
C           Now we check to see if - image is occulted
            ymin=0.5*(dist-sqrt(dist*dist+4.0d0*Re*Re))
            if(ymin.gt.Rl*Rsun)then
               Am=Amin(dist,Re)
            else
               Am=0.0d0
            endif
            Amp=Ap+Am  !total amplification

            PFlux=PFlux+zonearea*In*limb*Tr*Amp !sum up the flux

 11     continue
 10   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function aplus(y,Re)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This is equation 13
      implicit none
      double precision y,Re,ydRe

      if(Re.gt.0.0d0)then
         ydRe=y/Re
         aplus=ydRe*ydRe+2.0d0
         aplus=aplus/(2.0d0*ydRe*sqrt(ydRe*ydRe+4.0d0))+0.5
      else
         aplus=1.0d0
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function amin(y,Re)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This is equation 13
      implicit none
      double precision y,Re,ydRe

      if(Re.gt.0.0d0)then
         ydRe=y/Re
         amin=ydRe*ydRe+2.0d0
         amin=amin/(2.0d0*ydRe*sqrt(ydRe*ydRe+4.0d0))-0.5
      else
         amin=-1.0d0
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function Tflux(nrad,ntheta,nl,Pi)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Calculates flux from star with out transits or lensing.
C     Used to normalize the lightcurve
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
