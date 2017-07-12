CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function Zflux(nrad,ntheta,nl,Pi)
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

      ZFlux=0.0d0 !Initialize sum for Total Flux.
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

            if((nl(3).eq.0.0).and.(nl(4).eq.0.0))then
               limb=1.0d0-nl(1)*(1.0d0-u)-nl(2)*(1.0d0-u)*(1.0d0-u)
            else
               limb=0.0d0
               do 14 ii=1,4
                  limb=limb+nl(ii)*(1.0d0-u**dble(ii/2.0d0))
 14            continue
               limb=1.0d0-limb !set limb=1.0 to disable limb-darkening
            endif

            ZFlux=ZFlux+zonearea*In*limb

 11     continue
 10   continue

c      write(6,*) ZFlux

      return
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
            if((nl(3).eq.0.0).and.(nl(4).eq.0.0))then
               limb=1.0d0-nl(1)*(1.0d0-u)-nl(2)*(1.0d0-u)*(1.0d0-u)
            else
               limb=0.0d0
               do 14 ii=1,4
                  limb=limb+nl(ii)*(1.0d0-u**dble(ii/2.0d0))
 14            continue
               limb=1.0d0-limb !set limb=1.0 to disable limb-darkening
            endif

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
            if(yplus.ge.Rl*Rsun)then
               Ap=Aplus(dist,Re)
            else
               Ap=0.0d0
            endif
C           Now we check to see if - image is occulted
            ymin=0.5*(dist-sqrt(dist*dist+4.0d0*Re*Re))
            if(ymin.ge.Rl*Rsun)then
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

      ydRe=y/Re
      aplus=ydRe*ydRe+2.0d0
      aplus=aplus/(2.0d0*ydRe*sqrt(ydRe*ydRe+4.0d0))+0.5

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function amin(y,Re)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This is equation 13
      implicit none
      double precision y,Re,ydRe

      ydRe=y/Re
      amin=ydRe*ydRe+2.0d0
      amin=amin/(2.0d0*ydRe*sqrt(ydRe*ydRe+4.0d0))-0.5

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
