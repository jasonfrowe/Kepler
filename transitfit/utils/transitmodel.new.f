CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine transitmodel(npt,time,exptime,dtype,tmodel,nfit,sol)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Updated phi -> epoch for stability 
C
C     The transitmodel for double precision
C      -takes into account:
C       1. non-linear limb-darkening
C       2. reflected light from planet
C       3. circular orbits only
      implicit none
C     npt: the number of points for the model to return
C     nfit: the number of parameters that defines the model
C     nintg: used to find the average flux over the integration time
      integer npt,nfit,nintg,i,j,dtype(npt)
      parameter(nintg=11)
C     time: time (middle) for each model point (days)
C     exptime: the integration time for each model point (days)
C     tmodel: the model for each time point
C     sol: the input parameters for the model
      double precision time(npt),exptime(npt),tmodel(npt),sol(nfit)
C     Physical parameters of the model
      double precision M1,M2,Psec,asemi,R1,R2,Ag,c(4),zpt,norm,
     .  incl,Per,phi(nintg),eccn,w,ted,inclmin,epoch
C     Parameters that work towards MandelAgol
      double precision x1,y1,t(nintg),dnintg,dnintgm1,x2(nintg),
     .  y2(nintg),phase,tpi,angle,darea,b0(nintg),mu(nintg),mandelagol,
     .  mulimb0(nintg),mulimbf(nintg,5),dist(nintg),et(2),etest,phi0
C     Albedo routines
      double precision albedomod,zarea,xintold,y2pold,eclmod,y2pold2,
     .  xintold2,eclmod2,test,tides,tflux
C     RV parameters
C     K - amplitude of RV
C     voff - radial velocity offset (m/s)
C     vmodel - model velocities (m/s)
      integer nptv
      double precision K,voff
C     Non-circular orbit variables and routines
C     Eanom - Eccentric anomaly
C     Manom - Mean anomaly
C     kepler - solves Kepler equations
C     trueanomaly - calcaulates the true anomaly
C     Tanom - True anomaly 
C     arad - distance between star and planet
      double precision Eanom,Manom,trueanomaly,Tanom(nintg),arad(nintg),
     .  distance,Pid2
C     adhoc parameter to account for dilution.
      double precision dilute
C     Read in physical constants (Pi,G,Msun,..etc.)
      include "utils/physcons.f"

C     Adhoc dilute parameter
      dilute=sol(18)

C     Parameters for non-circular orbits
      tpi=2.0d0*Pi
      Pid2=Pi/2.0d0
C     Added eccentricity and w
c      eccn=sol(14)
c      w=Pi*sol(15)/180.0d0 !convert to radians
      eccn=sqrt(sol(14)*sol(14)+sol(15)*sol(15)) !eccentricity
      if(eccn.ge.1.0) eccn=0.99
      if(eccn.eq.0.0d0)then
        w=0.0d0
      else
        w=atan(sol(15)/sol(14))
      endif
      if((sol(14).gt.0.0d0).and.(sol(15).lt.0.0d0)) w=tPi+w
      if((sol(14).lt.0.0d0).and.(sol(15).gt.0.0d0)) w=Pi+w
      if((sol(14).lt.0.0d0).and.(sol(15).lt.0.0d0)) w=Pi+w
      Eanom=w !starting guess for Eanom
      
C     Parameters for Thermal eclipse (in micromag)
      ted=sol(16)*1.0d-6
      
C     We start by importing the model parameters.
C     Changing units and assigning more useless variable names         
      M1=sol(1)*Msun !kg ; mass of star
      M2=sol(2)*Mjup !kg ; mass of planet
      Psec=sol(5)*8.64d4 !sec ; period of planet
      Per=dble(sol(5)) !days  period of planet orbit in days
C     calculated semi-major axis
      asemi=(Psec*Psec*G*(M1+M2)/(4.0*Pi*Pi))**(1.0d0/3.0d0) !m
      R1=abs(sol(3)*Rsun)  !radius of star
      R2=abs(sol(4)*Rjup)  !radius of planet
C     non-linear limb darkening 
      c(1)=sol(10)
      c(2)=sol(11)
      c(3)=sol(12)
      c(4)=sol(13)   
      zpt=sol(8) !zero point
      
C     inclination angle
      inclmin=180.0*tan((R1+R2)/asemi)/pi
      inclmin=90.0-inclmin
      incl=sol(6)
      if((inclmin.ge.0).and.(inclmin.le.90.0))then
c        write(6,*) incl,inclmin
        if(incl.lt.inclmin) incl=inclmin
        if(incl.gt.90.0d0)incl=180.0-incl    
        if(incl.lt.inclmin) incl=inclmin !check a second time
      endif
      incl=Pi*(90.0d0-incl)/180.0d0 !radians
      
C     observer-star-planet angle
C     Find phase at centre of transit
      epoch=sol(7)
      phi0=tpi*(epoch/per-floor(epoch/per))
      call kepler(phi0,Eanom,eccn)
      Tanom(1)=trueanomaly(eccn,Eanom)
      phi0=Tanom(1)+w

C     RV initialization
      K=2.0*pi*G*M2**3*(sin(incl+Pid2))**3/
     .  (Psec*(1.0d0-eccn*eccn)**(3.0d0/2.0d0)*(M1+M2)*(M1+M2))
      K=K**(1.0d0/3.0d0)
      voff=sol(17) !velocity offset

c     normalization constant for light curve
      norm=Pi
      
      x1=0. !initialize center of star to 0,0 co-ordinates
      y1=0.
      dnintg=dble(nintg) !convert integer to double
      dnintgm1=2.0*dnintg-2.0
C     These next two variables are for eclmod
      xintold=0. !integration width initialization
      y2pold=0.0d0 !initialization of projected planet-star dist. 
c      write(6,*) "model start" 
      do 10 i=1,npt  !loop to get models points of each time
        do 11 j=1,nintg  !get array of times spanning exposure time
C           These times are centered on time(i)
            t(j)=time(i)+exptime(i)*(2.0*dble(j)-dnintg-1.0)/dnintgm1
            phi(j)=t(j)/per-floor(t(j)/per)
            phi(j)=phi(j)*tPi
            call kepler(phi(j),Eanom,eccn)
            Tanom(j)=trueanomaly(eccn,Eanom) !get true anonaly
            phi(j)=Tanom(j)+w-phi0  !compute phase relative to observer
            if(phi(j).gt.Pi) phi(j)=phi(j)-tPi            
            arad(j)=distance(asemi,eccn,Tanom(j))
            x2(j)=arad(j)*Sin(phi(j))
            y2(j)=arad(j)*Sin(incl)*Cos(phi(j))
cw            x2(j)=arad(j)*Sin(Tanom(j))
cw            y2(j)=arad(j)*Sin(incl)*Cos(Tanom(j))
c            write(0,*) t(j),Tanom,Manom
c            read(5,*)
c         write(6,*) t(j),x2(j)/R1
 500     format(28(F7.4,1X))
 11     continue
c        write(6,*) t(1),x2(1),y2(1)
c        read(5,*)
 
c        Manom=tpi*time(i)/Per+phi !mean anomaly
c        if(Manom.ge.tpi) Manom=Manom-tpi*int(Manom/tpi)
c        call kepler(Manom,Eanom,eccn)
c        Tanom=trueanomaly(eccn,Eanom)
 
        if(dtype(i).eq.0)then !case for photometric data
            zarea=0.0!initialize zarea
            do 12 j=1,nintg
C               ellipsidal component
                tflux=tides(M1,M2,R1,asemi,incl,tanom(j),eccn,Pi)
                zarea=zarea+tflux
C               flux from star + planet
                Ag=sol(9)*R1*R1/(arad(j)*arad(j))  !albedo of planet
                zarea=zarea+albedomod(Pi,t(j),Per,Ag,R1,R2,Tanom(j)+w) 
cw              zarea=zarea+albedomod(Pi,t(j),Per,Ag,R1,R2,Tanom(j))
 12         continue
            zarea=zarea/dnintg

C            determine what part 
c            phase=tPi*time(i)/Per+phi
            phase=Tanom(nintg/2+1)+w
cw           phase=Tanom(nintg/2+1)
c            write(6,*) time(i),Tanom(nintg/2+1)
            phase=phase/tPi-int(phase/tPi)+0.5
            if(phase.lt.0.0) phase=phase+1.0
            phase=phase-floor(phase)

C           Calculating the projected planet distance
c            do 13 j=1,nintg
c               angle=tPi*t(j)/Per+phi
c               x2(j)=asemi*Sin(angle)
c               y2(j)=asemi*Sin(incl)*Cos(angle)
c 13         continue
 
C       0.25 < phase < 0.75, then planet is in front
c            write(6,*) "hello3",i
            if((phase.gt.0.25).and.(phase.le.0.75))then
C     Now we calculate the flux change from planet transiting the star 
                darea=(Pi*mandelagol(nintg,R1,R2,x1,x2,y1,y2,c,b0,mu,
     .              mulimb0,mulimbf,dist)+zarea)/norm
C      otherwise the star is in front
            else
C     Now we calculate the flux change from the star blocking the planet
                darea=0.0
                do 14 j=1,nintg
c                   if( (t(j).gt.80.5).and.(t(j).lt.81.5) )
C                   My good version is not working.. ugh
c                   darea=darea+eclmod(R1,R2,x1,x2(j),y1,y2(j),zarea,
C    .                   norm,Pi,xintold,y2pold,ted)
C                   use approximation for now...
                    darea=darea+eclmod2(R1,R2,x1,x2(j),y1,y2(j),zarea,
     .                  norm,Pi,ted)
c                   darea=darea+(Pi+zarea)/norm
                    if(j.eq.1)then
                        xintold2=xintold
                        y2pold2=y2pold
                    endif
 14             continue
                darea=darea/dnintg
                xintold=xintold2
                y2pold=y2pold2
            endif
c           write(6,*) "hello4",i
C     Convert relative fluxes to magnitude to match observations
            tmodel(i)=zpt-2.5*log10(darea+(1.0-darea)*dilute)
            
        elseif(dtype(i).eq.1)then  !case for RV data

            tmodel(i)=0.0d0
            do 15 j=1,nintg
                tmodel(i)=tmodel(i)+K*(cos(Pid2+w+Tanom(j))+
     .              eccn*cos(Pid2+Tanom(j)))
 15         continue
            tmodel(i)=tmodel(i)/dnintg+voff
        endif
            
 10   continue
           
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function tides(M1,M2,R1,asemi,incl,tanom,eccn,Pi)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer l
      double precision M1,M2,R1,asemi,incl,tanom,eps,ra(3),ad(3),d,eccn,
     .  f(3),lambda(3),P(3),phi0,cos2incl,sin2incl,dJ(3),Pi,Inc
C     M1 - stellar mass
C     M2 - companion mass
C     R1 - stellar radius
C     asemi - semi-major axis
C     incl - orbital inclination
C     tanom - true anomoly
C     eccn - orbital eccentricity     

      tides=0.0d0  
      
      eps=M2/M1*(R1/asemi)**3.0d0
      d=asemi*(1.0d0-eccn*eccn)/(1.0d0+eccn*cos(tanom))
      phi0=0.0d0
      Inc=Pi-(incl-Pi/2.0d0)
      cos2incl=cos(Inc)*cos(Inc)
      sin2incl=sin(Inc)*sin(Inc)
            
      do 10 l=2,3
        ra(l)=(R1/asemi)**(l-2)
        ad(l)=(asemi/d)**(l+1)
        lambda(l)=dble(l)+2.0d0
 10   continue
c      lambda(2)=1.50d0
      f(2)=-1.3d1*(1.0d0+lambda(2)/4.0d0)/1.0d1
      f(3)=-5.0d0*(1.0d0+lambda(3)/1.0d1)/8.0d0
      P(2)=2.5d-1*(-(3.0d0*cos2incl-1.0d0)+
     .  3.0d0*sin2incl*cos(2.0d0*(tanom-phi0)))
      P(3)=1.25d-1*sin(Inc)*(-3.0d0*(5.0d0*cos2incl-1.0d0)*
     .  cos(tanom-phi0)+5.0d0*sin2incl*cos(3.0d0*(tanom-phi0)))
      
      do 11 l=2,3
        dJ(l)=ra(l)*ad(l)*f(l)*P(l)
        tides=tides+dJ(l)
 11   continue
      
      tides=tides*eps*Pi
c      tides=eps*cos(2.0d0*(tanom+Pi/2.0d0))
      
c      write(6,*) "TT:",tanom,tides
c      read(5,*)
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function albedomod(Pi,t,Per,ag,R1,R2,phi)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision Pi,t,Per,ag,R1,R2,phi,alpha,phase

      phi=phi+Pi
      if(phi.gt.2.0*Pi) phi=phi-2.0*Pi


      alpha=abs(phi)      
c      alpha=2.0*Pi*t/Per+phi
      alpha=alpha-2.0*Pi*int(alpha/(2.0*Pi))
      if(alpha.gt.Pi) alpha=abs(alpha-2.0*pi)
c      write(6,*) t,alpha
c      phase=(1.0d0+cos(alpha))/2.0d0
      phase=(sin(alpha)+(Pi-alpha)*cos(alpha))/Pi  !Lambertian Sphere
      
      albedomod=ag*Pi*R2*R2/(R1*R1)*phase
      
      return
      end