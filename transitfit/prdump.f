C     Used to dump data for online animations.
C     (c) Jason F. Rowe 2010
      program prdump
      implicit none
      integer iargc,nunit,nmax,npt,nfit,i,nsample,today(3)
      parameter(nmax=650000,nfit=18)
      integer dtype(nmax)
      double precision time(nmax),mag(nmax),merr(nmax),itime(nmax),
     .  MOSTtime,sol(nfit),serr(nfit,2),Dpvary(nfit),err(nfit),doe,toff,
     .  tmodel(nmax),meananom(nmax),tdur,transitdur,eccn,w,Pi,tPi,m1,m2,
     .  stime,dnintg,dnintgm1
      character*80 obsfile, inputsol
      
      Pi=acos(-1.d0)   !Pi
      tPi=2.0d0*Pi

      call idate(today)

      if(iargc().lt.1) goto 901 !check number of command line arguements
      call getarg(1,obsfile) !get filename for input solution
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,mag,merr,itime,MOSTtime)
      call readkeplc(nunit,nmax,npt,time,mag,merr,itime,MOSTtime)
      do 10 i=1,npt !central database of all data
        dtype(i)=0 !0 marks that we have photometric data
 10   continue
      close(nunit)!release unit number as we are done with the file.
      
      if(iargc().ge.2) then
        call getarg(2,inputsol) !get filename for input solution
        nunit=10 !unit number used for file input
        open(unit=nunit,file=inputsol,status='old',err=902)
C       We start by reading in solution from input file
        call getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff)
        close(nunit) !release unit number as we are done with file
      else
C     Import a default solution set to start with
        call defaultpars(nfit,sol,serr,Dpvary,err,doe,toff)
      endif
      
      eccn=sqrt(sol(14)*sol(14)+sol(15)*sol(15)) !eccentricity
      if(eccn.ge.1.0) eccn=0.99
      if(eccn.eq.0.0d0)then
        w=0.0d0
      else
        w=atan(sol(15)/sol(14))
        if((sol(14).gt.0.0d0).and.(sol(15).lt.0.0d0))then
            w=tPi+w
        elseif((sol(14).lt.0.0d0).and.(sol(15).ge.0.0d0))then 
            w=Pi+w
        elseif((sol(14).lt.0.0d0).and.(sol(15).lt.0.0d0))then
            w=Pi+w
        endif
      endif
      
      tdur=transitdur(nfit,sol)/86400.0d0
      write(0,*) "Tdur: ",tdur

      m1=(-tdur)/sol(5)-floor((-tdur)/sol(5))
      m1=m1*tPi+w
      if(m1.gt.tPi) m1=m1-tPi
      if(m1.lt.0.0d0) m1=m1+tPi

      m2=(tdur)/sol(5)-floor((tdur)/sol(5))
      m2=m2*tPi+w
      if(m2.gt.tPi) m2=m2-tPi
      if(m2.lt.0.0d0) m2=m2+tPi     
      
      write(0,*) "m1,m2:",m1,m2
      
      call transitmodel(npt,time,itime,dtype,meananom,tmodel,nfit,sol)
      
c      write(6,500) "# Kepler Observations, N=",npt
 500  format(A25,I6)
      write(6,501) today(2),today(1),2000+today(3)
 501  format('<!-- data for Kepler XXXb, generated ',2(i2.2,'/'),i4.4,
     .  ' -->')
      write(6,504) '<dataPoints>'
 504  format(A12)

      do 11 i=1,npt
c        write(0,*) meananom(i),10**(mag(i)/-2.5d0)
c        if((meananom(i).gt.m1).and.(meananom(i).lt.m2))
        if((meananom(i).gt.m1).or.(meananom(i).lt.m2))
     .      write(6,502) meananom(i),10**(mag(i)/-2.5d0)
 502        format('<pt><ma>',F10.8,'</ma><i>',F10.8,'</i></pt>')
 11   continue
      write(6,505) '</dataPoints>'
 505  format(A13)
      
      nsample=200
      dnintg=dble(nsample)
      dnintgm1=2.0*dnintg-2.0
      stime=tdur*2.0d0
      write(0,*) "stime:",stime,tdur
      do 12 i=1,nsample
        time(i)=sol(7)+stime*(2.0*dble(i)-dnintg-1.0)/dnintgm1
c        write(0,*) time(i),stime*(2.0*dble(i)-dnintg-1.0)/dnintgm1
 12   continue
      
      npt=nsample
      call transitmodel(npt,time,itime,dtype,meananom,tmodel,nfit,sol)
      
c      write(6,500) "# Model Observations,  N=",nsample
      write(6,506) '<curvePoints>'
 506  format(A13)
      do 13 i=1,npt
        if((meananom(i).gt.m1).or.(meananom(i).lt.m2))
     .      write(6,503) meananom(i),10**(tmodel(i)/-2.5d0)
 503        format('<pt><ma>',F10.8,'</ma><i>',F10.8,'</i></pt>')
 13   continue 
      write(6,507) '</curvePoints>'
 507  format(A14)      
       
      goto 999
 901  write(0,*) "Usage: modelbuild <obsfile> <modelpars>"
      write(0,*) " <obsfile>: file containing photometry"
      write(0,*) " <modelpars>: model parameters to adjust"
      write(0,*) "              omit to use defaults"
      goto 999
 902  write(0,*) "Cannot open ",inputsol
      goto 999
 903  write(0,*) "Cannot open ",obsfile
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function transitdur(nfit,sol)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This is for circular orbits only.
      implicit none
      integer nfit
      double precision sol(nfit),Psec,M1,M2,R1,R2,asemi,temp(4),incl
C     Read in physical constants (Pi,G,Msun,..etc.)      
      include "utils/physcons.f"
      
      M1=sol(1)*Msun !kg ; mass of star
      M2=sol(2)*Mjup !kg ; mass of planet
      Psec=sol(5)*24.0*60.0*60.0 !sec ; period of planet
      asemi=(Psec*Psec*G*(M1+M2)/(4.0*Pi*Pi))**(1.0/3.0) !m
      R1=sol(3)*Rsun  !radius of star
      R2=sol(4)*Rjup  !radius of planet
      incl=Pi*sol(6)/180.0d0
      
      temp(1)=Psec/Pi
      temp(2)=R1/asemi
      temp(3)=(1+(R2/R1))**2.0-((asemi/R1)*cos(incl))**2.0
      temp(4)=1-cos(incl)*cos(incl)     
      
      transitdur=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))
      
      return
      end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine transitmodel(npt,time,exptime,dtype,meananom,tmodel,
     .  nfit,sol)
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
     .  incl,Per,phi(nintg),eccn,w,ted,inclmin,epoch,phi1
C     Parameters that work towards MandelAgol
      double precision x1,y1,t(nintg),dnintg,dnintgm1,x2(nintg),
     .  y2(nintg),phase,tpi,angle,darea,b0(nintg),mu(nintg),mandelagol,
     .  mulimb0(nintg),mulimbf(nintg,5),dist(nintg),et(2),etest,phi0
C     Albedo routines
      double precision albedomod,zarea,xintold,y2pold,eclmod,y2pold2,
     .  xintold2,eclmod2,test,tides,tflux
      double precision meananom(npt)
C     RV parameters
C     K - amplitude of RV
C     voff - radial velocity offset (m/s)
C     vmodel - model velocities (m/s)
      integer nptv
      double precision K,voff,Kc,Cs,fDB,fT
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
      
      Cs=2.99792458e8 !Speed of light
      fDB=1.896 !Doppler Boosting factor
      fT=3.37 !ellisodial mass factor

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
        if((sol(14).gt.0.0d0).and.(sol(15).lt.0.0d0))then
            w=tPi+w
        elseif((sol(14).lt.0.0d0).and.(sol(15).ge.0.0d0))then 
            w=Pi+w
        elseif((sol(14).lt.0.0d0).and.(sol(15).lt.0.0d0))then
            w=Pi+w
        endif
      endif
c      write(0,*) "e,w",eccn,w
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
c        if(incl.lt.inclmin) incl=inclmin
        if(incl.gt.90.0d0)incl=180.0-incl    
c        if(incl.lt.inclmin) incl=inclmin !check a second time
      endif
      incl=Pi*(90.0d0-incl)/180.0d0 !radians
      
C     observer-star-planet angle
C     Find phase at centre of transit
      epoch=sol(7)
      Manom=w !mean anomaly
      call kepler(Manom,Eanom,eccn) !eccentric anomaly
      phi0=trueanomaly(eccn,Eanom)
c      write(6,*) eccn,w,phi0

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
c      open(unit=21,file="xy.dat")
c      write(21,*) 0.0,0.0
      do 10 i=1,npt  !loop to get models points of each time
        meananom(i)=((time(i)-epoch)/per-
     .      floor((time(i)-epoch)/per))*tPi+w
        if(meananom(i).gt.tPi) meananom(i)=meananom(i)-tPi
        if(meananom(i).lt.0.0d0) meananom(i)=meananom(i)+tPi
        do 11 j=1,nintg  !get array of times spanning exposure time
C           These times are centered on time(i)
            t(j)=time(i)+exptime(i)*(2.0*dble(j)-dnintg-1.0)/dnintgm1-
     .          epoch
            phi(j)=t(j)/per-floor(t(j)/per)
            phi(j)=phi(j)*tPi!-phi0      
c            if(phi(j).gt.tPi)   phi(j)=phi(j)-tPi
c            if(phi(j).lt.0.0d0) phi(j)=phi(j)+tPi 
            Manom=phi(j)+w
            if(Manom.gt.tPi) Manom=Manom-tPi
            if(Manom.lt.0.0d0) Manom=Manom+tPi
c            meananom(i)=meananom(i)+Manom
            call kepler(Manom,Eanom,eccn)
            Tanom(j)=trueanomaly(eccn,Eanom) !get true anonaly
c            phi(j)=Tanom(j)!+w-phi0  !compute phase relative to observer
            if(phi(j).gt.Pi) phi(j)=phi(j)-tPi            
            arad(j)=distance(asemi,eccn,Tanom(j))
c            x2(j)=arad(j)*Sin(phi(j))
c            y2(j)=arad(j)*Cos(phi(j))*Sin(incl)
c            write(21,*) arad(j)*Sin(phi(j)),arad(j)*Cos(phi(j)),Tanom(j)
            x2(j)=arad(j)*Sin(Tanom(j)-phi0)
            y2(j)=arad(j)*Cos(Tanom(j)-phi0)*Sin(incl)
c            write(21,*) arad(j)*Sin(Tanom(j)-phi0),
c     .          arad(j)*Cos(Tanom(j)-phi0),Tanom(j)
c            write(0,*) asemi,eccn,Tanom(j)
c            read(5,*)
cw            x2(j)=arad(j)*Sin(Tanom(j))
cw            y2(j)=arad(j)*Sin(incl)*Cos(Tanom(j))
c            write(0,*) t(j),Tanom,Manom
c            read(5,*)
c         write(6,*) t(j),x2(j)/R1
 500     format(28(F7.4,1X))
 11     continue
c        meananom(i)=meananom(i)/dnintg
c        write(6,*) t(1),x2(1),y2(1)
c        read(5,*)
 
c        Manom=tpi*time(i)/Per+phi !mean anomaly
c        if(Manom.ge.tpi) Manom=Manom-tpi*int(Manom/tpi)
c        call kepler(Manom,Eanom,eccn)
c        Tanom=trueanomaly(eccn,Eanom)
 
        if(dtype(i).eq.0)then !case for photometric data
            zarea=0.0!initialize zarea
            tflux=0.0
            do 12 j=1,nintg
C               Doppler Boosting
c                Kc=K*(cos(Pid2+Tanom(j)-phi0)+eccn*cos(w))
                Kc=-K*(cos(Pid2+Tanom(j)-phi0)+eccn*cos(w))
                tflux=tflux+fDB*Kc/Cs
c                write(0,*) "Kc: ",Kc,fDB*Kc/Cs
c                read(5,*)
            
C               ellipsidal component
c                tflux=tides(M1,M2,R1,asemi,incl,phi(j),eccn,Pi)
                tflux=tflux+
     .              tides(M1,M2,R1,asemi,incl,Tanom(j)-phi0,eccn,Pi,fT)
C               flux from star + planet
                Ag=sol(9)*R1*R1/(arad(j)*arad(j))  !albedo of planet
c                zarea=zarea+albedomod(Pi,t(j),Per,Ag,R1,R2,phi(j))
                zarea=zarea+
     .              albedomod(Pi,t(j),Per,Ag,R1,R2,Tanom(j)-phi0)
cw              zarea=zarea+albedomod(Pi,t(j),Per,Ag,R1,R2,Tanom(j))
 12         continue
            tflux=tflux/dnintg
            zarea=zarea/dnintg


            phase=Tanom(nintg/2+1)-phi0!phi(nintg/2+1)
            if(phase.gt.Pi) phase=phase-tPi
            if(phase.lt.-Pi) phase=phase+tPi
            if(abs(phase).lt.Pid2)then
C     Now we calculate the flux change from planet transiting the star 
                darea=(Pi*mandelagol(nintg,R1,R2,x1,x2,y1,y2,c,b0,mu,
     .              mulimb0,mulimbf,dist)+zarea+tflux)/norm
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
                darea=darea/dnintg+tflux/norm !average and add on tidal
                xintold=xintold2
                y2pold=y2pold2
            endif
c           write(6,*) "hello4",i
C     Convert relative fluxes to magnitude to match observations
            tmodel(i)=zpt-2.5*log10(darea+(1.0-darea)*dilute)
c            write(6,*) meananom(i),tmodel(i)
        elseif(dtype(i).eq.1)then  !case for RV data

            tmodel(i)=0.0d0
            do 15 j=1,nintg
                tmodel(i)=tmodel(i)+K*(cos(Pid2+Tanom(j)-phi0)+
     .              eccn*cos(w))
c                tmodel(i)=tmodel(i)+K*(cos(Pid2+phi(j))+
c     .              eccn*cos(w))
 15         continue
            tmodel(i)=tmodel(i)/dnintg+voff
        endif
            
 10   continue
c      close(21)
           
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function tides(M1,M2,R1,asemi,incl,tanom,eccn,Pi,
     .  fT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer l
      double precision M1,M2,R1,asemi,incl,tanom,eps,ra(3),ad(3),d,eccn,
     .  f(3),lambda(3),P(3),phi0,cos2incl,sin2incl,dJ(3),Pi,Inc,fT
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

c      write(0,*) incl,Inc
c      read(5,*)
      tides=tides*eps*Pi
c      tides=-fT*eps*cos(2.0d0*(tanom))*(sin(Pi/2.0d0-incl))**3

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