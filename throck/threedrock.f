      program threedrock
      implicit none
C     Computes the thermal evolution of a large 'rock' orbiting a star
      integer i,j,k,iargc,nunit
C     lmax is the order of the grid.
C     nmax defined the largest dimension on the grid (2**lmax)
C     niter - counts the number of models calculated  
      integer lmax,nmax,nstep,nplot
      parameter(lmax=5,nmax=2**lmax)
C     nlat and nlon define the grid size for the simulation
C     nR defines the number of shells 
      integer nlat,nlon,nR
      parameter(nlat=nmax,nlon=2*nmax,nR=nmax/2)
C     Plotting variables
      integer ia(nlon,nlat),naxes(2)
      real fitsdata(nlon,nlat),datamin,datamax,rtime,pts(nlon*nlat),
     .  px(nR),py(nR),pxold(nR),pyold(nr),py2(nR),pyold2(nr),py3(nR),
     .  pyold3(nr)
C     dtheta,dlam,lambda,theta - co-ordinates to map out grid
      double precision dtheta,dlam,lambda(nlon),theta(nlat)
C     dEn - amount of energy added or removed from an element [J]
C     T - temperature at each element [K]
C     time, dtime - time and time steps [s]
C     tmax - when time==tmax end the simulation [s]
C     sarea- surface area of elements across the surface [steradian]
C     rad0 - radius of the rock [m]
C     Manom - Mean anomaly [radian]
C     Eanom - Eccentric anomaly [radian]
C     Tanom - True anomaly [radian]
C     eccn - eccentricity (bound orbits only!!) [0-1]
C     dist - distance between planet and star [m]
C     asemi - semi-major axis of planet-star orbit [m]
C     Porb - orbital period of the planet [s]
C     Per - rotational period of the planet [s]
C     trueanomaly,distance - functions (see keplerian.f)
C     sangle - position of star in sky (lambda,theta) [radians]
C     obliquity- obliquity of planet rotation [degrees]
C     Ab - Bond albedo of the planet surface
C     tcond - thermal conductivity [W m-1 K-1]
C     rdepth - the total depth of the model
C     drad - depth of 1 layer [m]
C     Cp - heat capacity [J kg-1 K-1]
C     rho - density [kg m-3]
C     M - mass for each element [kg]
C     tdr - temperature change per unit depth to get core temp [K km-1]
C     Teq - equilibrium temperature of planet [K]
C     Tstar - effective temperature of star [K]
C     Tcore - core temeperature (base of model) [K]
      double precision dEn(nlat,nlon,nR),T(nlat,nlon,nR),dtime,time,
     .  tmax,sarea(nlat),rad(nR),Re,PorbD,PerD,asemiAU,dsa,sangle(2),
     .  alpha,alphad,obliquity,Lstar,L,Ab,Manom,Eanom,eccn,dist,asemi,
     .  Porb,trueanomaly,distance,Tanom,Per,Srot,drad,rho,
     .  M(nlat,nR),rad0,rdepth,tdr,Teq,Tstar,Tcore,Rsun,Eint,
     .  Cp(nlat,nlon,nR),tcond(nlat,nlon,nR),Rstar,R
C     Variables for calculating observed flux
      integer nKepler,nT
      parameter(nKepler=491,nT=3000)
      double precision Keplerlam(nKepler),Keplerpass(nKepler),planck,
     .  st_flux,Tgridmin,Tgridmax,Temps(nT),Fluxes(nT),yp1,yp2,yT(nT),
     .  pflux,obsflux,vangle(2),aveT
C     Constants
      double precision pi,fpi,tpi,REarth,AU,day2sec,dtr,Lsun,sbolt,pid2
C     parameter input file
      character*80 parfile
C     loop - logical variable to control integration loops
      logical loop

      if(iargc().lt.1) goto 901

C     Open plotting device
      call opengraphics(nlon,nlat,naxes)
      
C     Constants
      pi=acos(-1.d0)   !Pi
      tpi=2.0d0*pi !2Pi
      fpi=4.0d0*pi !4Pi
      pid2=pi/2.0d0
      dtr=pi/180.0d0 !radians/degree
      REarth=6.371*10d6 !Earth radius [m]
      AU=1.4959787069d11 !Astronomical Unit [m]
      day2sec=24.0d0*60.0d0*60.0d0 !sec/Day
      Lsun=3.839d26 !Solar luminosity [W]
      sbolt=5.670d-8 !Stefan-Boltzmann [W m-2 K-4]
      Rsun=696265.0*1000.0 !m  radius of Sun

      call getarg(1,parfile)
      nunit=10
      call getpars(nunit,parfile,dtime,tmax,Re,PorbD,PerD,asemiAU,eccn,
     .  obliquity,alphad,Lstar,Rstar,Ab,rdepth,rho,tdr,Tstar,Eint)

C     Change of Units
      rad0=REarth*Re !Convert Re to m
      Porb=PorbD*day2sec !Convert days to sec
      Per=PerD*day2sec !Convert days to sec
      asemi=asemiAU*AU !Convert AU to m
      alpha=alphad*dtr !Convert degress to radians
      Srot=1.0d0/Per-1.0d0/Porb !inverse siderial rotation rate [s-1]
      L=Lstar*Lsun !Convert Lsun to W
      R=Rstar*Rsun !Convert Rsun to m

C     Readin Kepler bandpass for Flux measurements
      call readKeplerbandpass(nKepler,Keplerlam,Keplerpass)
C     Calculate the flux contribution from the star
      st_flux=Planck(Tstar,nKepler,Keplerlam,Keplerpass)
C     1.0d-3 is to convert to W m^-2
      st_flux=st_flux*fPi*R*R*1.0d-3 ![W m^-2]

C     Precompute fluxes in expected temperature range, the interpolate
C     the results to increase computation speed.
      Tgridmin=30.0
      Tgridmax=4050.0 
      call genfluxes(nT,Tgridmin,Tgridmax,Temps,Fluxes,nKepler,
     .  Keplerlam,Keplerpass)
      yp1=1.0d30
      yp2=1.0d30
      call spline(Temps,Fluxes,nT,yp1,yp2,yT)

C     Starting position of star in sky
      sangle(1)=alpha !lambda
      sangle(2)=pi    !theta
      dsa=tpi*dtime*Srot !change in sun angle across nlon 
      
C     Angles across the surface
      dtheta = pi/(nlat+1)
      dlam = tpi/nlon      
      do 14 i=1,nlon
        lambda(i) = (i-1)*dlam
 14   continue
      do 15 i=1,nlat
        theta(i) = i*dtheta
 15   continue
      
C     Get surface area for each element as a function of latittude[ster]
      call surfacearea(nlat,theta,dlam,sarea,Pi)

      drad=rdepth/dble(nR-1)
      do 19 i=1,nR
        rad(i)=rad0-drad*real(i-1)
        if(rad(i).lt.0.0d0) write(0,*) "Hey... inverted radius?"
 19   continue

C     Get mass for each element as a function of latittude [kg]
      do 18 i=1,nlat
        do 22 k=1,nR
            M(i,k)=rho*sarea(i)*drad*rad(k)*rad(k)
 22     continue
 18   continue

C     Equilibrium temperature of the planet
      Teq=Tstar*(R/(2.0d0*asemi))**0.50d0*(1.0d0-Ab)**0.25d0
      write(0,*) "Teq: ",Teq
      Tcore=0.0!tdr*rdepth/1000.0 !core temperature [K]
 
C     initialize boundary conditions 
c      do 23 j=1,nlon
c        do 24 i=1,nlat
c            T(i,j,nR)=Tcore
c 24     continue
c 23   continue
            
C     initialize T to Teq
      do 12 i=1,nlat
        do 13 j=1,nlon
            do 20 k=1,nR
                T(i,j,k)=Teq
 20         continue
 13     continue
 12   continue
 
C     Initialize Eanom to first guess for solving Kepler
      Eanom=pi
C     If eccn=0, then dist==asemi
      dist=asemi
C     initialize nstep to counter number of models
      nstep=0
C     initialize nplot to count number of plots made
      nplot=0
C     initialize time for step 1
      time=0.0d0
C     initialize loop variable
      loop=.true.
C     the main loop
      do while(loop)
        nstep=nstep+1 !increase counter
C       Reset dEn to zero
        do 10 i=1,nlat
            do 11 j=1,nlon
                do 21 k=1,nR
                    dEn(i,j,k)=0.0
 21             continue
 11         continue
 10     continue

        if(eccn.gt.0.0)then
            Manom=tpi*(time)/Porb+pi !mean anomaly
            if (Manom.ge.tpi) Manom=Manom-tpi*int(Manom/tpi)
C           Now get distance and true anomaly for Keplerian orbits
            call kepler(Manom,Eanom,eccn)
            Tanom=trueanomaly(eccn,Eanom) !true anomaly 
            dist=distance(asemi,eccn,Tanom) !distance
        endif
        
C       Update position of star in sky.
        sangle(2)=sangle(2)-dsa
        if(sangle(2).gt.tpi) sangle(2)=sangle(2)-tpi
        if(sangle(2).lt.0.0) sangle(2)=sangle(2)+tpi       
C       Obliquity angle.  tilt angle of planet
        alphad=90.0d0+obliquity*sin(tpi*time/Porb-pi/4.0d0) 
        alpha=dtr*alphad
        sangle(1)=alpha
        
C       Update viewing angle
        vangle(2)=-tpi*(time/Per-int(time/per))
        vangle(1)=pid2

C       Calculate deposited energy from star
c        if(time.lt.864000.0)then
        call stellarrad(nlat,nlon,nR,lambda,theta,dEn,L,dist,Ab,pi,
     .      fPi,sarea,rad,dtime,sangle)
C       Calculate energy radiated away by the rock
        call radiate(nlat,nlon,nR,sbolt,dtime,sarea,rad,T,dEn)
c        endif
 
C       Update global values of tcond and Cp
        call kcp(nlat,nlon,nR,tcond,Cp,T,rho)
        !write(0,*) "Cp: ",Cp(1,1,1),tcond(1,1,1)
        
C       Calculate energy redistribution from conduction
        call conduction(nlat,nlon,nR,tcond,rad,drad,dtime,dtheta,dlam,
     .      theta,dEn,T,sarea,Eint)
     
        call heatcapacity(nlat,nlon,nR,Cp,T,dEn,M)
        
        if(mod(nstep,100).eq.0)then
            pflux=obsflux(nlat,nlon,nR,nT,vangle,theta,lambda,T,
     .          Tgridmin,Tgridmax,pi,pid2,tpi,Temps,Fluxes,yT,sarea,
     .          rad0,aveT)     
            write(6,500) time/day2sec,log10(pflux/st_flux),aveT
            write(0,501) time/day2sec,log10(pflux/st_flux),aveT
 500        format(1026(1PE16.9,1X))
 501        format(F9.3,1X,F6.2,1X,F9.3)
        endif
        
        if(mod(nstep,10).eq.0)then
            nplot=nplot+1
            datamin= 99.9e30
            datamax=-99.9e30
            do 16 i=1,nlat
                do 17 j=1,nlon
c                    fitsdata(j,i)=dEn(i,j,1)
                    fitsdata(j,i)=T(i,j,1)
                    datamin=min(fitsdata(j,i),datamin)
                    datamax=max(fitsdata(j,i),datamax)
 17             continue
 16         continue
            rtime=real(time) !convert real*8 to real*4
            call displayfits(nlon,nlat,ia,pts,fitsdata,naxes,datamin,
     .          datamax,rtime)
C           Plot on the log scale
            do 25 i=1,nlat
                do 26 j=1,nlon
                    fitsdata(j,i)=log10(fitsdata(j,i)-datamin+1.0)
 26             continue
 25         continue
            datamax=log10(datamax-datamin+1.0)
            datamin=0.0      
            call avetemp(nlat,nlon,nR,nplot,T,px,py,pxold,pyold,drad,
     .          Teq,Tcore,sarea,fpi,py2,pyold2,py3,pyold3)
        endif
        
c        read(5,*)
        
        time=time+dtime !increase time for next step
        if(time.gt.tmax) loop=.false.
      enddo
      
C     Close plotting device
      call closegraphics()
      
      goto 999
 901  write(0,*) "Usage: threerock <parsfile>"
      write(0,*) " <parsfile> : contains model parameters"
      goto 999
 999  end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine avetemp(nlat,nlon,nR,nplot,T,px,py,pxold,pyold,drad,
     .          Teq,Tcore,sarea,fpi,py2,pyold2,py3,pyold3)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nlat,nlon,nR,i,j,k,nplot
      real px(nR),py(nR),pxold(nR),pyold(nR),rn,bounds(4),py2(nR),
     .  pyold2(nR),py3(nR),pyold3(nR)
      double precision T(nlat,nlon,nR),meanT,drad,Teq,Tcore,sarea(nlat),
     .  fpi,Tmin,Tmax
    
      rn=real(fpi)
      
      do 10 k=1,nR
        meanT=0.0
        Tmin=T(1,1,k)
        Tmax=T(1,1,k)
        do 11 j=1,nlon
            do 12 i=1,nlat
                meanT=meanT+T(i,j,k)*sarea(i)
                Tmin=min(Tmin,T(i,j,k))
                Tmax=max(Tmax,T(i,j,k))
 12         continue
 11     continue
        px(k)=real(k-1)*real(drad) ![m]
        py(k)=real(meanT)/rn ![K]
        py2(k)=real(Tmin)
        py3(k)=real(Tmax)
 10   continue
 
      bounds(1)=px(1)
      bounds(2)=px(nR)
      bounds(3)=100.0d0
      bounds(4)=3000.0!max(1.25d0*Teq,1.25d0*Tcore)
      call pgsvp(0.10,0.45,0.65,0.95)
      call pgwindow(bounds(1),bounds(2),bounds(3),bounds(4))
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("depth (m)","Temp (K)","")
      if(nplot.gt.1)then
        call pgsci(0)
        call pgline(nR,pxold,pyold)
        call pgline(nR,pxold,pyold2)
        call pgline(nR,pxold,pyold3)
        call pgsci(1)
      endif
      call pgsci(5)
      call pgline(nR,px,py2)
      call pgsci(2)
      call pgline(nR,px,py3)
      call pgsci(1)
      call pgline(nR,px,py)
      
c      write(0,*) "Tb:",py(nR),py(1)
      
      do 13 i=1,nR
        pxold(i)=px(i)
        pyold(i)=py(i)
        pyold2(i)=py2(i)
        pyold3(i)=py3(i)
 13   continue
 
      return
      end
 
 
    
    
    
