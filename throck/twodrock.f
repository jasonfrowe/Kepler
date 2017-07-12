      program twodrock
      implicit none
C     Computes the thermal evolution of a large 'rock' orbiting a star
      integer i,j
C     lmax is the order of the grid.
C     nmax defined the largest dimension on the grid (2**lmax)
C     niter - counts the number of models calculated  
      integer lmax,nmax,nstep
      parameter(lmax=6,nmax=2**lmax)
C     nlat and nlon define the grid size for the simulation 
      integer nlat,nlon
      parameter(nlat=nmax,nlon=2*nmax)
C     Plotting variables
      integer ia(nlon,nlat),naxes(2)
      real fitsdata(nlon,nlat),datamin,datamax,rtime,pts(nlon*nlat)
C     dtheta,dlam,lambda,theta - co-ordinates to map out grid
      double precision dtheta,dlam,lambda(nlon),theta(nlat)
C     dEn - amount of energy added or removed from an element [J]
C     T - temperature at each element [K]
C     time, dtime - time and time steps [s]
C     tmax - when time==tmax end the simulation [s]
C     sarea- surface area of elements across the surface [steradian]
C     rad - radius of the rock [m]
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
C     drad - depth of 1 layer
C     Cp - heat capacity [J kg-1 K-1]
C     rho - density [kg m-3]
C     M - mass for each element [kg]
      double precision dEn(nlat,nlon),T(nlat,nlon),dtime,time,tmax,
     .  sarea(nlat),rad,Re,PorbD,PerD,asemiAU,dsa,sangle(2),alpha,
     .  alphad,obliquity,Lstar,L,Ab,Manom,Eanom,eccn,dist,asemi,Porb,
     .  trueanomaly,distance,Tanom,Per,Srot,tcond,drad,Cp,rho,M(nlat)
C     Constants
      double precision pi,fpi,tpi,REarth,AU,day2sec,dtr,Lsun,sbolt
C     loop - logical variable to control integration loops
      logical loop

C     Open plotting device
      call opengraphics(nlon,nlat,naxes)
      
C     Constants
      pi=acos(-1.d0)   !Pi
      tpi=2.0d0*Pi !2Pi
      fpi=4.0d0*Pi !4Pi
      dtr=pi/180.0d0 !radians/degree
      REarth=6.371*10d6 !Earth radius [m]
      AU=1.4959787069d11 !Astronomical Unit [m]
      day2sec=24.0d0*60.0d0*60.0d0 !sec/Day
      Lsun=3.839d26 !Solar luminosity [W]
      sbolt=5.670d-8 !Stefan-Boltzmann [W m-2 K-4]

C     Model parameters 
      dtime=30.0d0 !set time steps
      tmax=100.0*365.25*24.0*60.0*60.0 !stop model at this time
      Re=1.0d0 !radius of rock in Earth Radii
      PorbD=365.25d0 !orbital period of planet in days
      PerD=1.0d0 !rotational period of planet in days
      asemiAU=1.0d0 !semi-major axis [AU]
      eccn=0.016710219d0 !eccentricity 
      obliquity=23.439281d0 !obliquity of planet rotation
      alphad=90.0d0 !angle of rotated grid [degrees] - leave alone
      Lstar=1.0d0 !Luminosity of star [Lsun]
      Ab=0.367 !Bond Albedo of planet surface
      tcond=1.1d0 !thermal conductivity [W m-1 K-1]
      drad=1.0d0 !depth of layer [m]
      Cp=840.0d0 !heat capacity [J kg-1 K-1]
      rho=2700.0d0 !density [kg m-3]
C     Change of Units
      rad=REarth*Re !Convert Re to m
      Porb=PorbD*day2sec !Convert days to sec
      Per=PerD*day2sec !Convert days to sec
      asemi=asemiAU*AU !Convert AU to m
      alpha=alphad*dtr !Convert degress to radians
      Srot=1.0d0/Per-1.0d0/Porb !inverse siderial rotation rate [s-1]
      L=Lstar*Lsun !Convert Lsun to W

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
c        write(6,*) theta(i)
 15   continue
c      read(5,*)
      
C     Get surface area for each element as a function of latittude[ster]
      call surfacearea(nlat,theta,dlam,sarea,Pi)
C     Get mass for each element as a function of latittude [kg]
      do 18 i=1,nlat
        M(i)=rho*sarea(i)*drad*rad*rad
 18   continue


C     initialize T to zero
      do 12 i=1,nlat
        do 13 j=1,nlon
            T(i,j)=0.0d0
 13     continue
 12   continue
 
C     Initialize Eanom to first guess for solving Kepler
      Eanom=pi
C     If eccn=0, then dist==asemi
      dist=asemi
C     initialize nstep to counter number of models
      nstep=0
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
                dEn(i,j)=0.0
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
C       Obliquity angle.      
        alphad=90.0d0+obliquity*sin(tpi*time/Porb) !tilt angle of planet
        alpha=dtr*alphad
        sangle(1)=alpha

C       Calculate deposited energy from star
c        if(time.lt.864000.0)then
        call stellarrad(nlat,nlon,lambda,theta,dEn,L,dist,Ab,pi,fPi,
     .      sarea,rad,dtime,sangle)
     
C       Calculate energy radiated away by the rock
        call radiate(nlat,nlon,sbolt,dtime,sarea,rad,T,dEn)
c        endif
        
C       Calculate energy redistribution from conduction
        call conduction(nlat,nlon,tcond,rad,drad,dtime,dtheta,dlam,
     .      theta,dEn,T)
     
        call heatcapacity(nlat,nlon,Cp,T,dEn,M)
     
        if(mod(nstep,100).eq.0)then
            datamin= 99.9e30
            datamax=-99.9e30
            do 16 i=1,nlat
                do 17 j=1,nlon
c                    fitsdata(j,i)=dEn(i,j)
                    fitsdata(j,i)=T(i,j)
                    datamin=min(fitsdata(j,i),datamin)
                    datamax=max(fitsdata(j,i),datamax)
 17             continue
 16         continue
            rtime=real(time) !convert real*8 to real*4
            call displayfits(nlon,nlat,ia,pts,fitsdata,naxes,datamin,
     .          datamax,rtime)
        endif
        
c        read(5,*)
        
        time=time+dtime !increase time for next step
        if(time.gt.tmax) loop=.false.
      enddo
      
C     Close plotting device
      call closegraphics()
      
      end