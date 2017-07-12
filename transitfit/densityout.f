      program phys2geo
      implicit none
C     Convert a Physical space model to a geometric model
      integer nfitg,nfit,nunit,iargc,i
      parameter(nfitg=18,nfit=18)
      double precision solg(nfitg),serrg(nfitg,2),Dpvaryg(nfitg),
     .  errg(nfitg),doe,toff,sol(nfit),serr(nfit,2),err(nfit),Psec,
     .  M1,M2,R1,aConst,asemi,Pi,Pid2,G,sb,Msun,Mearth,Mjup,Rjup,Rsun,
     .  Lsun,incl,R2,eccn,tPi,K,AU,b,rdrs,adrs
      character*3 titles(nfit)
      character*80 inputsol      

      if(iargc().lt.1) goto 901
      
C     Constants
      Pi=acos(-1.d0)   !Pi
      Pid2=Pi/2.0d0
      tPi=Pi*2.0d0
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      sb=5.6704d-8 !W m^-2 K^-4 Stefan-Boltzman constant
      Msun=1.9891d30 !kg  mass of Sun
      Mearth=5.974d24 !kg mass of Earth
      Mjup=317.833d0*Mearth !kg  mass of Jupiter  
      Rjup=142980.0d0*1000.0d0/2.0d0 !m  radius of Jupiter     
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      Lsun=3.839d26 !W Solar Luminosity  
      AU=1.49598e11    
      
      call getarg(1,inputsol)
      nunit=10
      open(unit=nunit,file=inputsol,status='old',err=901)
      call getfitpars(nunit,nfitg,solg,serrg,Dpvaryg,errg,doe,toff)
      
      do 10 i=1,nfit
        serr(i,1)=0.0d0
        serr(i,2)=0.0d0
        err(i)=0.0d0
 10   continue
      do 11 i=1,6
        serr(i,2)=-1.0d0
 11   continue
      serr(15,2)=-1.0d0
      
      Psec=solg(5)*8.64d4 !sec ; period of planet
      M1=solg(1)*Msun !kg ; mass of star
      M2=solg(2)*Mjup !kg ; mass of planet
      R1=solg(3)*Rsun !m ; stellar radius
      R2=solg(4)*Rjup !m ; planet radius
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)
      asemi=(M1+M2)**(1.0d0/3.0d0)*Psec**(2.0d0/3.0d0)*aConst
c      write(6,*) "asemi: ",asemi/AU,asemi/R1
      incl=solg(6)
      if(incl.gt.90.0d0)incl=180.0-incl     
      incl=Pi*(incl)/180.0d0
      
      eccn=sqrt(solg(14)*solg(14)+solg(15)*solg(15)) !eccentricity
      if(eccn.ge.1.0) eccn=0.99
c      if(eccn.eq.0.0d0)then
c        w=0.0d0
c      else
c        w=atan(sol(8)/sol(7))
c        if((sol(14).gt.0.0d0).and.(sol(15).lt.0.0d0))then
c            w=tPi+w
c        elseif((sol(14).lt.0.0d0).and.(sol(15).gt.0.0d0))then 
c            w=Pi+w
c        elseif((sol(14).lt.0.0d0).and.(sol(15).lt.0.0d0))then
c            w=Pi+w
c        endif
c      endif

c      K=2.0*pi*G*M2**3*(sin(incl))**3/
c     .  (Psec*(1.0d0-eccn*eccn)**(3.0d0/2.0d0)*(M1+M2)*(M1+M2))
c      K=K**(1.0d0/3.0d0)
c      write(6,*) "K: ",K
      
      b=(asemi*cos(incl)/R1)**2.0
      rdrs=R2/R1
      adrs=asemi/R1
      
      write(6,500) b,b*0.30,rdrs,rdrs*errg(4)/solg(4),adrs,adrs*0.30
 500  format(2(F6.4,1X),2(F7.5,1X),2(F10.6,1X))
      
c      sol(1)=solg(7)  !Epoch
cc      write(6,*) sol(1)
c      sol(2)=solg(5)  !Period
cc      write(6,*) sol(2)
c      sol(3)=(asemi*cos(incl)/R1)**2.0 !b^2 
cc      write(6,*) sol(3) 
c      sol(4)=R2/R1  !R/R*
cc      write(6,*) sol(4)
c      sol(5)=asemi/R1*tPi/solg(5)/sqrt(1-sol(3))/(1-solg(15))*
c     .  sqrt(1-eccn*eccn) !z/R*
cc      sol(5)=asemi/R1*tPi/solg(5)/sqrt(1-sol(3))/(1-solg(15))
cc      write(6,*) sol(5)
c      sol(6)=solg(8) !zpt
cc      write(6,*) sol(6)
c      sol(7)=solg(14) !ecw
cc      write(6,*) sol(7)
c      sol(8)=solg(15) !esw
cc      write(6,*) sol(8)
cc      sol(9)=K !radial velocity amplitude
cc      write(6,*) sol(9)
cc      sol(10)=solg(17) !Velocity zero point
cc      write(6,*) sol(10)
c      sol(11)=solg(10)
c      sol(12)=solg(11)
c      sol(13)=solg(12)
c      sol(14)=solg(13)
cc      write(6,*) sol(11),sol(12),sol(13),sol(14)
c      sol(15)=solg(16) !TED
cc      write(6,*) sol(15)
c      sol(16)=0.0 !ellipsoidal variations 
c      sol(17)=0.0 !Phase changes of planet (Albedo)
c      sol(18)=solg(18) !Dilution
cc      write(6,*) sol(16)
cc      write(6,*) sol(17)
cc      write(6,*) sol(18)
      
      goto 999
 901  write(0,*) "Cannot open" ,inputsol
      goto 999
 999  end
      