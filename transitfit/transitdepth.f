      program transitdepth
      implicit none
C     Convert a Physical space model to a geometric model
      integer nfitg,nfit,nunit,iargc,i,npt
      parameter(nfitg=18,nfit=18)
      double precision solg(nfitg),serrg(nfitg,2),Dpvaryg(nfitg),
     .  errg(nfitg),doe,toff,sol(nfit),serr(nfit,2),err(nfit),Psec,
     .  M1,M2,R1,aConst,asemi,Pi,Pid2,G,sb,Msun,Mearth,Mjup,Rjup,Rsun,
     .  Lsun,incl,R2,eccn,tPi,K,AU,b,rdrs,adrs,tdepth,time(1),
     .  exptime(1),dtype(1),tmodel(1)
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
      
      npt=1
      time(1)=solg(7)
      exptime(1)=1765.5/86400.0d0
      dtype(1)=0
      call transitmodel(npt,time,exptime,dtype,tmodel,nfit,solg)
      tdepth=(1.0d0-10**((tmodel(1)-solg(8))/-2.5d0))*1.0d6
      
      write(6,*) int(tdepth+0.5)
 500  format(F6.4)
      
      
      goto 999
 901  write(0,*) "Cannot open" ,inputsol
      goto 999
 999  end
      