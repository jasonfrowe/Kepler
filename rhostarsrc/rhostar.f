      program rhostar
C     Matches Teff/log(g)/rhostar/[Fe/H]/L to Y^2 Stellar models
C     Jason Rowe - jasonfrowe@gmail.com
      implicit none
      integer iargc,nunit,i,nline,inda,indz,indm,nmodel,nmodelmax,seed,
     .  now(3),niter,flag,j,naccept,ng,nfrho,nsol,nfitm,nbuffer
      parameter(nmodelmax=1000,nsol=7,nfitm=3,nbuffer=100)
      double precision adrsi,adrsierr(6),massi,radi,Teff,Tefferr,Z,
     .  Zerr(2),
     .  agei,getage,tage(nmodelmax),tTeff(nmodelmax),tlogL(nmodelmax),
     .  dmass,dage,trad(nmodelmax),G,Pi,rhoi,Psec,dumr,ran2,sol(nsol),
     .  gasdev,dsol(3),y2(nmodelmax),rhoierr(9),yprho(9),rhoin(9),
     .  yp1,ypn,sol2(nsol),accrate,rchi,trho(nmodelmax),
     .  tdrhodt(nmodelmax),rhoc,rhocerr,dZ,L,Lerr,logg,loggerr,
     .  tmcore(nmodelmax),aotu,tdloggdt(nmodelmax)
      integer ngcor(nfitm),ngcorsub(nfitm),ngprob(nfitm),
     .  ngprobsub(nfitm),nb,nmov,nup,npars,nupdate,nacor,nacorsub,
     .  naprob,naprobsub
      double precision gscale(nfitm),buffer(nfitm,nbuffer),corscale
      character*80 inputparms,workdir,cline
      parameter (nline=150)
      parameter (inda=3)
      parameter (indz=11)
      parameter (indm=35)
      integer a_one,z_one,m_one
      real*8 avalue(inda)
      real*8 zvalue(indz)
      real*8 xmass(indm)
      common /grid/avalue,zvalue,xmass
      common /single/a_one,z_one,m_one
      data rhoin/-1.0d2,-3.0d0,-2.0d0,-1.0d0,0.0d0,1.0d0,2.0d0,3.0d0,
     .  1.0d2/
      
!     simple check to see if a Y^2 track is present.
      open(unit=10,file="a0o2/x53z08/m04x53z08.track1",
     .  status='old',err=904)
!     if the file opened without an error, close it and move on.
      close(10)
      
      aotu=14.0d0 !Age Of The Universe.
      Pi=acos(-1.d0)!define Pi and 2*Pi
      G=6.674d-11 !N m^2 kg^-2  Gravitation constamt
      
C     Initialize Gibbs scaler
      do 18 i=1,nfitm
        gscale(i)=1.0d0
        ngcor(i)=0
        ngcorsub(i)=0
        ngprob(i)=0
        ngprobsub(i)=0
 18   continue
      corscale=1.0d0
      nacor=0
      nacorsub=0
      naprob=0
      naprobsub=0
      nb=0
      
      if(iargc().lt.2) goto 901
      call getarg(1,inputparms)
      call getarg(2,cline)
      read(cline,*) niter
      nunit=10
      open(unit=nunit,file=inputparms,status='old',err=902)
      call getinputpars(nunit,nfrho,adrsi,adrsierr,massi,radi,agei,Teff,
     .  Tefferr,Z,Zerr,dZ,L,Lerr,Psec,dmass,dage,workdir,rhoc,
     .  rhocerr,logg,loggerr,aotu)
      close(nunit)
c      write(0,*) adrsi,(adrsierr(i),i=1,6)
c      write(0,*) Teff,Tefferr
c      write(0,*) Z,Zerr
c      write(0,*) Psec/(24.0d0*60.0d0*60.0d0)
c      write(0,*) massi,radi
c      write(0,*) dmass,dage
c      write(0,*) workdir
      
c1      write(0,*) "Teff: ",Teff,Tefferr
c1      call gettrack(massi,Z,nmodelmax,nmodel,tage,tTeff,tlogL,trad,trho,
c1     .  tdrhodt,tmcore)      
c      do 10 i=1,nmodel
c        write(6,*) tage(i),tTeff(i),10**tlogL(i),trad(i)
c 10   continue
      rhoi=adrsi**3.0*Pi*3.0d0/(Psec*Psec*G)-rhoc
      write(0,*) "Rho: ",rhoi,adrsi
      write(0,*) "Rhoc:",rhoc,rhocerr
      rhoierr(1)=-rhoi
      do 15 i=2,4
        rhoierr(i)=(adrsi+adrsierr(i-1))**3.0*Pi*3.0d0/(Psec*Psec*G)
     .      -rhoi-rhoc
        rhoierr(i)=sqrt(rhoierr(i)**2+(dble(5-i)*rhocerr)**2)*
     .      abs(rhoierr(i))/rhoierr(i)
        
        write(0,*) rhoierr(i),adrsierr(i-1)
 15   continue
      rhoierr(5)=0.0d0
      do 16 i=6,8
        rhoierr(i)=(adrsi+adrsierr(i-2))**3.0*Pi*3.0d0/(Psec*Psec*G)
     .      -rhoi-rhoc
        rhoierr(i)=sqrt(rhoierr(i)**2+(dble(i-5)*rhocerr)**2)
        write(0,*) rhoierr(i),adrsierr(i-2)
 16   continue
      rhoierr(9)=rhoierr(8)*10.0d0!20000.0-rhoi
      yp1=1.0e30
      ypn=1.0e30
c1      call spline(rhoierr,rhoin,9,yp1,ypn,yprho)
c1      agei=getage(massi,Teff,Tefferr,rhoi,rhoierr,rhoin,yprho,nmodelmax,
c1     .  nmodel,tTeff,tage,trad,flag,aotu)
c1      if(flag.eq.1) goto 903  
c      write(6,*) "Agei: ",Agei

C     We should do a chi-square minimization here!
C     call rhomin() 
      
C     Initialize random number generator      
      call itime(now)
      seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
      dumr=ran2(-seed)
      
 11   sol(1)=gasdev(seed)*dmass+massi
      if(sol(1).lt.0.4) goto 11
      if(sol(1).gt.5.2) goto 11
      dsol(1)=dmass
 12   sol(2)=gasdev(seed)*dage+agei
      if(sol(2).lt.0.0) goto 12
      if(sol(2).gt.aotu) goto 12
      dsol(2)=dage
 13   sol(3)=gasdev(seed)*dZ+Z
      if(sol(3).lt.0.00001) goto 13
      if(sol(3).gt.0.08) goto 13
      dsol(3)=dZ
      
c      write(6,*) (sol(i),i=1,3)
     
C     initializing counters for MCMC
      nmov=0
      nup=0
      npars=1
      nupdate=0
      
      naccept=0 !calculate accept rate
      j=1
      do 14 while (j.lt.niter)!j=1,niter+1
        call mcmc(seed,nsol,sol,sol2,dsol,nmodelmax,tage,tTeff,tlogL,
     .      trad,y2,Teff,Tefferr,Z,Zerr,L,Lerr,logg,loggerr,nfrho,rhoi,
     .      rhoierr,rhoin,yprho,flag,ng,rchi,nfitm,gscale,ngcor,
     .      ngcorsub,ngprob,ngprobsub,nbuffer,nb,nmov,buffer,
     .      nup,npars,nupdate,corscale,nacor,nacorsub,naprob,naprobsub,
     .      tmcore,aotu,tdloggdt)
        if(flag.eq.0)then
            do 17 i=1,nsol
                sol(i)=sol2(i)
 17         continue
            if((j.gt.1).and.(nmov.ge.nbuffer)) naccept=naccept+1
        endif
        if((j.gt.1).and.(nmov.ge.nbuffer)) 
     .      write(6,500) (sol(i),i=1,nsol),rchi,flag,ng
 500    format(8(1X,1PE17.10),2(1X,I2)) 
        if(nmov.ge.nbuffer) j=j+1
 14   enddo
      accrate=dble(naccept)/dble(niter)
      write(0,*) "Accept Rate: ",accrate
      
      goto 999
 901  write(0,*) "Usage: rhostar filename.rho niter"
      write(0,*) " filename.rho contains stellar parameters"
      write(0,*) " niter: number of MCMC iterations"
      goto 999
 902  write(0,*) "Cannot open ",inputparms
      goto 999
 903  write(0,*) "Cannot place age on initial Teff,Mass pair"
      goto 999
 904  write(0,*) "Y^2 Models seems to be missing"
      write(0,*) "a0o2,a2o2,a4o2 directories need to be in working dir"
      goto 999
 999  end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine msmasstemp(Teff,mass)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nms,i
      parameter(nms=61)
      double precision msmass(nms),mstemp(nms),Teff,mass
      data msmass/0.08,0.09,0.10,0.13,0.15,0.20,0.25,0.30,0.35,0.40,
     .  0.45,0.50,0.55,0.60,0.70,0.80,0.90,1.00,1.10,1.20,1.30,1.40,
     .  1.50,1.60,1.70,1.80,1.90,2.00,2.10,2.20,2.30,2.40,2.50,2.60,
     .  2.80,3.00,3.20,3.40,3.60,3.80,4.00,4.20,4.40,4.60,4.80,5.00,
     .  5.25,5.50,5.75,6.00,6.50,7.00,7.50,8.00,8.50,9.00,9.50,10.0,
     .  11.0,12.0,13.0/
      data mstemp/2258.2,2787.4,2874.9,3060.5,3126.3,3271.6,3405.7,
     .  3524.2,3624.8,3723.6,3822.2,3925.2,4046.6,4200.8,4599.3,5012.1,
     .  5362.6,5649.3,5883.0,6076.7,6242.9,6399.4,6554.0,6711.5,6884.0,
     .  7129.2,7392.7,7634.3,7913.4,8188.1,8457.9,8723.0,8984.6,9243.0,
     .  9750.8,10247.4,10733.3,11208.7,11674.1,12128.7,12571.2,13000.3,
     .  13417.0,13824.2,14222.1,14610.2,15084.1,15546.9,15999.2,16441.3,
     .  17296.5,18113.2,18889.5,19621.4,20317.4,20987.7,21633.7,22256.4,
     .  23427.7,24516.3,25533.8/
      
      if(Teff.le.mstemp(1))then
        mass=msmass(1)
      elseif(Teff.gt.mstemp(nms))then
        mass=msmass(nms)
      else
        do 10 i=1,nms-1
            if((Teff.gt.mstemp(i)).and.(Teff.le.mstemp(i+1)))then
                mass=(mstemp(i+1)-Teff)/(mstemp(i+1)-mstemp(i))*
     .              (msmass(i+1)-msmass(i))+msmass(i)
            endif
 10     continue
      endif
      
c      mass=1.4
      
      
      return
      end
