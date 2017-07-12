CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program transitmcmc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Fits for stellar density (a/R*)
      implicit none
      integer iargc,nfit,npt,nmax,i,nunit,nptv,npta,j,naccept,ng
      parameter(nfit=18,nmax=2000000)
      integer dtype(nmax),niter,flag
      double precision sol(nfit),time(nmax),dt,tmodel(nmax),
     .  flux(nmax),ferr(nmax),exptime(nmax),Keplertime,serr(nfit,2),
     .  err(nfit),vtime(nmax),vel(nmax),verr(nmax),vetime(nmax),
     .  aT(nmax),aM(nmax),aE(nmax),aIT(nmax),kmag,kerr,ran2,dumr,
     .  sol2(nfit),rchi,bchi,gasdev,accrate,dil(2)
C     Random number generation
      integer now(3),seed
      character*80 inputsol,obsfile,rvfile,cline,chseed
      
      if(iargc().lt.4) goto 901

c      call defaultpars(nfit,sol)

cC     Make up some fake time stamps
c      npt=2000
c      dt=30.0/1440.0d0 !30 minute cadence 
c      do 10 i=1,npt
c        time(i)=dt*dble(i)
c 10   continue

C     Parse the name of the observations data file from the commandline
      call getarg(1,obsfile)
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,flux,ferr,exptime,Keplertime)
      call readkeplc(nunit,nmax,npt,time,flux,ferr,exptime,Keplertime)
      close(nunit)!release unit number as we are done with the file.
C     apply a boxcar filter
c      boxbin=1.0 !filter width (days)
c      call boxcar(npt,time,mag,merr,boxbin)

      kmag=0.0d0
      if(iargc().ge.5)then
        call getarg(5,cline)
        read(cline,*) kmag
      endif
      if(kmag.gt.0.0)then !estimate of expected scatter
C       quadratic version
c        kerr=(4.0d0*(kmag-8.0d0)**2.0+10.0d0)*1.0d-6
C       Powerlaw version
         kerr=(10.0d0**(0.23*kmag-0.9))*1.0d-6
      else
        kerr=100.0*1.0d-6
      endif 
      write(0,*) "KMAG,KERR:",kmag,kerr
      
      do 11 i=1,npt
        ferr(i)=kerr
 11   continue

      call getarg(2,rvfile) !get filename for RV data
      if(rvfile.eq.'null')then
        nptv=0
      else
        nunit=10
        open(unit=nunit,file=rvfile,status='old',err=904)
        call readrv(nunit,nmax,nptv,vtime,vel,verr,vetime)
        close(nunit)
c        write(6,*) "rv:",(vtime(i),i=1,nptv)
      endif


      call getarg(3,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
      call getfitpars(nunit,nfit,sol,serr,err)
      close(nunit) !release unit number as we are done with file
    
C     Store all the observations in the master file
      do 17 i=1,npt !central database of all data
        dtype(i)=0 !0 marks that we have photometric data
        aT(i)=time(i)
        aM(i)=flux(i)
        aE(i)=ferr(i)
        aIT(i)=exptime(i)
 17   continue
      npta=npt
      do 18 i=1,nptv
        dtype(npt+i)=1 !mark data as RV
        aT(npt+i)=vtime(i)
        aM(npt+i)=vel(i)
        aE(npt+i)=verr(i)
        aIT(npt+i)=vetime(i)
 18   continue
      npta=npta+nptv
      
      write(0,*) "Initialization of random number"
C     Random number initialization(so we get a different seed each time)
      if(iargc().ge.6)then
        call getarg(6,chseed)
        read(chseed,*) seed
      else
        call itime(now)
        seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
      endif
      write(0,*) "Seed: ",seed
      dumr=ran2(-seed)
      
      call transitmodel(nfit,sol,npta,aT,aIT,tmodel,dtype)
      bchi=0.0d0
      do 12 i=1,npta
        bchi=bchi+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i))
 12   continue
      write(0,*) "bchi",bchi,npta-1
      bchi=dble(npta-1)/bchi
      write(0,*) "bchi_f:",bchi
c      bchi=1.0d0

C     Save dilution parameters
      dil(1)=serr(18,1)
      dil(2)=serr(18,2)

C     Initialization of Markov chain
      do 13 i=1,nfit-1
c        sol(i)=serr(i,1)+serr(i,2)*(2.0d0*ran2(seed)-1.0d0)
 22     sol(i)=gasdev(seed)*serr(i,2)+sol(i)
        if((i.eq.3).and.(sol(3).lt.0.0d0)) goto 22
        if((i.eq.3).and.(sol(3).gt.1.0d0)) goto 22
        if((i.eq.4).and.(sol(4).lt.0.0d0)) goto 22
 13   continue
 21   sol(18)=gasdev(seed)*dil(2)+dil(1)
      if(sol(18).lt.0.0d0) goto 21

      naccept=0 !calculate accept rate
      call getarg(4,cline)
      read(cline,*) niter
      do 19 j=1,niter
        call mcmc(npta,aT,aM,aE,aIT,dtype,nfit,sol,sol2,serr,tmodel,r
     .      chi,seed,flag,bchi,ng,dil)
        if(flag.eq.0)then
            do 20 i=1,nfit
                sol(i)=sol2(i)
 20         continue
            naccept=naccept+1
        endif
        if(j.gt.1) write(6,500) (sol(i),i=1,nfit),rchi,flag,ng
 500    format(19(1X,1PE17.10),2(1X,I2)) 
 19   continue
      accrate=dble(naccept)/dble(niter)
      write(0,*) "Accept Rate: ",accrate
      
      goto 999
 901  write(0,*) "Usage: transitmcmc <photfile> <rvfile> <fitpars> <nite
     .r>"
      goto 999
 902  write(0,*) "Error opening ",inputsol
      goto 999
 903  write(0,*) "Error opening ",obsfile
      goto 999
 904  write(0,*) "Error opening ",rvfile
      goto 999
 999  end
 
