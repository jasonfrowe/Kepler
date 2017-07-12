      program transitboot
      implicit none
      integer nmax,iargc,nunit,npt,nfit,npars,itmax,it,i,j,iboot,nprob,
     .  nunitp
      parameter(nmax=600000,nfit=18,iboot=10000)
      integer ipars(nfit),ia(nfit)
      double precision time(nmax),mag(nmax),merr(nmax),tbin,
     .  etime(nmax),MOSTtime,sol(nfit),serr(nfit,2),Dpvary(nfit),toff,
     .  tmodel(nmax),p(nfit+1,nfit),pvary(nfit),doe,serr2(nfit,2),
     .  dchistop,kmag,kerr,err(nfit),inclmin,mass(nmax),massprob(nmax),
     .  rad(nmax),radprob(nmax)
C     RV Parameters
      integer nptv
      double precision vtime(nmax),vel(nmax),verr(nmax),vetime(nmax)
C     Master data array
      integer npta,dtype(nmax)
      double precision aT(nmax),aM(nmax),aE(nmax),aIT(nmax)
c     MRQMIN varibles
      double precision alpha(nfit,nfit),covar(nfit,nfit),chisquared,
     .  rchi
C     Bootstrap parameters
      integer now(3),seed
      double precision sol2(nfit),Dpvary2(nfit),dumr,ran2,
     .  dlk,lk,lkold,likelihood,tmodel2(nmax)
      character dumc
      character*80 obsfile,inputsol,chseed,rvfile,mrprobfile
      logical loop
      common /Fitting/ ipars,npta,aIT,sol2,pvary,inclmin

      
      if(iargc().lt.4) goto 901 !check number of command line arguements

C     Parse the name of the observations data file from the commandline
      call getarg(1,obsfile)
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,mag,merr,etime,MOSTtime)
      call readkeplc(nunit,nmax,npt,time,mag,merr,etime,MOSTtime)
      close(nunit)!release unit number as we are done with the file.

c      kmag=0.0d0
c      if(kmag.gt.0.0)then !estimate of expected scatter
cC       quadratic version
cc        kerr=(4.0d0*(kmag-8.0d0)**2.0+10.0d0)*1.0d-6
cC       Powerlaw version
c         kerr=(10.0d0**(0.23*kmag-0.9))*1.0d-6
c      else
c        kerr=100.0*1.0d-6
c      endif 
c      kerr=2.5*log10(1.0d0+kerr) !convert flux error to mag
cc      kerr=1.6d-4
c      write(0,*) "KMAG,KERR:",kmag,kerr
      
      do 11 i=1,npt
        merr(i)=kerr
 11   continue


      call getarg(2,rvfile) !get filename for RV data
      if(rvfile.eq.'null')then
        nptv=0
      else
        write(0,*) "Reading RV data"
        nunit=10
        open(unit=nunit,file=rvfile,status='old',err=904)
        call readrv(nunit,nmax,nptv,vtime,vel,verr,vetime)
        close(nunit)
c        write(6,*) "rv:",(vtime(i),i=1,nptv)
      endif

      write(0,*) "Reading input solution"
      call getarg(3,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
      call getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff)
      close(nunit) !release unit number as we are done with file
      
      write(0,*) "Reading in mass/radius probabilities"
      call getarg(4,mrprobfile) !get filename for probablities
      nunitp=11
      open(unit=nunitp,file=mrprobfile,status='old',err=905)
      read(nunitp,*) dumc
C     If reading in probs..
c      call getprobs(nunit,nmax,nprob,mass,massprob,rad,radprob)
c      close(nunitp)
c      write(0,*) "nprob:",nprob

C     Store all the observations in the master file
      do 17 i=1,npt !central database of all data
        dtype(i)=0 !0 marks that we have photometric data
        aT(i)=time(i)
        aM(i)=mag(i)
        aE(i)=merr(i)
        aIT(i)=etime(i)
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
      if(iargc().ge.5)then
        call getarg(5,chseed)
        read(chseed,*) seed
      else
        call itime(now)
        seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
      endif
      write(0,*) "Seed: ",seed
      dumr=ran2(-seed)
      
      call transitmodel(npta,aT,aIT,dtype,tmodel2,nfit,sol)
      
      write(0,*) "Start Bootstrapping"
      do 10 j=1,iboot
C       Create bootstrapped data      
        call bootcopy(npt,time,mag,merr,etime,
     .      nptv,vtime,vel,verr,vetime,npta,aT,aM,aE,aIT,
     .      seed,nfit,sol,Dpvary,sol2,Dpvary2,tmodel2,serr,serr2,
     .      nprob,mass,massprob,rad,radprob,nunitp)
      
C       We can bin the data
        tbin=10.0 !*** move this to pars file***
c        call bindt(npt2,time2,mag2,merr2,etime2,tbin)
      
C       Now we auto-explore for the scale length for minimization
        itmax=1 !maximum number of iterations to maximize likelihood
        it=0 !iteration counter
        lkold=likelihood(npta,tmodel,aT,aM,aE,aIT,dtype,nfit,sol2,serr2)
        loop=.true.
        dchistop=1.0e-3 !fit criteria 
        do while(loop) 
            it=it+1
c            write(0,*) "Iteration #:",it  
            call fittransitmodel2(npta,aT,aM,aE,dtype,nfit,sol2,serr2,
     .          npars,ipars,p,pvary,Dpvary2,dchistop,alpha,covar,ia,
     .          inclmin,err)
            rchi=chisquared(npta,tmodel,aT,aM,aE,
     .          aIT,dtype,nfit,sol2)/dble(npt-1)
c            write(0,*) "RChi:",rchi
           lk=likelihood(npta,tmodel,aT,aM,aE,aIT,dtype,nfit,sol2,serr2)
            dlk=(lkold-lk)/lkold
c           write(0,*) "lk: ",lk,lkold,dlk
c           lkold=lk
c           if((abs(dlk).lt.1.0e-8).and.(abs(dlk).gt.0.0d0)) loop=.false.
c           call exportfit(nfit,sol,serr,Dpvary,toff,titles)
            if(it.gt.itmax)loop=.false.
        enddo
c        if(serr(6,2).ne.0)

C        next three lines are for probability distributions
c        if((sol2(1).ge.mass(1)).and.(sol2(1).le.mass(nprob)).and.
c     .      (sol2(3).ge.rad(1)).and.(sol2(3).le.rad(nprob)))
c     .   write(6,500) (sol2(i),i=1,nfit),rchi !write bootstrap to STDOUT
        write(6,500) (sol2(i),i=1,nfit),rchi !write bootstrap to STDOUT
 500    format(19(1X,1PE17.10))      
 10   continue
    
      goto 999
 901  write(0,*) "Usage: transitboot2.0 <obsfile> <rvfile> <inputsol> <p
     .robfile>"
      write(0,*) "<rvfile>=null if no RV available"
      goto 999
 902  write(0,*) "Cannot open ",inputsol
      goto 999
 903  write(0,*) "Cannot open ",obsfile
      goto 999
 904  write(0,*) "Cannot open ",rvfile
      goto 999
 905  write(0,*) "Cannot open ",mrprobfile
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bootcopy(npt,time,mag,merr,etime,
     .      nptv,vtime,vel,verr,vetime,npta,aT,aM,aE,aIT,
     .      seed,nfit,sol,Dpvary,sol2,Dpvary2,tmodel,serr,serr2,
     .      nprob,mass,massprob,rad,radprob,nunitp)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nptv,npta,seed,i,j,nfit,nprob,npm,nunitp
      double precision time(npt),mag(npt),merr(npt),etime(npt),
     .  vtime(nptv),vel(nptv),verr(nptv),vetime(nptv),
     .  aT(npta),aM(npta),aE(npta),aIT(npta),serr2(nfit,2),serr(nfit,2),
     .  ran2,sol(nfit),Dpvary(nfit),sol2(nfit),Dpvary2(nfit),
     .  tmodel(npta),mass(nprob),massprob(nprob),rad(nprob),gasdev,
     .  radprob(nprob),ranprobm,massshuf,radshuf,rmass,rrad,ranprobr
     
C     ***If using probs file***
c      massshuf=(mass(nprob)-mass(1))/dble(nprob-1)
c      radshuf=(mass(nprob)-mass(1))/dble(nprob-1)
c      ranprobm=ran2(seed)
c      ranprobr=ran2(seed)
c      do 13 i=2,nprob
c        if((ranprobm.ge.massprob(i-1)).and.
c     .      (ranprobm.lt.massprob(i)))then
c                rmass=mass(i-1)+ran2(seed)*massshuf
c        endif
c        if((ranprobr.ge.radprob(i-1)).and.
c     .      (ranprobr.lt.radprob(i)))then
c                rrad=rad(i-1)+ran2(seed)*radshuf
c        endif
c 13   continue
c      if(ranprobm.ge.1.0d0) rmass=mass(nprob-1)+ran2(seed)*massshuf
c      if(ranprobr.ge.1.0d0) rrad =rad (nprob-1)+ran2(seed)*radshuf
      
      read(nunitp,*,end=901) rmass,rrad


C     resample data with replacement
      do 10 i=1,npt
        j=ran2(seed)*real(npt)+1.0
        aT(i) =time(i)
        aM(i) =mag(j)-tmodel(j)+tmodel(i)
        aE(i) =merr(j)
        aIT(i)=etime(j)
 10   continue
      do 11 i=1,nptv
        j=ran2(seed)*real(nptv)+1.0
        aT(npt+i) =vtime(i)
        aM(npt+i) =vel(j)-tmodel(npt+j)+tmodel(npt+i)
        aE(npt+i) =verr(j)
        aIT(npt+i)=vetime(j)
 11   continue
      npta=npt+nptv
      
C     Make a copy of the solution parameters      
      do 12 i=1,nfit
        sol2(i)=sol(i)
        Dpvary2(i)=Dpvary(i)
        serr2(i,1)=serr(i,1)
        serr2(i,2)=serr(i,2)
 12   continue
 
 13   if(serr(18,2).gt.0.0d0)sol2(18)=gasdev(seed)*serr(18,2)+serr(18,1)
      if(sol2(18).lt.0.0d0) goto 13
c      write(6,*) "DIL:",sol2(18)
c      read(5,*)
 
C     tag on new fixed mass and rad
      sol2(1)=rmass
      sol2(3)=rrad
      serr2(1,2)=0.0d0
      serr2(3,2)=0.0d0 
      
      goto 999
 901  write(0,*) "End of MCMC data"
      close(nunitp)
      stop
 999  return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getprobs(nunit,nmax,nprob,mass,massprob,rad,radprob)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nmax,nprob,i
      double precision mass(nmax),massprob(nmax),rad(nmax),
     .  radprob(nmax),mnorm,rnorm

C     Read in probability distribution      
      i=1
 10   read(nunit,*,end=11) mass(i),massprob(i),rad(i),radprob(i)
        i=i+1
        goto 10
 11   continue
      nprob=i-1
      
C     sum up all probabilities so we can normalize sum==1
      mnorm=0.0
      rnorm=0.0
      do 12 i=1,nprob
        mnorm=mnorm+massprob(i)
        rnorm=rnorm+radprob(i)
        if(i.gt.1) massprob(i)=massprob(i)+massprob(i-1)
        if(i.gt.1) radprob(i)=radprob(i)+radprob(i-1)
 12   continue
      
C     normalization
      do 13 i=1,nprob
        massprob(i)=massprob(i)/mnorm
        radprob(i)=radprob(i)/rnorm
 13   continue
 
      return
      end