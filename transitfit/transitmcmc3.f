      program transitboot
      implicit none
      integer nmax,iargc,nunit,npt,nfit,i,iboot,nunitp,ndmp,naccept,
     .  niter,j,flag,ng
      parameter(nmax=600000,nfit=18,iboot=10000,ndmp=100000)
      integer ipars(nfit)
      double precision time(nmax),mag(nmax),merr(nmax),accrate,
     .  etime(nmax),MOSTtime,sol(nfit),serr(nfit,2),Dpvary(nfit),toff,
     .  tmodel(nmax),pvary(nfit),doe,rchi,tmp,
     .  kmag,kerr,err(nfit),inclmin,
     .  mdmp(ndmp),rdmp(ndmp),bchi,gasdev,
     .  dil(2)
C     RV Parameters
      integer nptv
      double precision vtime(nmax),vel(nmax),verr(nmax),vetime(nmax)
C     Master data array
      integer npta,dtype(nmax)
      double precision aT(nmax),aM(nmax),aE(nmax),aIT(nmax)
cc     MRQMIN varibles
c      double precision alpha(nfit,nfit),covar(nfit,nfit),chisquared,
c     .  rchi
C     Bootstrap parameters
      integer now(3),seed
      double precision sol2(nfit),dumr,ran2
      character dumc
      character*80 obsfile,inputsol,chseed,rvfile,mrprobfile,cline
      common /Fitting/ ipars,npta,aIT,sol2,pvary,inclmin

c      write(0,*) iargc()
      if(iargc().lt.5) goto 901 !check number of command line arguements

C     Parse the name of the observations data file from the commandline
      call getarg(1,obsfile)
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,mag,merr,etime,MOSTtime)
      call readkeplc(nunit,nmax,npt,time,mag,merr,etime,MOSTtime)
      close(nunit)!release unit number as we are done with the file.

      kmag=0.0d0
      if(iargc().ge.6)then
        call getarg(6,cline)
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
      kerr=2.5*log10(1.0d0+kerr) !convert flux error to mag
c      kerr=1.6d-4
      write(0,*) "KMAG,KERR:",kmag,kerr
      
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
c      read(nunitp,*) dumc
      do 19 i=1,ndmp
c        write(0,*) i
c        read(nunitp,*) mdmp(i),dumr,dumr,rdmp(i)
        read(nunitp,*) mdmp(i),rdmp(i)
 19   continue

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
      if(iargc().ge.7)then
        call getarg(7,chseed)
        read(chseed,*) seed
      else
        call itime(now)
        seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
      endif
      write(0,*) "Seed: ",seed
      dumr=ran2(-seed)
      
      call transitmodel(npta,aT,aIT,dtype,tmodel,nfit,sol)
      bchi=0.0d0
      do 12 i=1,npta
        bchi=bchi+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i))
 12   continue
      write(0,*) "bchi",bchi,npta-1
      bchi=dble(npta-1)/bchi
      write(0,*) "bchi_f:",bchi
      
C     Save dilution parameters
      dil(1)=serr(18,1)
      dil(2)=serr(18,2)

C     Initialization of Markov chain
 
      do 13 i=1,nfit-1
c        sol(i)=serr(i,1)+serr(i,2)*(2.0d0*ran2(seed)-1.0d0)
 14     continue
        tmp=gasdev(seed)*serr(i,2)+sol(i)
        if(i.eq.4)tmp=sol2(3)*(sol(4)/sol(3)+gasdev(seed)*serr(4,2))
c        if(i.eq.6)tmp=sol2(4)/sol2(3)*(sol(6)/(sol(4)/sol(3))+
c     .      gasdev(seed)*serr(6,2))
        if((i.eq.1).and.(tmp.lt.0.0)) goto 14
        if((i.eq.2).and.(tmp.lt.0.0)) goto 14
        if((i.eq.4).and.(tmp.lt.0.0)) goto 14
        if((i.eq.6).and.(tmp.gt.90.0)) goto 14
        sol2(i)=tmp
 13   continue
 21   sol(18)=gasdev(seed)*dil(2)+dil(1)
      if(sol(18).lt.0.0d0) goto 21
      do 15 i=1,nfit-1
        sol(i)=sol2(i)
 15   continue
      
      naccept=0 !calculate accept rate
      call getarg(5,cline)
      read(cline,*) niter
      do 22 j=1,niter
        call mcmc(npta,aT,aM,aE,aIT,dtype,nfit,sol,sol2,serr,tmodel,
     .      rchi,seed,flag,bchi,ng,dil,ndmp,mdmp,rdmp)
        if(flag.eq.0)then
            do 20 i=1,nfit
                sol(i)=sol2(i)
 20         continue
            naccept=naccept+1
        endif
        if(j.gt.1) write(6,500) (sol(i),i=1,nfit),rchi,flag,ng
 500    format(19(1X,1PE17.10),2(1X,I2)) 
 22   continue
      accrate=dble(naccept)/dble(niter)
      write(0,*) "Accept Rate: ",accrate

      close(nunitp)
    
      goto 999
 901  write(0,*) "Usage: transitboot2.0 <obsfile> <rvfile> <inputsol> <p
     .robfile> <niter>"
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

