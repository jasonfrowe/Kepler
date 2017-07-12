C     New and improved photometry transit fitting program
C     All calculations done in <b>double precision!<\b>
      program transitfit
C     (c) Jason Rowe 2011 (jasonfrowe@gmail.com)
C     Version 3.0 (July 07/2009)
C     +moved from amoeba to mrqmin
C     Version 2.1.0 (April 03/2009)
C     +added non-circular orbits
C     +an autoexplore to define search scale length for amoeba
C     Version 2.0.2 (March 26/2009)
C     +fixed segmentation fault in mandelagol
C     +fitting appear to work for chi and bayesian
C     +reworked ploting surface for better panel shapes 
C     -fitting only converges for binned data
C     -no convergence for reflected light
C     +use of support program modelbuild really helps
C     +switched input variables to match modelbuild logic
C     Version 2.0.1 (March 24/2009)
C     +reads in MOST data as real*8
C     +toff is now a phase-offset (easier to use -maybe call poff now?)
C     +can overlay model plots now
C     Version 2.0.0 (March 23/2009)
C     +model parameters now read from external file
C     +Transit model is now 100% real*8
C     +Integration time accounted for in transit model
      implicit none
C     nfit: defines the number of parameters fit by the model.
C     nmodel: defines which model solution to start with
C     iargc: gives the number of input arguements
C     nunit: unit number assigned for reading in data
C     nmax: maximum number of data points that can be read in
C     npt: number of data points read in
      integer nfit,iargc,nunit,nmax,npt,itmax,it,i,ip,iper,npta,ndt
C     There are 13 parameters that define the lightcurve shape
      parameter(nfit=18,nmax=600000)
C     px,py,pz: PGPLOT only accepts real data, so these are temp arrays
      real px(nmax),py(nmax),pz(nmax)
C     time:  time at middle of exposure [d]
C     mag:   instrumental magnitude [mag]
C     merr:  observational error [mag]  
C     itime: integration time [d]
C     MOSTtime: HJD offset
      double precision time(nmax),mag(nmax),merr(nmax),itime(nmax),
     .  MOSTtime,phase(nmax),tbin,boxbin,kmag,kerr,ptry(4),pold,pchi,
     .  epoch,time2(nmax),mag2(nmax),merr2(nmax),doe,ph
C     sol: contains the solution
C     serr: contains the errors on the solution (see getfitpars)
C     toff: optional parameter to move transit to phase 0.75
C     bb: time-series plotting bounds
      double precision sol(nfit),serr(nfit,2),toff,bb(4),bbp(4),bbp2(4),
     .  avgitime,mean,chisquared,likelihood,lk,lkold,dlk,dchistop,
     .  err(nfit),inclmin,rchi,lchi,pi,tpi,sol2(nfit),serr2(nfit,2),
     .  err2(nfit),bbp3(4),avgvtime,jrv
C     RV parameters
      integer nptv
      double precision vtime(nmax),vel(nmax),verr(nmax),vetime(nmax)
C     Data holders (holds RV+Photometric data)
      integer dtype(nmax),dtype2(nmax)
      double precision aT(nmax),aM(nmax),aE(nmax),aIT(nmax)
C     ipars: labels for fitted parameters passed in fittransitmodel
C     p: array fed to mrqmin to search for local minimum
C     pvary: used to construct p, to search parameter space
C     Dpvary: how much to vary solution when no priors are available
C     Alpha,Covar,ia: Used by mrqmin
      integer ipars(nfit),ia(nfit)
      double precision p(nfit),pvary(nfit),
     .  Dpvary(nfit),alpha(nfit,nfit),covar(nfit,nfit)
C     Variable for funk
      integer npars
      double precision tmodel(nmax),soltemp(nfit),temp
C     inputsol: filename for input model parameters
C     obsfile:  filename containing observations to model
      character*80 inputsol,obsfile,cline,rvfile
C     for exporting the new solution
      character*3 titles(nfit)
      logical loop
      common /Fitting/ ipars,npta,aIT,sol,pvary,inclmin
c      data Dpvary /0.2d0,0.2d0,0.2d0,0.2d0,1.0d-5,1.0,0.05,1.0d-5,0.5d0,
c     .  0.01d0,0.01d0,0.01d0,0.01d0/
      
      jrv=2.5 !RV stellar jitter (need to move to f1.dat)
      Pi=acos(-1.d0)!define Pi and 2*Pi
      tPi=2.0d0*Pi   
      
      if(iargc().lt.3) goto 901 !check number of command line arguements

C     Parse the name of the observations data file from the commandline
      call getarg(1,obsfile)
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,mag,merr,itime,MOSTtime)
      call readkeplc(nunit,nmax,npt,time,mag,merr,itime,MOSTtime)
      close(nunit)!release unit number as we are done with the file.
C     apply a boxcar filter
c      boxbin=1.0 !filter width (days)
c      call boxcar(npt,time,mag,merr,boxbin)

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
      call getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff)
      close(nunit) !release unit number as we are done with file
      
      kmag=0.0d0
      if(iargc().ge.4)then
        call getarg(4,cline)
        read(cline,*) kmag
      endif
      if(kmag.gt.0.0)then !estimate of expected scatter
C       quadratic version
c        kerr=(4.0d0*(kmag-8.0d0)**2.0+10.0d0)*1.0d-6
C       Powerlaw version
         kerr=(10.0d0**(0.23*kmag-0.9))*1.0d-6
      else
        kerr=100.0*1.0d-6!*sqrt(30.0d0)
      endif 
      kerr=2.5*log10(1.0d0+kerr) !convert flux error to mag
c      kerr=1.6d-4
      write(0,*) "KMAG,KERR:",kmag,kerr
      
      do 11 i=1,npt
        merr(i)=kerr
 11   continue
   
C     We can bin the data
      tbin=30.0
c      call bindt(npt,time,mag,merr,itime,tbin)

C     Store all the observations in the master file
      do 17 i=1,npt !central database of all data
        dtype(i)=0 !0 marks that we have photometric data
        aT(i)=time(i)
        aM(i)=mag(i)
        aE(i)=merr(i)
        aIT(i)=itime(i)
 17   continue
      npta=npt
      do 18 i=1,nptv
        dtype(npt+i)=1 !mark data as RV
        aT(npt+i)=vtime(i)
        aM(npt+i)=vel(i)
        aE(npt+i)=sqrt(verr(i)*verr(i)+jrv*jrv)
        aIT(npt+i)=vetime(i)
 18   continue
      npta=npta+nptv

c      goto 21 !skip fitting and just output residuals

c      call opengraphics()
      call pgopen('/null')
      
c      goto 16
      
      iper=0
 15   dchistop=1.0d-10 !criteria for good fit.
      iper=iper+1
      itmax=1 !just in case...
      it=0
      lkold=likelihood(npta,tmodel,aT,aM,aE,aIT,dtype,nfit,sol,serr)
c      write(6,*) "lkold: ",lkold
      loop=.true.
      do while(loop) 
        it=it+1
        write(0,*) "Iteration #:",it
        temp=serr(6,2)
        serr(6,2)=0.0
        serr(6,2)=temp      
        call fittransitmodel2(npta,aT,aM,aE,dtype,nfit,sol,serr,npars,
     .      ipars,p,pvary,Dpvary,dchistop,alpha,covar,ia,inclmin,err)
        rchi=chisquared(npta,tmodel,aT,aM,aE,
     .      aIT,dtype,nfit,sol)/dble(npt-1)
        write(0,*) "RChi:",rchi
        lk=likelihood(npta,tmodel,aT,aM,aE,aIT,dtype,nfit,sol,serr)
        dlk=(lkold-lk)/lkold
c        write(0,*) "lk: ",lk,lkold,dlk
c        lkold=lk
c        if((abs(dlk).lt.1.0e-8).and.(abs(dlk).gt.0.0d0)) loop=.false.
c        call exportfit(nfit,sol,serr,Dpvary,toff,titles)
        if(it.gt.itmax)loop=.false.
 10   enddo
c      goto 125 !skip even odd
     
      write(0,*) "Checking period multiples" 
      
      pold=sol(5)
c      epoch=sol(5)*(0.5-sol(7)/tpi)
      ip=0
      do 14 i=1,0
        if(i.eq.1) ptry(i)=pold/3.0d0
        if(i.eq.2) ptry(i)=pold/2.0d0
        if(i.eq.3) ptry(i)=pold*2.0d0
        if(i.eq.4) ptry(i)=pold*3.0d0
        sol(5)=ptry(i)
c        sol(7)=pi-2.0*pi*(epoch/sol(5)-int(epoch/sol(5)))
        pchi=chisquared(npta,tmodel,aT,aM,aE,
     .      aIT,dtype,nfit,sol)/dble(npt-1)
        write(0,*) ptry(i),sol(6),pchi
        if(pchi.lt.rchi)then
            ip=i
            rchi=pchi
        endif
 14   continue
      if(ip.eq.0)then
        sol(5)=pold
      else
        sol(5)=ptry(ip)
      endif
c      sol(7)=pi-2.0*pi*(epoch/sol(5)-int(epoch/sol(5)))
c      if(sol(7).lt.0.0)sol(7)=sol(7)+2.0*pi
      if((ip.ne.0).and.(iper.lt.4)) goto 15
      
      write(0,*) "Computing odd-even.."
      
 16   call oddeven(npta,aT,aM,aE,dtype,nfit,sol,serr,npars,
     .  ipars,p,pvary,Dpvary,dchistop,alpha,covar,ia,inclmin,err,phase,
     .  time2,mag2,merr2,sol2,serr2,doe,err2,dtype2)    
      write(0,*) "Done odd-even"
     
      lchi=0.0
      do 13 i=1,npt
        lchi=lchi+(sol(8)-mag(i))*(sol(8)-mag(i))/(merr(i)*merr(i))
 13   continue
      lchi=lchi/dble(npt)
      
 125  open(unit=12,file="fitstats.dat")
      write(12,*) (sol(i),i=1,nfit),(err(i),i=1,nfit),rchi,lchi
      close(12)
 12   continue

      write(0,*) "Exporting fit"     
      call exportfit(nfit,sol,serr,Dpvary,err,doe,toff,titles)
c      
C     Plot of the observations
      call plotdata(npt,time,mag,merr,-1.0d0,-1.0d0,1,MOSTtime,nmax,px,
     .  py,pz,bb)
c     Phase plot of the observations
      call plotph(npt,time,mag,merr,sol(5),0.0d0,1.0d0,2,-1,nmax,px,
     .  py,pz,phase,toff,1,bbp)
     
      call plotph(npt,time,mag,merr,sol(5),0.700d0,0.800d0,3,-1,nmax,px,
     .  py,pz,phase,toff,0,bbp2)
     
      call plotrvph(nptv,vtime,vel,verr,sol(5),0.0d0,1.0d0,4,17,nmax,
     .  px,py,pz,phase,toff,0,bbp3)

C     Plotting the solution     
      avgitime=mean(npt,itime)
      avgvtime=mean(nptv,vetime)
      call modelplot(nfit,sol,nmax,phase,bb,bbp,bbp2,bbp3,avgitime,
     .  avgvtime,toff)
c      call modelplot(nfit,sol,nmax,phase,bb,bbp,bbp2,avgitime,toff)      
      
      call closegraphics()
      
 21   call transitmodel(npta,aT,aIT,dtype,tmodel,nfit,sol)
c      do 19 i=1,npt
c        write(6,*) time(i)+54900.0d0-0.5d0,
c     .      10.0**((mag(i)-tmodel(i))/-2.5d0)-1.0d0,merr(i)
cc        write(6,*) time(i),mag(i),tmodel(i)!,merr(i)
c 19   continue
      ndt=0
      do 20 i=1,npta
        if(dtype(i).eq.ndt)then !select RV or phot
            if(ndt.eq.1) write(6,*) aT(i)+4900.0d0,aM(i)-tmodel(i),aE(i)
c            if(ndt.eq.0) write(6,*) aT(i)+54900.0d0-0.5d0,
c     .          10.0**((aM(i)-tmodel(i))/-2.5d0)-1.0d0,aE(i)
c            ph=(aT(i)-sol(7))/sol(5)-int((aT(i)-sol(7))/sol(5))
c            if(ph.lt.0.0) ph=ph+1.0
c            if((ph.lt.0.44).or.(ph.gt.0.55))then
c               if(ndt.eq.0) write(6,500) aT(i)+54900.d0-0.5d0,
c     .             10.0**(aM(i)/-2.5d0)-1.0d0,aE(i)
c            endif
c            if(ndt.eq.0) write(6,*) 10.0**((tmodel(i))/-2.5d0)-1.0
c            if(ndt.eq.0) write(6,*) aT(i),aM(i),tmodel(i),aE(i)
C            For Steve B...
            if(ndt.eq.0) write(6,500) aT(i),
     .          10.0**(aM(i)/-2.5d0)-1.0d0,
     .          10.0**(tmodel(i)/-2.5d0)-1.0d0!,aE(i)
        endif
  20  continue
 500  format(19(1X,1PE17.10)) 
      
      goto 999
 901  write(0,*) "Usage: transitfit datafile.dat fitpars.dat"
      write(0,*) "  datafile.dat : contains MOST photometry"
      write(0,*) "  rvdata.dat   : contains RV (enter 'null' if none)"
      write(0,*) "  fitspars.dat : contains fitting information"
      goto 999
 902  write(0,*) "Cannot open ",inputsol
      goto 999
 903  write(0,*) "Cannot open ",obsfile
      goto 999
 904  write(0,*) "Cannot open ",rvfile
      goto 999
 999  end
 

      
