C     New and improved photometry transit fitting program
C     All calculations done in <b>double precision!<\b>
      program transitfit
C     (c) Jason Rowe 2009 (jasonfrowe@gmail.com)
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
      integer nfit,iargc,nunit,nmax,npt,i,itmax,it
C     There are 13 parameters that define the lightcurve shape
      parameter(nfit=17,nmax=600000)
C     px,py,pz: PGPLOT only accepts real data, so these are temp arrays
      real px(nmax),py(nmax),pz(nmax),dtype(nmax)
C     time:  time at middle of exposure [d]
C     mag:   instrumental magnitude [mag]
C     merr:  observational error [mag]  
C     itime: integration time [d]
C     MOSTtime: HJD offset
      double precision time(nmax),mag(nmax),merr(nmax),itime(nmax),
     .  MOSTtime,phase(nmax),tbin
C     sol: contains the solution
C     serr: contains the errors on the solution (see getfitpars)
C     toff: optional parameter to move transit to phase 0.75
C     bb: time-series plotting bounds
      double precision sol(nfit),serr(nfit,2),toff,bb(4),bbp(4),bbp2(4),
     .  avgitime,mean,chisquared,likelihood,lk,lkold,dlk,dchistop,
     .  err(nfit),doe
C     ipars: labels for fitted parameters passed in fittransitmodel
C     p: array fed to amoeba to search for local minimum
C     pvary: used to construct p, to search parameter space
C     Dpvary: how much to vary solution when no priors are available
      integer ipars(nfit),tp(nfit+1)
      double precision p(nfit+1,nfit),pvary(nfit),y(nfit+1),x(nfit),
     .  Dpvary(nfit)
C     Variable for cfunk
      integer npars
      double precision sol2(nfit),tmodel(nmax),soltemp(nfit)
C     inputsol: filename for input model parameters
C     obsfile:  filename containing observations to model
      character*80 inputsol,obsfile
C     for exporting the new solution
      character*3 titles(nfit)
      logical loop
      common /Fitting/ npars,ipars,npt,time,mag,merr,itime,sol,serr,
     .  sol2,tmodel
c      data Dpvary /0.2d0,0.2d0,0.2d0,0.2d0,1.0d-5,1.0,0.05,1.0d-5,0.5d0,
c     .  0.01d0,0.01d0,0.01d0,0.01d0/
      
      if(iargc().lt.2) goto 901 !check number of command line arguements

C     Parse the name of the observations data file from the commandline
      call getarg(1,obsfile)
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,mag,merr,itime,MOSTtime)
      call readkeplc(nunit,nmax,npt,time,mag,merr,itime,MOSTtime)
      close(nunit)!release unit number as we are done with the file.

      call getarg(2,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
      call getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff)
      close(nunit) !release unit number as we are done with file
      
C     We can bin the data
      tbin=30.0
c      call bindt(npt,time,mag,merr,itime,tbin)

      call opengraphics()
     
      dchistop=1.0e-8 !criteria for good fit.
      itmax=100
      it=0
      lkold=likelihood(npt,tmodel,time,mag,merr,itime,nfit,sol,serr)
      loop=.true.
      do while(loop) 
        it=it+1
        write(0,*) "Iteration #:",it
        call autoexplore(nmax,npt,time,mag,merr,itime,tmodel,nfit,
     .      sol,serr,Dpvary,soltemp)      
        call fittransitmodel(nfit,sol,serr,npars,ipars,p,pvary,Dpvary,
     .      x,y,tp,dchistop)
        write(0,*) "RChi:",chisquared(npt,tmodel,time,mag,merr,
     .      itime,nfit,sol)/dble(npt-1)
        lk=likelihood(npt,tmodel,time,mag,merr,itime,nfit,sol,serr)
        dlk=(lkold-lk)/lkold
        write(0,*) "lk: ",lk,lkold,dlk
        lkold=lk
        if((abs(dlk).lt.1.0e-8).and.(abs(dlk).gt.0.0d0)) loop=.false.
        call exportfit(nfit,sol,serr,Dpvary,err,doe,toff,titles)
 10   enddo
     
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

C     Plotting the solution     
      avgitime=mean(npt,itime)
      call modelplot(nfit,sol,nmax,phase,bb,bbp,bbp2,avgitime,toff)      
      
      call closegraphics()
      
      goto 999
 901  write(0,*) "Usage: transitfit datafile.dat fitpars.dat"
      write(0,*) "  datafile.dat : contains MOST photometry"
      write(0,*) "  fitspars.dat : contains fitting information"
      goto 999
 902  write(0,*) "Cannot open ",inputsol
      goto 999
 903  write(0,*) "Cannot open ",obsfile
 999  end
 

      