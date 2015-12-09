CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program transitfit5
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Fits for mean stellar density (p_star) for multi-planet candidates
C     gfortran for OS/X fails with minpack as -mcmodel=large is broken
C     for the default assemler. This version falls back to the slower
C     Numerical Recipes version for compatiblity.
      implicit none
      integer iargc,nfit,npt,nmax,i,nunit,nptv,npta,nplanet,ndt,
     .  nplanetmax
      parameter(nfit=108,nmax=2000000,nplanetmax=10)
      integer ia(nfit),dtype(nmax),ntt(nplanetmax)
      double precision sol(nfit),time(nmax),dt,tmodel(nmax),
     .  flux(nmax),ferr(nmax),itime(nmax),Keplertime,serr(nfit,2),
     .  err(nfit),covar(nfit,nfit),alpha(nfit,nfit),dchistop,
     .  vtime(nmax),vel(nmax),verr(nmax),vetime(nmax),aT(nmax),aM(nmax),
     .  aE(nmax),aIT(nmax),kmag,kerr,chisq
C     minpack
      integer iwa(nfit),lwa
      parameter(lwa=nmax*nfit+5*nmax*nfit)
      double precision fvec(nmax),wa(lwa),sol2(nfit)
C     rhostar chi-square?
      integer nfrho
      double precision rhoierr(9),rhoi,chifac
C     TT variations
      double precision tobs(nplanetmax,nmax),omc(nplanetmax,nmax)
      character*3 titles(15)
      character*80 inputsol,obsfile,rvfile,cline,ttfile
#ifdef gfortran
      common /Fitting/ npta,nplanet,aIT,ntt,tobs,omc
#else
      common /Fitting2/ nfrho,nplanet,dtype,aT,aM,aE,aIT,sol,serr,
     .  rhoi,rhoierr,chifac,ntt,tobs,omc
#endif
      
      nfrho=1
      
      if(iargc().lt.3) goto 901

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
c      call readdata(nunit,nmax,npt,time,flux,ferr,itime,Keplertime)
      call readkeplc(nunit,nmax,npt,time,flux,ferr,itime,Keplertime)
      close(nunit)!release unit number as we are done with the file.
C     apply a boxcar filter
c      boxbin=1.0 !filter width (days)
c      call boxcar(npt,time,mag,merr,boxbin)

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
        do 11 i=1,npt
            ferr(i)=kerr
 11     continue
c      else
c        kerr=100.0*1.0d-6
      endif 
      write(0,*) "KMAG,KERR:",kmag,kerr

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
      write(0,*) "reading in input solution"
      call getfitpars(nunit,nfit,nplanet,sol,serr,err)
      write(0,*) "done reading input solution"
      close(nunit) !release unit number as we are done with file
    
C     Store all the observations in the master file
      do 17 i=1,npt !central database of all data
        dtype(i)=0 !0 marks that we have photometric data
        aT(i)=time(i)
        aM(i)=flux(i)
        aE(i)=ferr(i)
        aIT(i)=itime(i)
c        write(0,*) "itime",itime(i)
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

      do 21 i=1,nplanet
        if(iargc().ge.4+i)then
            call getarg(4+i,ttfile)

            if(ttfile.eq.'null')then
                ntt(i)=0
            else
                nunit=10
                open(unit=nunit,file=ttfile,status='old',err=905)
                call readttfile(nunit,nplanetmax,nmax,i,ntt,tobs,omc)
                close(nunit)
            endif
            
        else
            ntt(i)=0
        endif
c        write(0,*) "ntt",i,ntt(i)
 21   continue

C     Skip fitting and just exit..
C      goto 20

C     Use common block Fitting or Fitting2 as required
      dchistop=1.0d-5
#ifdef gfortran
      call fittransitmodel(npta,aT,aM,aE,dtype,nfit,nplanet,sol,serr,ia,
     .  covar,alpha,dchistop,err)
#else
      call fittransitmodel2(npta,nfit,sol,sol2,serr,fvec,iwa,lwa,wa)
#endif

      call exportfit(nfit,nplanet,sol,serr,err,titles)

 20   continue 

      call transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npta,aT,aIT,
     .  ntt,tobs,omc,tmodel,dtype)

      chisq=0.0d0
      do 22 i=1,npta
         chisq=chisq+((aM(i)-tmodel(i))/aE(i))**2.0d0
 22   continue
      write(0,*) "Chisq: ",chisq

 
      ndt=0
      do 10 i=1,npta
        if(dtype(i).eq.ndt)then !select RV or phot
            if(ndt.eq.1) write(6,*) aT(i)+4900.0d0,aM(i)-tmodel(i),aE(i)
            if(ndt.eq.0) then
                write(6,*) aT(i)+54900.0d0-0.5d0,
     .              aM(i),tmodel(i)
c                write(6,*) aT(i)+54900.0d0-0.5d0,aM(i),tmodel(i)
c               write(6,*) aT(i),aM(i),tmodel(i)
            endif
        endif
 10   continue
      
      goto 999
 901  write(0,*) "Usage: transitfit5 <photfile> <rvfile> <fitpars> [kmag
     .] [ttfiles]"
      goto 999
 902  write(0,*) "Error opening ",inputsol
      goto 999
 903  write(0,*) "Error opening ",obsfile
      goto 999
 904  write(0,*) "Error opening ",rvfile
      goto 999
 905  write(0,*) "Error opening ",ttfile
      goto 999
 999  end
