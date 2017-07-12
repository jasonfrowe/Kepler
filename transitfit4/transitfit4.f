CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program transitfit4
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Fits for stellar density (a/R*)
      implicit none
      integer iargc,nfit,npt,nmax,i,nunit,nptv,npta,ntt,ndt
      parameter(nfit=18,nmax=2000000)
      integer ia(nfit),dtype(nmax)
      double precision sol(nfit),time(nmax),dt,tmodel(nmax),
     .  flux(nmax),ferr(nmax),itime(nmax),Keplertime,serr(nfit,2),
     .  err(nfit),covar(nfit,nfit),alpha(nfit,nfit),dchistop,
     .  vtime(nmax),vel(nmax),verr(nmax),vetime(nmax),aT(nmax),aM(nmax),
     .  aE(nmax),aIT(nmax),kmag,kerr,tobs(nmax),omc(nmax),ttcor(nmax),
     .  phase,chi
      character*3 titles(nfit)
      character*80 inputsol,obsfile,rvfile,cline,ttfile
      common /Fitting/ npta,aIT
      
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
      
      if(iargc().ge.4)then
        call getarg(4,ttfile)
        if(ttfile.eq.'null')then
            ntt=0
        else
            nunit=10
            open(unit=nunit,file=ttfile,status='old',err=905)
            call readttfile(nunit,nmax,ntt,tobs,omc)
            close(nunit)
        endif
      endif
      
c      open(unit=21,file="junk.dat")
C     If we have TTV info then use it
      if(ntt.gt.0)then
        do 21 i=1,npt
            call lininterp(tobs,omc,ntt,time(i),ttcor(i))
            time(i)=time(i)-ttcor(i)
c            write(21,*) time(i),ttcor,flux(i)
 21     continue
      endif
c      close(21)
c      goto 999
    
C     Store all the observations in the master file
      do 17 i=1,npt !central database of all data
        dtype(i)=0 !0 marks that we have photometric data
        aT(i)=time(i)
        aM(i)=flux(i)
        aE(i)=ferr(i)
        aIT(i)=itime(i)
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

C     Uncomment to skip fitting and just output current fit
c      goto 20

      dchistop=1.0d-10
      call fittransitmodel(npta,aT,aM,aE,dtype,nfit,sol,serr,ia,covar,
     .  alpha,dchistop,err)
     
      call exportfit(nfit,sol,serr,err,titles)

 20   call transitmodel(nfit,sol,npta,aT,aIT,tmodel,dtype)
 
      chi=0.0d0
      do 22 i=1,npt
        chi=chi+(aM(i)-tmodel(i))*(aM(i)-tmodel(i))/(aE(i)*aE(i))
 22   continue
 
      write(0,*) "Chi-squared :",chi
 
 
      ndt=0
      do 10 i=1,npta
c        if(dtype(i).eq.0) write(6,*) aT(i),aM(i),tmodel(i)
        if(dtype(i).eq.ndt)then !select RV or phot
            if(ndt.eq.1) write(6,*) aT(i)+4900.0d0,aM(i)-tmodel(i)
            if(ndt.eq.0) write(6,*) aT(i)+54900.0d0-0.5d0,
     .          aM(i)-1.0d0,tmodel(i)
        endif
c        phase=aT(i)/sol(2)-int(aT(i)/sol(2))
c        if((phase.lt.0.316).or.(phase.gt.0.365))then
c            if(dtype(i).eq.0) write(6,*) 54900.0d0+aT(i)-0.50d0,
c     .          aM(i)-1.0d0,aE(i)
c        endif
c        write(6,*) aT(i),aM(i)-tmodel(i),aE(i)
 10   continue
      
      goto 999
 901  write(0,*) "Usage: transitfit4 <photfile> <rvfile> <fitpars> <ttfi
     .le> <kmag>"
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
 
