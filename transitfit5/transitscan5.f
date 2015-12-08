CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program transitfitscan5
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Fits for mean stellar density (p_star) for multi-planet candidates
      implicit none
      integer iargc,nfit,npt,nmax,i,nunit,nptv,npta,nplanet,ndt,k,
     .  niter1,niter2,nplanetmax
      parameter(nfit=108,nmax=2000000,nplanetmax=10)
      integer ia(nfit),dtype(nmax),now(3),seed,niter,j,nf,nfrho
      double precision sol(nfit),time(nmax),dt,tmodel(nmax),
     .  flux(nmax),ferr(nmax),exptime(nmax),Keplertime,serr(nfit,2),
     .  err(nfit),dchistop,rhoread(8),rhoi,rhoierr(9),
     .  vtime(nmax),vel(nmax),verr(nmax),vetime(nmax),aT(nmax),aM(nmax),
     .  aE(nmax),aIT(nmax),kmag,kerr,ran2,chisq,dumr,eccn,chifac,
     .  dniter1,dniter2
C     TTV
      integer ntt(nplanetmax)
      double precision tobs(nplanetmax,nmax),omc(nplanetmax)
      character*80 ttfile
C     NR
c      double precision covar(nfit,nfit),alpha(nfit,nfit)
C     minpack
      integer iwa(nfit),lwa
      parameter(lwa=nmax*nfit+5*nmax*nfit)
      double precision fvec(nmax),wa(lwa),sol2(nfit),dniter,sol3(nfit),
     .  serr2(nfit,2),rhoin(9),drho,dsig
      character*3 titles(15)
      character*80 inputsol,obsfile,rvfile,cline,rhofile
c      common /Fitting/ npta,nplanet,aIT
      common /Fitting2/ nfrho,nplanet,dtype,aT,aM,aE,aIT,sol2,serr2,
     .  rhoi,rhoierr,chifac,ntt,tobs,omc
      data rhoin/-1.0d2,-3.0d0,-2.0d0,-1.0d0,0.0d0,1.0d0,2.0d0,3.0d0,
     .  1.0d2/

      do 5 i=1,nfit
        serr(i,1)=0.0d0
        serr(i,2)=0.0d0
        sol2(i)=0.0d0
        sol3(i)=0.0d0
 5    continue

      
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
      if(iargc().ge.6)then
        call getarg(6,cline)
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
      else
        kerr=100.0*1.0d-6
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
c      call getfitpars(nunit,nfit,sol,serr,err)
      close(nunit) !release unit number as we are done with file
      nf=nplanet*10+8
C     Store all the observations in the master file
      do 17 i=1,npt !central database of all data
        dtype(i)=0 !0 marks that we have photometric data
        aT(i)=time(i)
        aM(i)=flux(i)
        aE(i)=ferr(i)
        aIT(i)=exptime(i)
c        write(0,*) "exptime",exptime(i)
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

C     Get density constraints..
      nfrho=1
      if(iargc().ge.5)then
        call getarg(5,rhofile)
        if(rhofile.eq.'null') goto 25
        open(unit=nunit,file=rhofile,status='old',err=905)
        read(nunit,*) (rhoread(i),i=1,8)
        close(nunit)
        do 26 i=1,8
            rhoread(i)=rhoread(i)*1000.0d0
 26     continue
        rhoi=rhoread(1)
        rhoierr(1)=-rhoi
        rhoierr(2)=rhoread(8)
        rhoierr(3)=rhoread(6)
        rhoierr(4)=rhoread(4)
        rhoierr(5)=0.0d0
        rhoierr(6)=rhoread(3)
        rhoierr(7)=rhoread(5)
        rhoierr(8)=rhoread(7)
        rhoierr(9)=rhoierr(8)*10.0d0
        if(rhoread(2).gt.0.0d0) then 
            nfrho=0
            write(0,*) "density constraints are ON"
        endif
 25     continue
      endif
      if(nfrho.eq.1) write(0,*) "density constraints are OFF"

      do 27 i=1,nplanet
        if(iargc().ge.6+i)then
            call getarg(6+i,ttfile)

            if(ttfile.eq.'null')then
                ntt(i)=0
            else
                nunit=10
                open(unit=nunit,file=ttfile,status='old',err=906)
                call readttfile(nunit,nplanetmax,nmax,i,ntt,tobs,omc)
                close(nunit)
            endif
            
        else
            ntt(i)=0
        endif
c        write(0,*) "ntt",i,ntt(i)
 27   continue


cc     Initialize random seed
c      call itime(now)
c      seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
c      write(0,*) "Seed: ",seed
c      dumr=ran2(-seed)

      call transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npta,aT,aIT,
     .  ntt,tobs,omc,tmodel,dtype)

c     Compute chi^2
      chisq=0.0d0
      do 24 i=1,npta
        chisq=chisq+(tmodel(i)-aM(i))*(tmodel(i)-aM(i))/(aE(i)*aE(i))
 24   continue
      write(0,*) "chisq:",chisq
      chifac=dble(npta-1)/chisq

c     Iterate (change niter to input parameter)
      call getarg(4,cline)
      read(cline,*) niter
      
      dniter=dble(niter)
c     Update parameters
c      do 20 i=1,nf
      do 20 i=1,1
        do 23 j=1,nf
            sol2(j)=sol(j)
            if(serr(j,2)-serr(j,1).gt.0.0d0) then
                serr2(j,2)=-1.0d0
            else
                serr2(j,2)=0.0d0
            endif
c            write(0,*) serr(j,2)-serr(j,1),serr2(j,2)
 23     continue
cc        read(5,*) 
        serr2(i,2)=0.0d0
        if(serr(i,2)-serr(i,1).gt.0.0d0)then
c            niter1=int(dniter*(sol(i)-serr(i,1))/(serr(i,2)-serr(i,1)))
            niter1=int(dniter*(serr(i,2)-sol(i))/(serr(i,2)-serr(i,1)))
            dniter1=dble(niter1)
            
c            write(0,*) "hello",niter1
            do 19 j=1,niter1               
c                sol2(i)=(serr(i,2)-serr(i,1))*dble(j)/dniter+serr(i,1)
                sol2(i)=(serr(i,2)-sol(i))*dble(j)/dniter1+sol(i)

                call fittransitmodel2(npta,nfit,sol2,sol3,serr2,fvec,
     .            iwa,lwa,wa)
     
                call transitmodel(nfit,nplanet,nplanetmax,sol2,nmax,
     .           npta,aT,aIT,ntt,tobs,omc,tmodel,dtype)


c               Compute chi^2
                chisq=0.0d0
                do 21 k=1,npta
                    chisq=chisq+(tmodel(k)-aM(k))*(tmodel(k)-aM(k))/
     .                  (aE(k)*aE(k))
 21             continue
                drho=1.0d3*sol2(1)-rhoi
                call getrhosig(rhoierr,rhoin,9,drho,dsig)
                chisq=chisq*chifac+dsig*dsig

c               Output chi^2, sol2(i) 
                write(0,22) chisq,i,j,sol2(i)
                write(6,22) chisq,i,j,(sol2(k),k=1,nf)
 22             format(1PE17.10,2(1X,I3),1X,108(1PE17.10,1X))

 19         continue
        endif
 20   continue

     
      goto 999
 901  write(0,*) "Usage: transitscan5 <photfile> <rvfile> <fitpars> <nit
     .er> [rhofile] [kmag] [ttfiles]"
      goto 999
 902  write(0,*) "Error opening ",inputsol
      goto 999
 903  write(0,*) "Error opening ",obsfile
      goto 999
 904  write(0,*) "Error opening ",rvfile
      goto 999
 905  write(0,*) "Error opening ",rhofile
      goto 999
 906  write(0,*) "Error opening ",ttfile
      goto 999
 999  end
 
