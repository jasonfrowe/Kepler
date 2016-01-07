      program transittiming5
      implicit none
      integer nfit,nmax,iargc,nunit,npt,i,nptv,nplanet,npta,np,
     .   nplanetmax,j,k,flag
      parameter(nfit=108,nmax=2000000,nplanetmax=10)
      integer dtype(nmax),ntt(nplanetmax),ia(nfit)
      double precision time(nmax),flux(nmax),ferr(nmax),itime(nmax),
     .  Keplertime,kmag,kerr,vtime(nmax),vel(nmax),verr(nmax),
     .  vetime(nmax),sol(nfit),serr(nfit,2),err(nfit),aT(nmax),
     .  aM(nmax),aE(nmax),aIT(nmax),tdur,transitdur,Tmin,Tmax,T0,
     .  tobs(nplanetmax,nmax),omc(nplanetmax,nmax),ttcor,ttold,
     .  tcor(nmax),Ts,Te,Ts2,Te2,sol2(nfit),dchistop,ttold2,
     .  covar(nfit,nfit),alpha(nfit,nfit)
      character*80 obsfile,cline,rvfile,inputsol
      common /Fitting/ npta,nplanet,aIT,ntt,tobs,omc

      flag=0 !if flag=0 do not predict next time, flag=1 predict next
      
      if(iargc().lt.3) goto 901
      
C     Parse the name of the observations data file from the commandline
      call getarg(1,obsfile)
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,flux,ferr,itime,Keplertime)
      call readkeplc(nunit,nmax,npt,time,flux,ferr,itime,Keplertime)
      close(nunit)!release unit number as we are done with the file.

C     get range of photometric observations
      Tmin=time(1)
      Tmax=time(1)
      do 14 i=2,npt
        Tmin=min(time(i),Tmin)
        Tmax=max(time(i),Tmax)
 14   continue

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
      if(nplanet.gt.1) goto 905 !only want 1 planet at this point
      
C     Store all the observations in the master file
c      do 17 i=1,npt !central database of all data
c        dtype(i)=0 !0 marks that we have photometric data
c        aT(i)=time(i)
c        aM(i)=flux(i)
c        aE(i)=ferr(i)
c        aIT(i)=itime(i)
cc        write(0,*) "itime",itime(i)
c 17   continue
c      npta=npt
c      do 18 i=1,nptv
c        dtype(npt+i)=1 !mark data as RV
c        aT(npt+i)=vtime(i)
c        aM(npt+i)=vel(i)
c        aE(npt+i)=verr(i)
c        aIT(npt+i)=vetime(i)
c 18   continue
c      npta=npta+nptv
      
      do 10 i=1,nplanet
         ntt(i)=0
 10   continue
      
      np=1
      tdur=transitdur(nfit,sol,np)/8.64d4
c      write(0,*) "Tdur: ",tdur
      tdur=max(tdur,2.55d0/24.0) !make sure duration is at least 2.5 hours
      write(0,*) "Tdur: ",tdur*24.0
c      read(5,*)

C     I think there is a bug with this line.. changed 0.5 to 0.0
      T0=sol(9)+int((Tmin-sol(9))/sol(10)+0.0d0)*sol(10)
      write(0,*) "T0, Per:",T0,sol(10)
      write(0,*) "Tmax: ",Tmax
c      read(5,*)
      if(sol(10).le.0.0)then
        write(0,*) "Period is non-sense.. exiting.."
        goto 999
      endif

      ttold=0.15 !use timing offset from previous measurement
c       ttold=-0.5
c      ttold=1.5

      do while(T0.lt.Tmax)
c         if(T0.gt.1340.0) ttold=-0.01*sol(2)
c         if(T0.gt.1470.0) ttold=-0.02*sol(2)
c        if((T0.gt.300.0).and.(T0.lt.400.0)) ttold=0.1
c         write(6,*) "ttold ",ttold
c        write(0,*) T0
C     Section of data we want to fit:
        Ts=T0-2.0*tdur+ttold
        Te=T0+2.0*tdur+ttold
        Ts2=T0-0.5*tdur+ttold
        Te2=T0+0.5*tdur+ttold

C       We only fit the epoch
        do 12 i=1,nfit
            serr(i,2)=0.0d0 !disable all parameter fits, except epoch
 12     continue
        serr(8,2)= -1.0d0
        serr(9,2)= -1.0d0
C       Copy all-fit solution
        do 13 i=1,nfit
            sol2(i)=sol(i)
 13     continue
        sol2(9)=sol(9)+ttold

C     Store all the observations in the master file
        j=0
        k=0 !count points inside expected transit.
        do 17 i=1,npt !central database of all data
            if((time(i).ge.Ts).and.(time(i).le.Te))then
                if((time(i).ge.Ts2).and.(time(i).le.Te2)) k=k+1
                j=j+1
                dtype(j)=0 !0 marks that we have photometric data
                aT(j)=time(i)
                aM(j)=flux(i)
                aE(j)=ferr(i)
                aIT(j)=itime(i)
            endif
 17     continue
        npta=j

        if(k.ge.4)then
            write(0,*) "k:",k
            dchistop=1.0d-10
c            call fittransitmodel(npta,aT,aM,aE,dtype,nfit,sol2,serr,ia,
c     .          covar,alpha,dchistop,err)
            call fittransitmodel(npta,aT,aM,aE,dtype,nfit,nplanet,sol2,
     .         serr,ia,covar,alpha,dchistop,err)
            ttold2=sol2(9)-sol(9)-ttold
            ttold=sol2(9)-sol(9)
            write(6,*) T0,sol2(9)-sol(9),err(9)!,sol(4),err(4)
 500  format(19(1X,1PE17.10))
        else
            ttold=ttold+ttold2
        endif

        if(flag.eq.0) ttold=0.0d0  !check if we are using predictive
        T0=T0+sol(10)


      enddo
      
      goto 999
 901  write(0,*) "Usage: transittiming5 photometry RVs n1.dat <Kmag>"
      write(0,*) ""
      write(0,*) "photometry: Kepler photometry"
      write(0,*) "RVs : radial velocity measurements"
      write(0,*) "n1.dat : transit model solution"
C      write(0,*) "data.tt : transit timing measurements"
      write(0,*) "Kmag : Kepler magnitude (optional)"
      goto 999
 902  write(0,*) "Cannot open ",inputsol
      goto 999
 903  write(0,*) "Cannot open ",obsfile
      goto 999
 904  write(0,*) "Cannot open ",rvfile
      goto 999
 905  write(0,*) "This code can only handle one planet at a time. Sorry"
      write(0,*) " "
      goto 901
 999  end
