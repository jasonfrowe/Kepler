CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program transittiming
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Fits for individual transit times
      implicit none
      integer iargc,nfit,npt,nmax,i,nunit,nptv,npta,j,k,flag
      parameter(nfit=17,nmax=2000000)
      integer ia(nfit),dtype(nmax)
      double precision sol(nfit),time(nmax),dt,tmodel(nmax),
     .  flux(nmax),ferr(nmax),itime(nmax),Keplertime,serr(nfit,2),
     .  err(nfit),covar(nfit,nfit),alpha(nfit,nfit),dchistop,
     .  vtime(nmax),vel(nmax),verr(nmax),vetime(nmax),aT(nmax),aM(nmax),
     .  aE(nmax),aIT(nmax),kmag,kerr,sol2(nfit),tdur,T0,Tmin,Tmax,Ts,Te,
     .  Ts2,Te2,ttold,ttold2,Pi,tPi,eccn,b,Psec,adrs,cincl,temp(4)
      character*3 titles(nfit)
      character*80 inputsol,obsfile,rvfile,cline
      common /Fitting/ npta,aIT
      
      flag=1 !if flag=0 do not predict next time, flag=1 predict next
      
      Pi=acos(-1.d0)   !Pi
      tPi=2.0d0*Pi   !Pi
      
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
      
      eccn=sqrt(sol(7)*sol(7)+sol(8)*sol(8)) !eccentricity
      if(eccn.ge.1.0) eccn=0.99
      b=sqrt(sol(3))
      Psec=sol(2)*8.64d4 !sec ; period of planet
      adrs=sol(5)*sol(2)/tpi*sqrt((1.0d0+sol(4))**2.0d0-sol(3))*
     .  (1+sol(8))/sqrt(1-eccn*eccn)
      cincl=b/adrs !cos(i)
        
      temp(1)=Psec/Pi
      temp(2)=1.0d0/adrs
      temp(3)=(1+sol(4))**2.0-sol(3)
      temp(4)=1-cincl*cincl
C     Transit duration in hours
      tdur=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))/8.64d4
c      tdur=2.0d0/sol(5) !transit duration
c      write(0,*) "tdur: ",tdur
c      read(5,*)

      tdur=max(tdur,2.0d0/24.0) !make sure duration is at least 2.0hours
      write(0,*) "Tdur: ",tdur*24.0

C     I think there is a bug with this line.. changed 0.5 to 0.0      
      T0=sol(1)+int((Tmin-sol(1))/sol(2)+0.0d0)*sol(2)
      write(0,*) "T0, Per:",T0,sol(2)
      write(0,*) "Tmax: ",Tmax
c      read(5,*)
      if(sol(2).le.0.0)then
        write(0,*) "Period is non-sense.. exiting.."
        goto 999
      endif
      
      ttold=0.0 !use timing offset from previous measurement
c       ttold=-0.1
c      ttold=1.5

      do while(T0.lt.Tmax)
c         if(T0.gt.1340.0) ttold=-0.01*sol(2)
c         if(T0.gt.1470.0) ttold=-0.02*sol(2)
c         write(6,*) "ttold ",ttold
c        write(0,*) T0
C     Section of data we want to fit:      
        Ts=T0-2.0*tdur+ttold
        Te=T0+2.0*tdur+ttold
        Ts2=T0-0.5*tdur+ttold
        Te2=T0+0.5*tdur+ttold
c       write(0,*) "ttold:",ttold
c       read(5,*)
      
C       We only fit the epoch
        do 12 i=2,nfit
            serr(i,2)=0.0d0 !disable all parameter fits, except epoch
 12     continue
        serr(1,2)= -1.0d0
        serr(6,2)= -1.0d0
C       Copy all-fit solution
        do 13 i=1,nfit
            sol2(i)=sol(i)
 13     continue
        sol2(1)=sol(1)+ttold
 
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
c        write(0,*) "Number of points ",npta
        do 18 i=1,nptv
            dtype(npt+i)=1 !mark data as RV
            aT(npt+i)=vtime(i)
            aM(npt+i)=vel(i)
            aE(npt+i)=verr(i)
            aIT(npt+i)=vetime(i)
 18     continue
        npta=npta+nptv

        if(k.ge.4)then
            write(0,*) "k:",k
            dchistop=1.0d-10
c            write(0,*) (serr(i,2),i=1,nfit)
c            read(5,*)
            call fittransitmodel(npta,aT,aM,aE,dtype,nfit,sol2,serr,ia,
     .          covar,alpha,dchistop,err)
            ttold2=sol2(1)-sol(1)-ttold
            ttold=sol2(1)-sol(1)            
            write(6,*) T0,sol2(1)-sol(1),err(1)!,sol(4),err(4)
 500  format(19(1X,1PE17.10)) 
        else
            ttold=ttold+ttold2
        endif
        
        if(flag.eq.0) ttold=0.0d0  !check if we are using predictive
        T0=T0+sol(2)
      enddo
     
     
c      call exportfit(nfit,sol,serr,err,titles)
c
c      call transitmodel(nfit,sol,npta,aT,aIT,tmodel,dtype)
c      do 10 i=1,npta
c        write(6,*) aT(i),aM(i),tmodel(i)
c 10   continue
      
      goto 999
 901  write(0,*) "Usage: transitfit4 <photfile> <rvfile> <fitpars>"
      goto 999
 902  write(0,*) "Error opening ",inputsol
      goto 999
 903  write(0,*) "Error opening ",obsfile
      goto 999
 904  write(0,*) "Error opening ",rvfile
      goto 999
 999  end
 
