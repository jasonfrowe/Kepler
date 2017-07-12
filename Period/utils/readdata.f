CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readdata(kID,nmax,npt,time,mag,merr,Keplertime,props)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer kID,npt,nunit,nbzi,nT,nstar,i,nmax,loop
      parameter(nT=476,nstar=52496)
      integer rkID(nstar),Teff(nstar),t(nstar),s(nstar),c(nstar),
     .  m(nstar),o(nstar),xcoo(nstar),ycoo(nstar)
      double precision time(nmax),mag(nmax),merr(nmax),Keplertime,
     .  props(5),mintime
      real rtime(nT),flux(nT,nstar),median(nstar),mad(nstar),
     .  madiff(nstar),Kmag(nstar),crowd(nstar),logg(nstar),rad(nstar)
      character*80 filename

      Keplertime=54900.0d0

      nunit=10
      filename="/home/rowe/Kepler/lcdata/nflux.dat"
      
C     Open the binary data file
      nbzi=nT*(nstar+1)
      open(unit=nunit,file=filename,status='old',form='unformatted',
     .  access='direct',recl=nbzi,err=901)
      read(nunit,rec=1) rtime,flux
      close(nunit)
      
      filename="/home/rowe/Kepler/lcdata/cdpp.idse"
      open(unit=nunit,file=filename,status='old',err=901)
      do 10 i=1,nstar
        read(nunit,*) rkID(i),t(i),Kmag(i),s(i),c(i),m(i),o(i),xcoo(i),
     .      ycoo(i),crowd(i),Teff(i),logg(i),rad(i)
 10   continue
      close(nunit)

      filename="/home/rowe/Kepler/lcdata/medMadMadDiffNflux.dat"      
C     Open the binary data file
      nbzi=nstar*3
      open(unit=nunit,file=filename,status='old',form='unformatted',
     .  access='direct',recl=nbzi,err=901)
      read(nunit,rec=1) median,mad,madiff
      close(nunit)    
      
      loop=0
      i=1
      do while(loop.eq.0)
        if(i.gt.nstar)loop=-1
        if(kID.eq.rkID(i))loop=i
        i=i+1
      enddo
      
      if(loop.eq.-1)then
        write(0,*) "KeplerID not found in database"
        stop
      endif
     
      mintime=99.9e30
      do 11 i=1,nT
        time(i)=real(rtime(i))
c  reset times form arbitrary start of 0.0 to HJD of field center
c  add 53.038152 to MJD center 1st exp
        time(i)=time(i)+53.038152
c  add HJD 5/2/09 adjustment for field center
        time(i)=time(i)+0.50020
c  add time dependent change for field center
        time(i)=time(i)+4.1d-5*(time(i)-53.0)
        mintime=min(mintime,time(i))
        mag(i)=real(flux(i,loop))!*median(loop)+median(loop)
        merr(i)=1.0
 11   continue
      npt=nT
      
      Keplertime=Keplertime+mintime  !correct time=0 time
      do 30 i=1,npt
         time(i)=time(i)-mintime
 30   continue
      
      props(1)=dble(Kmag(loop))
      props(2)=dble(Teff(loop))
      props(3)=dble(logg(loop))
      props(4)=dble(rad(loop))
      props(5)=dble(t(loop))
     
     
      goto 999
 901  write(0,*) "Cannot open ",filename
      stop
 999  return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readkeplc(nunit,nmax,npt,dtime,mag,merr,itime,
     .  Keplertime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nmax,npt,i
      double precision dtime(nmax),mag(nmax),merr(nmax),itime(nmax),
     .  Keplertime,sec2day,mintime
            
C     if time=0, then time=MOSTtime (that does not mean drink up)
      Keplertime=54900.0

C     number of seconds in a day
      sec2day=86400.0d0
      
      i=1
      
      mintime=99.9d30
  9   continue
 10   read(nunit,*,err=9,end=20) dtime(i),mag(i),merr(i)
cc  reset times form arbitrary start of 0.0 to HJD of field center
cc  add 53.038152 to MJD center 1st exp
c        dtime(i)=dtime(i)+53.038152
cc  add HJD 5/2/09 adjustment for field center
        dtime(i)=dtime(i)+0.5000-Keplertime
cc  add time dependent change for field center
c        dtime(i)=dtime(i)+4.1d-5*(dtime(i)-53.0)
        mintime=min(mintime,dtime(i))
c        mag(i)=-2.5*log10(mag(i)+1.0d0)
c        merr(i)=0.00005
        itime(i)=1800.0/sec2day
c        merr(i)=0.0002
        i=i+1
      goto 10
 20   continue
        
      npt=i-1
      write(0,*)   "-------------------------"
      write(0,500) "Observations read: ",npt
      write(0,*) "Mintime: ",mintime
      write(0,*)   "-------------------------"
 500  format(1X,A19,I6)
 
      Keplertime=Keplertime+mintime  !correct time=0 time
      do 30 i=1,npt
         dtime(i)=dtime(i)-mintime
 30   continue
 
      return
      end
  