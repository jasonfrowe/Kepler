      program detrend
C     (c) Jason Rowe (2010) jasonfrowe@gmail.com
      implicit none
      integer nmax,iargc,nunit,npt,i,nfit
      parameter(nmax=650000)
      double precision time(nmax),flux(nmax),ferr(nmax),boxbin,sig,
     .  rmavg
      character*80 filename,cline
      
      if(iargc().lt.2) goto 901
      call getarg(1,filename)
      call getarg(2,cline)
      read(cline,*) boxbin
      if(boxbin.lt.0) goto 901
      
      nunit=10
      open(unit=nunit,file=filename,status='old',err=902)
      call readkeplc(nunit,nmax,npt,time,flux,ferr)
      close(nunit)
      if(npt.eq.0) goto 999
      
      nfit=3
      call polydetrend(npt,time,flux,ferr,nfit)   
c      goto 9   
      
      if(boxbin.gt.0) call boxcar(npt,time,flux,ferr,boxbin)

      sig=3.0
c      call sigup(npt,time,flux,ferr,sig,rmavg)
      
 9    do 10 i=1,npt
        write(6,*) time(i),flux(i),
     .      ferr(i)
 10   continue
      
      goto 999
 901  write(0,*) "Usage: detrend <filein> <boxbin>"
      write(0,*) "<boxbin> width of boxcar in days"
      goto 999
 902  write(0,*) "Cannot open ",filename
 999  end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine sigup(npt,time,mag,merr,sig,rmavg)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nmax,ntmp,i,nclip,niter,nitmax
      parameter (nmax=650000,nitmax=50)
      double precision time(npt),mag(npt),merr(npt),sig,std,stdev,
     .  tmp1(nmax),tmp2(nmax),tmp3(nmax),mean,dum,tmp4(nmax),omean,rmavg


C     watch out for infinite loops
      niter=0

C     find mean of data set
      ntmp=0
      mean=0.
      do 4 i=1,npt
c         if(mag(i).lt.90.0) then
            ntmp=ntmp+1
            mean=mean+mag(i)
c         endif
 4    continue
c      write(0,*) mean,ntmp
      mean=mean/dble(ntmp)
C     find standard dev. of data set.
      std=stdev(npt,mag,mean)

      write(0,*) "sigclip:",npt,mean,std

C     count number of clipped points
 6    nclip=0
C     count number of new points
      ntmp=0
      do 10 i=1,npt
c         dum=abs(mag(i)-mean)
         dum=mag(i)-mean !only clip in + direction
         if(dum.lt.sig*std) then
            ntmp=ntmp+1
            tmp1(ntmp)=time(i)
            tmp2(ntmp)=mag(i)
            tmp3(ntmp)=merr(i)
         else
            nclip=nclip+1
         endif
 10   continue
C     if nothing is clipped, no point in continuing
      if(nclip.eq.0) goto 15
      
C     save copy of old mean
      omean=mean

C     find mean of new data set
      mean=0.
      do 5 i=1,ntmp
         mean=mean+tmp2(i)
 5    continue
      mean=mean/dble(ntmp)
C     if mean doesn't change, we're done.
      if(mean.eq.omean) goto 15

C     if mean doesn't change, we're done
      if(nclip.gt.0) goto 15 
C     find st. dev. of new data set
      std=stdev(ntmp,tmp2,mean)
C     now restart clipping
      niter=niter+1
C     if we loop too much.. get out. (probably an infinite loop)
      if(niter.gt.nitmax) goto 15
      goto 6

 15   do 20 i=1,ntmp
         time(i)=tmp1(i)
         mag(i)=tmp2(i)-mean
         merr(i)=tmp3(i)
 20   continue
      npt=ntmp
C     update zero point
      rmavg=rmavg+mean

      return
      end