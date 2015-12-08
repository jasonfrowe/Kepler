      program frequencyremover
      implicit none
      integer addsub,nunit,npt,nmax,nfreq,nunit2,iargc,ma,nstar,i
      parameter(nmax=650000)
      double precision time(nmax),mag(nmax),merr(nmax),
     .   freq(nmax),amp(nmax),ph(nmax),ztime,a(nmax),sns(nmax),
     .  mintime
      character*80 filename,cline,outname,freqname
      
      if(iargc().lt.4)goto 900
      call getarg(1,filename)
      call getarg(2,cline)
      read(cline,*) addsub
      if((addsub.lt.0).or.(addsub.gt.1))goto 900
      call getarg(3,freqname)
      call getarg(4,outname)
      
      
      nunit=10 !set unit number for file read
      
      open(unit=nunit,file=filename,status='old',err=901)
      
C     open file with frequencies, amplitudes and phases
C     from a cosine series
      nunit2=12
      open(unit=nunit2,file=freqname,status='old',err=901)
    
C     read in photometry data  
      call readdata(nunit,nmax,npt,time,mag,merr)
C     we have photometry, can close the file now.
      close(nunit)
      
      mintime=time(1)
      do 11 i=1,npt !convert to ppt (~mmag)
        mag(i)=mag(i)*1000.0
        merr(i)=merr(i)*1000.0
        mintime=min(time(i),mintime)
 11   continue
 
      do 13 i=1,npt
        time(i)=time(i)-mintime
 13   continue
      
C     read in frequency information.
      nfreq=nmax
      call readfreqs(nunit2,ma,nstar,a,sns,ztime)
C     have frequencies, can close file now.
      close(nunit2)
      
      call freqaddsub(addsub,npt,time,mag,ma,nstar,a,ztime)
      
      do 12 i=1,npt
        mag(i)=mag(i)/1000.0
        merr(i)=merr(i)/1000.0
        time(i)=time(i)+mintime
 12   continue
      call exportdata(npt,time,mag,merr,outname)
      
      goto 999
 900  write(6,*) "Usage: freqremove <fname> <addsub> <freqdat> <outdat>"
        write(6,*) "<fname>: Kepler Photometry datafile"
        write(6,*) "<addsub>: 0 - add frequencies"
        write(6,*) "          1 - remove frequencies"
        write(6,*) "<freqdat>: Frequency output from MOSTPer"
        write(6,*) "<outdat>: Output name for photometry"
        goto 999
 901  write(6,*) "Cannot open ",filename
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine exportdata(npt,time,mag,merr,outname)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i
      double precision time(npt),mag(npt),merr(npt)
      character*80 outname
      
      open(unit=11,file=outname)
      
      do 10 i=1,npt
c        mag(i)=10.0**(mag(i)/-2.5d0)
        write(11,*) time(i),mag(i),merr(i)
 10   continue
 
      close(11)
 
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readdata(nunit,nmax,npt,time,mag,merr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,npt,nmax,i
      double precision time(nmax),mag(nmax),merr(nmax)
      
      i=1
 10   read(nunit,*,end=11) time(i),mag(i),merr(i)
c        mag(i)=-2.5*log10(mag(i)+1.0d0)
c        merr(i)=0.0001
        i=i+1
        goto 10
 11   continue
      npt=i-1
      
      return
      end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine freqaddsub(addsub,npt,time,mag,ma,nstar,a,ztime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer addsub,npt,i,j,ma,nstar,i2,k,cc
      double precision time(npt),mag(npt),Pi,tPi,ztime,ttime,temp,
     .  a(npt),arg
      Pi=3.141592654
      tPi=2.0*Pi
      

      j=(ma-nstar)/nstar
      do 10 i2=1,npt
         temp=0.0
         ttime=time(i2)-ztime
         do 20 k=1,nstar
            cc=(j+1)*(k-1)+2
            do 21 i=cc,cc+j-2,2
               arg=dble((i-cc+2)/2)*a(cc-1)*ttime+a(i+1)
               temp=temp+a(i)*cos(arg)
 21         continue
 20      continue
         if(addsub.eq.0) then
            mag(i2)=mag(i2)+temp
         elseif(addsub.eq.1) then
            mag(i2)=mag(i2)-temp
         endif
 10   continue
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readfreqs(nunit,ma,nstar,a,sns,ztime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,i,ma,nstar,k,j,cc,nmax
      parameter(nmax=650000)
      double precision ztime,a(nmax),twopi,pi,sns(nmax)
      
      pi=3.141592654
      twopi=2.0*pi
      
      read(nunit,*) ztime
      read(nunit,*) ma,nstar
      j=(ma-nstar)/nstar
      do 10 k=1,nstar
         cc=(j+1)*(k-1)+2
         read(nunit,*) a(cc-1),(a(i),i=cc,cc+j-2,2),
     .        (a(i+1),i=cc,cc+j-2,2),sns(k)
c         write(6,500) a(cc-1),(a(i),i=cc,cc+j-2,2),
c     .        (a(i+1),i=cc,cc+j-2,2),sns(k) 
         a(cc-1)=a(cc-1)*twopi
 500     format(20(E14.7,1X))
 10   continue
      
      return
      end
      
      