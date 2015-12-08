      program detrend
C     (c) Jason Rowe (2010) jasonfrowe@gmail.com
      implicit none
      integer nmax,iargc,nunit,npt,i,nfit,nfitp,ngap,ftype
      parameter(nmax=650000,nfit=18)
      integer ts(nmax),tflag(nmax),nx(nmax)
      double precision time(nmax),flux(nmax),ferr(nmax),boxbin,sig,
     .  rmavg,t(nmax),f(nmax),fe(nmax),sol(nfit),serr(nfit,2),
     .  Dpvary(nfit),err(nfit),doe,toff,eoff,phase(nmax),x(nmax),
     .  y(nmax),z(nmax),gaps(nmax),offset(nmax),work(nmax)
      character*80 filename,cline,fitparsfile
      
      if(iargc().lt.3) goto 901
      call getarg(1,cline)
      read(cline,*,err=901) ftype
      call getarg(2,filename)
      call getarg(3,cline)
      read(cline,*) boxbin
      if(boxbin.lt.0) goto 901
            
      nunit=10
      open(unit=nunit,file=filename,status='old',err=902)
      call readkeplc(nunit,nmax,npt,time,flux,ferr)
      close(nunit)
      if(npt.eq.0) goto 999
      
      do 11 i=1,npt
        tflag(i)=0
 11   continue
      
      do 12 i=4,14 !upto 10 planets
        if(iargc().ge.i) then
            write(0,*) "Planet ",i-3
            call getarg(i,fitparsfile)
            nunit=10
            open(unit=nunit,file=fitparsfile,status='old',err=904)
           call getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff,eoff)
            close(nunit)
            call marktransit(npt,phase,time,flux,tflag,nfit,sol)
        endif
 12   continue
      
      nfitp=4
      
      if(ftype.eq.0)then
        call polyfilter(npt,time,flux,ferr,ts,tflag,boxbin,x,y,z,ngap,
     .      gaps,offset,nfitp,work)
 9      do 10 i=1,npt
c            if(tflag(i).eq.0) write(6,*) time(i)-0.5d0+54900.0d0,
c     .          flux(i),ferr(i)
            write(6,*) time(i)-0.5d0+54900.0d0,
     .          flux(i),ferr(i)
 10     continue
      else
        call boxfilter(npt,time,flux,ferr,ts,ngap,gaps,tflag,x,y,z,
     .      boxbin,nx)
      endif
      

      
      goto 999
 901  write(0,*) "Usage: detrend <ftype> <filein> <boxbin> [f1.dat]"
      write(0,*) "<ftype> 0=polyfilter 1=boxcar"
      write(0,*) "<filein> Kepler photometry"
      write(0,*) "<boxbin> width of filter in days"
      write(0,*) "[f1.dat] transit solution (optional)"
      goto 999
 902  write(0,*) "Cannot open ",filename
      goto 999
 904  write(6,*) "Cannot open ",fitparsfile
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine boxfilter(npt,time,flux,ferr,ts,ngap,gaps,tflag,x,y,z,
     .  boxbin,nx)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,ts(npt),tflag(npt),ngap,i,nc,nx(npt),j,k,nfitp
      double precision time(npt),flux(npt),ferr(npt),gaps(npt),x(npt),
     .  y(npt),z(npt),cadence,boxbin,t1,t2
      
      nc=15.0d0 !number of cadences that need to be missing to mark a gap
      cadence=0.02*nc
      
C     sort the data by time
      call rqsort(npt,time,ts)
      
      ngap=0
      do 12 i=1,npt-1
        if(time(ts(i+1))-time(ts(i)).gt.cadence)then
            ngap=ngap+1
            gaps(ngap)=(time(ts(i+1))+time(ts(i+1)))/2.0d0
c            write(0,*) ngap,gaps(ngap)
        endif
 12   continue
 
      do 13 i=1,ngap+1
        if(i.eq.1)then
            t1=time(ts(1))
            if(ngap.eq.0)then
                t2=time(ts(npt))+0.001
            else
                t2=gaps(1)
            endif
        elseif((i.eq.ngap+1).and.(ngap.gt.0))then
            t1=gaps(ngap)
            t2=time(ts(npt))+0.001
        else
            t1=gaps(i-1)
            t2=gaps(i)
        endif    
        j=0
        do 14 k=1,npt
            if((time(k).ge.t1).and.(time(k).lt.t2))then
                j=j+1
                x(j)=time(k)
                y(j)=flux(k)
                z(j)=ferr(k)
                nx(j)=tflag(k)
            endif
 14     continue
c        write(0,*) "j",j,t1,t2
        nfitp=3
        if(j.gt.nfitp+1) call polydetrend2(j,x,y,z,nfitp)
        call boxcar(j,x,y,z,boxbin,tflag,nx)
        do 10 k=1,j
            write(6,*) x(k)-0.5d0+54900.0d0,
     .          y(k),z(k)
 10     continue
 13   continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine polydetrend2(npt,time,mag,merr,nfit)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nfit
      integer i,ma,j,npc
      parameter(npc=10)
      integer ia(npc)
      real*8 time(npt),mag(npt),merr(npt),ans(npc),covar(npc,npc),
     .     chisq,T,meanT

      if(nfit.gt.npc) pause "increase npc in detrend"
      do 36 i=1,nfit
         ia(i)=1
 36   continue
 
      meanT=0.0d0
      do 11 i=1,npt
        meanT=meanT+time(i)
 11   continue
      meanT=meanT/dble(npt)

      do 12 i=1,npt
        time(i)=time(i)-meanT
 12   continue

      ma=nfit
      call lfit(time,mag,merr,npt,ans,ia,ma,covar,npc,chisq)
      write(0,*) (ans(i),i=1,nfit)

      do 10 i=1,npt
         T=0.
         DO 33 J = 1, NFIT
            T = ANS(J)*((time(i))**REAL(J-1)) +  T
 33      CONTINUE
         mag(i)=mag(i)-T
 10   continue

      do 13 i=1,npt
        time(i)=time(i)+meanT
 13   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine phasept(npt,time,phase,period,toff)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     npt - number of data points
C     time - times of observations
C     mag - the observations
C     phase - returned phase of observation
C     period - fixed period for data
      implicit none
      integer npt
      double precision time(npt),phase(npt),period,toff

      integer i
      double precision temp

      do 10 i=1,npt
         temp=time(i)
C        Get the phase
         phase(i)=temp/period-int(temp/period)
C        apply optional phase offset to make plot pretty
         phase(i)=phase(i)+toff
C        make sure phase is between 0 and 1
         if(phase(i).lt.0.0) phase(i)=phase(i)+1.0
         if(phase(i).gt.1.0) phase(i)=phase(i)-1.0
 10   continue
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readkeplc(nunit,nmax,npt,dtime,flux,ferr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nmax,npt,i
      double precision dtime(nmax),flux(nmax),ferr(nmax),Keplertime

      Keplertime=54900.0d0

      i=1
      
  9   continue
 10   read(nunit,*,err=9,end=20) dtime(i),flux(i),ferr(i)
         dtime(i)=dtime(i)+0.5d0-Keplertime
c        mag(i)=-2.5*log10(mag(i)+1.0d0)
c        ferr(i)=0.00005
        i=i+1
      goto 10
 20   continue
        
      npt=i-1
      write(0,*)   "-------------------------"
      write(0,500) "Observations read: ",npt
      write(0,*)   "-------------------------"
 500  format(1X,A19,I6)
 
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,nunit,i
      double precision sol(nfit),serr(nfit,2),Dpvary(nfit),toff,
     .  err(nfit),doe
      character*160 command
      
c      write(0,*) "----------------------------"
c      write(0,*) "Reading in starting solution"
c      write(0,*) "----------------------------"
      
      i=0 !initialize line counter
      toff=0.0 !initialize toff to zero
C     Start of loop to read in file
  10  read(nunit,500,end=11) command
 500  format(A160) !command much be contained within first 80 characters
        i=i+1 !increase line counter
        if(command(1:1).eq."#") then
            write(0,*) "Comment line on line ",i
        elseif(command(1:3).eq."SMA") then
            read(command(5:160),*) sol(1),serr(1,1),serr(1,2),Dpvary(1),
     .          err(1)
c           write(0,501) "SMA: ",sol(1),serr(1,1),serr(1,2),Dpvary(1),
c     .          err(1)
        elseif(command(1:3).eq."PMA") then
            read(command(5:160),*) sol(2),serr(2,1),serr(2,2),Dpvary(2),
     .          err(2)
c           write(0,501) "PMA: ",sol(2),serr(2,1),serr(2,2),Dpvary(2),
c     .          err(2)
        elseif(command(1:3).eq."SRA") then
            read(command(5:160),*) sol(3),serr(3,1),serr(3,2),Dpvary(3),
     .          err(3)
c           write(0,501) "SRA: ",sol(3),serr(3,1),serr(3,2),Dpvary(3),
c     .          err(3)
        elseif(command(1:3).eq."PRA") then
            read(command(5:160),*) sol(4),serr(4,1),serr(4,2),Dpvary(4),
     .          err(4)
c           write(0,501) "PRA: ",sol(4),serr(4,1),serr(4,2),Dpvary(4),
c     .          err(4)
        elseif(command(1:3).eq."PER") then
            read(command(5:160),*) sol(5),serr(5,1),serr(5,2),Dpvary(5),
     .          err(5)
c           write(0,501) "PER: ",sol(5),serr(5,1),serr(5,2),Dpvary(5),
c     .          err(5)
        elseif(command(1:3).eq."INC") then
            read(command(5:160),*) sol(6),serr(6,1),serr(6,2),Dpvary(6),
     .          err(6)
c           write(0,501) "INC: ",sol(6),serr(6,1),serr(6,2),Dpvary(6),
c     .          err(6)
        elseif(command(1:3).eq."EPO") then
            read(command(5:160),*) sol(7),serr(7,1),serr(7,2),Dpvary(7),
     .          err(7)
c           write(0,501) "EPO: ",sol(7),serr(7,1),serr(7,2),Dpvary(7),
c     .          err(7)
        elseif(command(1:3).eq."ZPT") then
            read(command(5:160),*) sol(8),serr(8,1),serr(8,2),Dpvary(8),
     .          err(8)
c           write(0,501) "ZPT: ",sol(8),serr(8,1),serr(8,2),Dpvary(8),
c     .          err(8)
        elseif(command(1:3).eq."ALB") then
            read(command(5:160),*) sol(9),serr(9,1),serr(9,2),Dpvary(9),
     .          err(9)
c           write(0,501) "ALB: ",sol(9),serr(9,1),serr(9,2),Dpvary(9),
c     .          err(9)
        elseif(command(1:3).eq."NL1") then
            read(command(5:160),*) sol(10),serr(10,1),serr(10,2),
     .          Dpvary(10),err(10)
c         write(0,501) "NL1: ",sol(10),serr(10,1),serr(10,2),Dpvary(10),
c     .          err(10)
        elseif(command(1:3).eq."NL2") then
            read(command(5:160),*) sol(11),serr(11,1),serr(11,2),
     .          Dpvary(11),err(11)
c         write(0,501) "NL2: ",sol(11),serr(11,1),serr(11,2),Dpvary(11),
c     .          err(11)
        elseif(command(1:3).eq."NL3") then
            read(command(5:160),*) sol(12),serr(12,1),serr(12,2),
     .          Dpvary(12),err(12)
c         write(0,501) "NL3: ",sol(12),serr(12,1),serr(12,2),Dpvary(12),
c     .          err(12)
        elseif(command(1:3).eq."NL4") then
            read(command(5:160),*) sol(13),serr(13,1),serr(13,2),
     .          Dpvary(13),err(13)
c         write(0,501) "NL4: ",sol(13),serr(13,1),serr(13,2),Dpvary(13),
c     .          err(13)
        elseif(command(1:3).eq."ECW") then
            read(command(5:160),*) sol(14),serr(14,1),serr(14,2),
     .          Dpvary(14),err(14)
c         write(0,501) "ECN: ",sol(14),serr(14,1),serr(14,2),Dpvary(14),
c     .          err(14)
        elseif(command(1:3).eq."ESW") then
            read(command(5:160),*) sol(15),serr(15,1),serr(15,2),
     .          Dpvary(15),err(15)
c         write(0,501) "WWW: ",sol(15),serr(15,1),serr(15,2),Dpvary(15),
c     .          err(15)
        elseif(command(1:3).eq."TED") then
            read(command(5:160),*) sol(16),serr(16,1),serr(16,2),
     .          Dpvary(16),err(16)
c         write(0,501) "TED: ",sol(16),serr(16,1),serr(16,2),Dpvary(16),
c     .          err(16)
        elseif(command(1:3).eq."VOF") then
            read(command(5:160),*) sol(17),serr(17,1),serr(17,2),
     .          Dpvary(17),err(17)
        elseif(command(1:3).eq."DIL") then
            read(command(5:160),*) sol(18),serr(18,1),serr(18,2),
     .          Dpvary(18),err(18)
        elseif(command(1:3).eq."OFF") then
            read(command(5:160),*,err=12,end=12) toff,doe
            goto 13
 12         doe=0.0
            read(command(5:160),*) toff
 13     continue
c            write(0,501) "OFF: ",toff,doe
        endif
 501    format(A5,5(1X,1PE17.10))
      goto 10 !loop back to read statement
 11   continue !end loop when EOF is reached
c      write(0,*) "----------------------------"
             
      return
      end 

 
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
