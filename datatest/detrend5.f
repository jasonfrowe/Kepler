      program detrend5
C     (c) Jason Rowe (2010) jasonfrowe@gmail.com
      implicit none
      integer nmax,iargc,nunit,npt,i,nfit,nfitp,ngap,ftype,nplanetmax,
     .  nplanet
      parameter(nmax=2000000,nfit=108,nplanetmax=10)
      integer ts(nmax),tflag(nmax),nx(nmax),ntt(nplanetmax)
      double precision time(nmax),flux(nmax),ferr(nmax),boxbin,sig,
     .  rmavg,t(nmax),f(nmax),fe(nmax),sol(nfit),serr(nfit,2),
     .  Dpvary(nfit),err(nfit),doe,toff,eoff,phase(nmax),x(nmax),
     .  y(nmax),z(nmax),gaps(nmax),offset(nmax),work(nmax),
     .  tobs(nplanetmax,nmax),omc(nplanetmax,nmax)
      character*80 filename,cline,inputsol,ttfile
      
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

      nplanet=0
      if(iargc().ge.4) then
        call getarg(4,inputsol) !get filename for input solution
        nunit=10 !unit number used for file input
        open(unit=nunit,file=inputsol,status='old',err=904)
C       We start by reading in solution from input file
c        write(0,*) "reading in input solution"
        call getfitpars(nunit,nfit,nplanet,sol,serr,err)
c        write(0,*) "done reading input solution"
        close(nunit) !release unit number as we are done with file
      endif

      do 21 i=1,nplanet
        write(0,*) "i: ",i
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
      
      do 22 i=1,nplanet
        call marktransit(i,npt,phase,time,tflag,nfit,sol,ntt,tobs,omc)  
 22   continue

      
  
c      do 12 i=4,14 !upto 10 planets
c        if(iargc().ge.i) then
c            write(0,*) "Planet ",i-3
c            call getarg(i,fitparsfile)
c            nunit=10
c            open(unit=nunit,file=fitparsfile,status='old',err=904)
c           call getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff,eoff)
c            close(nunit)
c            call marktransit(npt,phase,time,flux,tflag,nfit,sol)
c        endif
c 12   continue
      
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
 904  write(6,*) "Cannot open ",inputsol
      goto 999
 905  write(0,*) "Error opening ",ttfile
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
      
      nc=5.0d0 !number of cadences that need to be missing to mark a gap
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
 500  format(1X,A19,I7)
 
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getfitpars(nunit,nfit,nplanet,sol,serr,err)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,nunit,i,nplanet,np,j
      double precision sol(nfit),serr(nfit,2),err(nfit)
      character*160 command
      
c      write(0,*) "----------------------------"
c      write(0,*) "Reading in starting solution"
c      write(0,*) "----------------------------"
      
      nplanet=0 !initialize number of planets
      
      i=0 !initialize line counter
C     Start of loop to read in file
  10  read(nunit,500,end=11,err=901) command
 500  format(A160) !command much be contained within first 80 characters
        i=i+1 !increase line counter
        if(command(1:1).eq."#") then
c            write(0,*) "Comment line on line ",i
        elseif(command(1:3).eq."RHO") then
            read(command(5:160),*) sol(1),serr(1,1),serr(1,2),err(1)
c            write(0,*) "WTF!",sol(1),serr(1,1)
        elseif(command(1:3).eq."NL1") then
            read(command(5:160),*) sol(2),serr(2,1),serr(2,2),err(2)
        elseif(command(1:3).eq."NL2") then
            read(command(5:160),*) sol(3),serr(3,1),serr(3,2),err(3)
        elseif(command(1:3).eq."NL3") then
            read(command(5:160),*) sol(4),serr(4,1),serr(4,2),err(4)
        elseif(command(1:3).eq."NL4") then
            read(command(5:160),*) sol(5),serr(5,1),serr(5,2),err(5)
        elseif(command(1:3).eq."DIL") then
            read(command(5:160),*) sol(6),serr(6,1),serr(6,2),err(6)
        elseif(command(1:3).eq."VOF") then
            read(command(5:160),*) sol(7),serr(7,1),serr(7,2),err(7)
        elseif(command(1:3).eq."ZPT") then
            read(command(5:160),*) sol(8),serr(8,1),serr(8,2),err(8)
        elseif(command(1:2).eq."EP") then
            read(command(3:3),*,err=901) np
            if(np.gt.nplanet)nplanet=np
C           j=planet parameters*(np-1)+8 initial parameters
            j=10*(np-1)+8+1
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."PE") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+2
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."BB") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+3
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."RD") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+4
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."EC") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+5
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."ES") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+6
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."KR") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+7
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."TE") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+8
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."EL") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+9
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."AL") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+10
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        endif
 501    format(A5,5(1X,1PE17.10))
c        write(0,*) command
      goto 10 !loop back to read statement
 11   continue !end loop when EOF is reached
c      write(0,*) "----------------------------"

      goto 999
 901  write(0,*) "Error on line ",i+1
      pause
      goto 999       
 999  return
      end 

 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine sigup(npt,time,mag,merr,sig,rmavg)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nmax,ntmp,i,nclip,niter,nitmax
      parameter (nmax=2000000,nitmax=50)
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
