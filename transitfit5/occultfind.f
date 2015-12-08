      program occultfind
      implicit none
      integer nunit,nmax,npt,nfit,nplanetmax,nplanet,i,koi
      parameter(nfit=108,nmax=2000000,nplanetmax=10)
      integer ntt(nplanetmax),tflag(nmax),p(nmax)
      double precision time(nmax),flux(nmax),ferr(nmax),itime(nmax),
     .   Keplertime,sol(nfit),serr(nfit,2),err(nfit),phase(nmax),
     .   tobs(nplanetmax,nmax),omc(nplanetmax,nmax),tcor(nmax),
     .   pts(nmax),maxocc,maxoccsig,maxoccph,omean(nmax),omeanerr(nmax),
     .   x(nmax),y(nmax),dstd
      character*80 obsfile,inputsol,ttfile,cline

      if(iargc().lt.2) goto 901

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

      call getarg(2,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
      write(0,*) "reading in input solution"
      call getfitpars(nunit,nfit,nplanet,sol,serr,err)
      write(0,*) "done reading input solution"
      close(nunit) !release unit number as we are done with file

      koi=0
      if(iargc().ge.3)then
         call getarg(3,cline)
         read(cline,*) koi
      endif


      do 21 i=1,nplanet
        if(iargc().ge.3+i)then
            call getarg(3+i,ttfile)

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

C     Note, transitdur is found in marktransit5.f
      do 23 i=1,nplanet
c      do 23 i=1,1
         call findoccult(i,nplanetmax,npt,nmax,phase,time,flux,ferr,
     .      tflag,nfit,sol,ntt,tobs,omc,tcor,p,pts,maxocc,maxoccsig,
     .      maxoccph,omean,omeanerr,x,y,dstd)
         write(0,500) koi,".0",i,-maxocc*1.0d6,-maxoccsig,maxoccph,dstd
 500     format(I4,A2,I1,1X,F8.2,1X,F5.1,1X,F5.3,1X,F5.1)
c         read(5,*)
 23   continue

      goto 999
 901  write(0,*) "Usage: occultfind klc.dat n1.dat"
      goto 999
 902  write(0,*) "Error opening ",inputsol
      goto 999
 903  write(0,*) "Error opening ",obsfile
      goto 999
 905  write(0,*) "Error opening ",ttfile
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine findoccult(np,nplanetmax,npt,nmax,phase,time,flux,ferr,
     .   tflag,nfit,sol,ntt,tobs,omc,tcor,p,pts,maxocc,maxoccsig,
     .   maxoccph,omean,omeanerr,x,y,dstd)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer np,npt,nplanetmax,nmax,tflag(npt),nfit,ntt(nplanetmax),
     .   col,i,p(nmax),j,k,ii,jj,kk,bins,minpt
      double precision phase(npt),time(npt),sol(nfit),
     .   tobs(nplanetmax,nmax),omc(nplanetmax,nmax),tdur,transitdur,
     .   tcor(nmax),epo,period,toff,ttcor,flux(npt),ferr(npt),stdevf,
     .   std,mean,sigcut,ph,phdif,omean(nmax),omeanerr(nmax),pts(nmax),
     .   stdev,maxocc,maxoccsig,maxoccph,diff,differr,x(nmax),y(nmax),
     .   dmean,dstd,ph1,ph2

C     sigma clipping factor
      sigcut=4.0d0

C     Get transit-duration
      tdur=transitdur(np,nfit,sol)/86400.0d0+0.03

      col=10*(np-1)
      epo=sol(9+col)
      period=sol(10+col)

C     Add in transit-timing corrections..
      do 24 i=1,npt
         call lininterp(tobs,omc,nplanetmax,nmax,np,ntt,time(i),ttcor)
         ttcor=0.0d0
         tcor(i)=time(i)-ttcor
c         write(0,*) tcor(i),ttcor
 24   continue

      toff=-(epo/period-int(epo/period))
      if(toff.lt.0.0)toff=toff+1.0
      call phasept(npt,tcor,phase,period,toff)

C     Compute sigma clipping
      std=stdevf(npt,flux,mean,tflag)

      j=0
      do 15 i=1,npt
         if((tflag(i).eq.0).and.
     .    (abs(flux(i)-mean).lt.sigcut*std))then
            j=j+1
         endif
 15   continue
      minpt=int(dble(j)*tdur/period/2.0) !should include a good sample
c      write(0,*) minpt

C     Sort based on phase
      call rqsort(npt,phase,p)

      ph=tdur/period/2.0d0 !phase covered by 1/2 transitduration
      do 10 i=1,npt
         j=1
         k=2
         phdif=phase(p(i))-phase(p(j))
         do while(phdif.gt.ph)
            if(j.lt.i)then
               j=j+1
               phdif=phase(p(i))-phase(p(j))
            else
               phdif=ph
            endif
         enddo

         phdif=phase(p(k))-phase(p(i))
         do while(phdif.lt.ph)
            if(k.lt.npt)then
               k=k+1
               phdif=phase(p(k))-phase(p(i))
            else
               phdif=ph
            endif
         enddo

         omean(p(i))=0.0d0
         jj=0
         do 11 ii=j,k
            if((tflag(p(ii)).eq.0).and.
     .       (abs(flux(p(ii))-mean).lt.sigcut*std))then
               jj=jj+1
               omean(p(i))=omean(p(i))+flux(p(ii))
               pts(jj)=flux(p(ii))
            endif
 11      continue
         if(jj.gt.minpt)then
c            omean=omean/dble(jj)
            omeanerr(p(i))=stdev(jj,pts,omean(p(i)))/sqrt(dble(jj))
         else
            omean(p(i))=0.0d0
            omeanerr(p(i))=0.0d0
         endif
 10   continue

c      boxbin=3.0*tdur/period
c      call boxcar(npt2,x,y,z,boxbin,tf)


      kk=0 !counter
      ph=tdur/period*1.0 !phase covered by transitduration
c      write(0,*) "ph:",ph
      maxoccsig=0.0d0 !return maximum occultation found.
      do 12 i=1,npt
         j=1
         k=2
         phdif=phase(p(i))-phase(p(j))
         do while(phdif.gt.ph)
            if(j.lt.i)then
               j=j+1
               phdif=phase(p(i))-phase(p(j))
            else
               phdif=ph
            endif
         enddo

         phdif=phase(p(k))-phase(p(i))
         do while(phdif.lt.ph)
            if(k.lt.npt)then
               k=k+1
               phdif=phase(p(k))-phase(p(i))
            else
               phdif=ph
            endif
         enddo

         if(omean(p(i)).gt.0.0)then
            if((omean(p(j)).gt.0.0).and.(omean(p(k)).gt.0.0))then
               diff=omean(p(i))-(omean(p(j))+omean(p(k)))/2.0d0
               differr=sqrt(omeanerr(p(i))**2.0+
     .           (omeanerr(p(j))/2.0d0)**2.0+
     .           (omeanerr(p(k))/2.0d0)**2.0)
            elseif((omean(p(j)).eq.0.0).and.(omean(p(k)).gt.0.0))then
               diff=omean(p(i))-omean(p(k))
               differr=sqrt(omeanerr(p(i))**2.0+omeanerr(p(k))**2.0)
            elseif((omean(p(j)).gt.0.0).and.(omean(p(k)).eq.0.0))then
               diff=omean(p(i))-omean(p(j))
               differr=sqrt(omeanerr(p(i))**2.0+omeanerr(p(j))**2.0)
            else
               diff=0.0
               differr=1.0
            endif
            kk=kk+1
            x(kk)=phase(p(i))
            y(kk)=diff/differr

            write(6,*) phase(p(i)),diff,differr
c            write(6,*) phase(p(i)),omean(p(i)),omeanerr(p(i))
            if(diff/differr.lt.maxoccsig)then
               maxoccsig=diff/differr
               maxocc=diff
               maxoccph=phase(p(i))
            endif
         endif

 12   continue
c      write(0,501) "Planet#",np," done"
 501  format(A7,I2,A5)

      j=0
      ph=tdur/period*1.0
      ph1=maxoccph-ph
      ph2=maxoccph+ph
      do 14 i=1,kk
         if((x(i).lt.ph1).or.(x(i).gt.ph2))then
            j=j+1
            y(j)=y(i)
         endif
 14   continue
      dstd=stdev(j,y,dmean)
c      write(0,*) dstd
c      maxoccsig=maxoccsig/dstd

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine boxcar(npt,time,mag,merr,boxbin,tflag)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Robust Boxcar filter for unevenly sampled data.
C
C     This routine assumes that the data are sorted in time from
C     earliest to latest.
      implicit none
      integer npt,nmax,i,j,k,ii,dcsum,dcsumold,it,itmax
      parameter(nmax=2000000)
      integer dc(nmax),nd(nmax),tflag(nmax)
      double precision time(nmax),mag(nmax),merr(nmax),boxbin,diffm,
     .  diffp,bbd2,avgm(nmax),avge(nmax),abin(nmax),std,adif(nmax),mean,
     .  stdev,sigcut,pts(nmax)

      sigcut=3.0
      itmax=10
      it=0

      do 3 i=1,npt
        dc(i)=tflag(i)
 3    continue
      dcsum=0

 4    continue
      dcsumold=dcsum
C     Initialize arrays to zero.
      do 5 i=1,npt
        avgm(i)=0.0d0
        avge(i)=0.0d0
        abin(i)=0.0d0
c        dc=0
 5    continue

      bbd2=boxbin/2.0d0 !half of boxcar width
      i=1 !initialize arrays that determine which elements are inside
      j=1 !the boxcar width centered on time(k)
      do 10 k=1,npt
        diffp=time(j)-time(k)
        do while (diffp.lt.bbd2)
            if(j.ge.npt) then
                diffp=bbd2
                j=npt
            else
                j=j+1
                diffp=time(j)-time(k)
            endif
        enddo
        diffm=time(k)-time(i)
        do while (diffm.gt.bbd2)
            if(i.ge.k) then
                diffm=bbd2
                i=k
            else
                i=i+1
                diffm=time(k)-time(i)
            endif
        enddo
        if(diffp.gt.bbd2)j=j-1
        if(diffm.gt.bbd2)i=i-1
c        write(6,*) k,i,j,time(k)-time(i),time(j)-time(k)
c        read(5,*)
        do 11 ii=i,j
            if(dc(ii).eq.0)then
                avgm(k)=avgm(k)+mag(ii)/merr(ii)
                avge(k)=avge(k)+1.0
                abin(k)=abin(k)+1.0/merr(ii)
            endif
 11     continue
        if(abin(k).eq.0)then
            do 16 ii=i,j
                avgm(k)=avgm(k)+mag(ii)/merr(ii)
                avge(k)=avge(k)+1.0
                abin(k)=abin(k)+1.0/merr(ii)
 16         continue
        endif
        avgm(k)=avgm(k)/abin(k)
C        Find median
         call rqsort(j,mag,nd)
         avge(k)=mag(nd(j/2))
c        avge(k)=(avge(k)**0.5)/abin(k)
 10   continue

      k=0
      do 14 i=1,npt
        if(dc(i).eq.0)then
            k=k+1
            adif(i)=abs(mag(i)-avgm(i))
            pts(k)=adif(i)
        endif
 14   continue
      std=stdev(k,pts,mean)
      dcsum=0
      do 15 i=1,npt
        if((adif(i).gt.mean+sigcut*std).or.(tflag(i).eq.1)) then
            dc(i)=1
c            write(11,*) time(i),mag(i)
        else
            dc(i)=0
        endif
c        write(0,*) "adif:",adif(i),sigcut*std
c        read(5,*)
        dcsum=dcsum+dc(i)
 15   continue
c      write(0,*) dcsum,dcsumold
      it=it+1
      if((dcsum.ne.dcsumold).and.(it.le.itmax)) goto 4


c      i=1
c      diffm=time(i)-time(1)
c      do while(diffm.lt.bbd2)
c        i=i+1
c        diffm=time(i)-time(1)
c      enddo
c      j=npt
c      diffp=time(npt)-time(j)
c      do while(diffp.lt.bbd2)
c        j=j-1
c        diffp=time(npt)-time(j)
c      enddo
c
c      npt=0
c      do 12 k=i,j
c        npt=npt+1
c        time(npt)=time(k)
c        mag(npt)=mag(k)-avgm(k)
c        merr(npt)=sqrt(merr(k)*merr(k)+avge(k)*avge(k))
c 12   continue
      do 13 k=1,npt
        mag(k)=mag(k)-avgm(k)
        merr(k)=sqrt(merr(k)*merr(k)+avge(k)*avge(k))
 13   continue

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function stdev(npt,pts,mean)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Calculates standard deviation of data set given the mean.
      implicit none

      integer npt,i
      double precision pts(npt),mean,s,ep,adev,p,var,sdev,mean2

      s=0.
      do 11 i=1,npt
         s=s+pts(i)
 11   continue
      mean=s/npt

      ep=0.
      var=0.
      do 10 i=1,npt
         s=pts(i)-mean
         ep=ep+s
         p=s*s
         var=var+p
 10   continue
      var=(var-ep**2/npt)/(npt-1)
      sdev=sqrt(var)

      stdev=sdev

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function stdevf(npt,pts,mean,tflag)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Calculates standard deviation of data set given the mean.
      implicit none

      integer npt,i,tflag(npt),j
      double precision pts(npt),mean,s,ep,adev,p,var,sdev,mean2

      s=0.
      j=0
      do 11 i=1,npt
         if(tflag(i).eq.0)then
            s=s+pts(i)
            j=j+1
         endif
 11   continue
      mean=s/dble(j)

      ep=0.
      var=0.
      j=0
      do 10 i=1,npt
         if(tflag(i).eq.0)then
            s=pts(i)-mean
            ep=ep+s
            p=s*s
            var=var+p
            j=j+1
         endif
 10   continue
      var=(var-ep**2/dble(j))/dble(j-1)
      sdev=sqrt(var)

      stdevf=sdev

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
