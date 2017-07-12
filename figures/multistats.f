      program multistats
      implicit none
      integer nmax,i,nunit,j,k,nall,nTeq,nbin,nbinmax,nEB,nflag,kidp,
     .   nplanetmax,koip
      parameter(nmax=10000,nbinmax=500,nplanetmax=10)
      integer npt,kid(nmax),pflag(nmax),flag(nmax),cf(nmax),sflag(nmax),
     .   np(nmax),koimax,koinum,n(nmax),np2(nmax),kideb(nmax),dumi,
     .   koimat(nmax,nplanetmax)
      real realp(nmax),bdatax(nbinmax),bdatay(nbinmax),errs(6),ave,std,
     .   datamin2,datamax2
      real*8 koi(nmax),kepmag(nmax),Rp(nmax),T0(nmax),Per(nmax),
     .   Teff(nmax),logg(nmax),Rstar(nmax),FeH(nmax),adrs(nmax),
     .   rdr(nmax),b(nmax),asemi(nmax),Teq(nmax),gmag(nmax),RA(nmax),
     .   DEC(nmax),tdepth(nmax),SN(nmax),Tdur(nmax),chisq(nmax),
     .   dp(nmax),work(nmax),dp2(nmax),PEB(nmax),dumr
      character*80 filename,dumc,header,title,filenameEB

      nunit=10
      filename="koi_multi.20120926.dat"

      i=0
      open(unit=nunit,file=filename,status='old',err=901)
      read(nunit,*,err=902) dumc
      read(nunit,*,err=902) dumc
      i=1 !count number of lines read
 10   read(nunit,*,end=11,err=15) koi(i),kid(i),kepmag(i),Rp(i),T0(i),
     .   Per(i),Teff(i),logg(i),Rstar(i),FeH(i),pflag(i),adrs(i),rdr(i),
     .   b(i),asemi(i),Teq(i),gmag(i),RA(i),DEC(i),Tdepth(i),SN(i),
     .   Tdur(i),chisq(i),flag(i)!,cf(i),sflag(i)
         if(koi(i).gt.3059.0) goto 10
c         if(cf(i).eq.2) goto 10
         i=i+1
      goto 10
 15   write(0,*) "Error on line: ",i+2
      goto 10
 11   continue
      close(nunit)
      npt=i-1
      write(0,*) "Number of lines read: ",npt

      i=1
      filenameEB="EB2.list"
      open(unit=nunit,file=filenameEB,status='old',err=903)
 34   read(nunit,*,end=35) kidEB(i),PEB(i)
         nflag=0
         do 36 j=1,npt
            if(kidEB(i).eq.kid(j))then
               nflag=1
            endif
 36      continue
         if(nflag.eq.0) i=i+1
      goto 34
 35   continue
      close(nunit)
      nEB=i
      write(0,*) "Number of EBs: ",nEB
      do 45 i=1,nEB
         if(kidEB(i).lt.10)then
            write(6,505) "0000000",kidEB(i),"||TPS_EXCLUDE|||"
 505        format(A7,I1,A16)
         elseif(kidEB(i).lt.100)then
            write(6,506) "000000",kidEB(i),"||TPS_EXCLUDE|||"
 506        format(A6,I2,A16)
         elseif(kidEB(i).lt.1000)then
            write(6,507) "00000",kidEB(i),"||TPS_EXCLUDE|||"
 507        format(A5,I3,A16)
         elseif(kidEB(i).lt.10000)then
            write(6,508) "0000",kidEB(i),"||TPS_EXCLUDE|||"
 508        format(A4,I4,A16)
         elseif(kidEB(i).lt.100000)then
            write(6,509) "000",kidEB(i),"||TPS_EXCLUDE|||"
 509        format(A3,I5,A16)
         elseif(kidEB(i).lt.1000000)then
            write(6,510) "00",kidEB(i),"||TPS_EXCLUDE|||"
 510        format(A2,I6,A16)
         elseif(kidEB(i).lt.10000000)then
            write(6,511) "0",kidEB(i),"||TPS_EXCLUDE|||"
 511        format(A1,I7,A16)
         else
            write(6,512) kidEB(i),"||TPS_EXCLUDE|||"
 512        format(I8,A16)
         endif
 45   continue

      i=1
      filenameEB="eb2.list"
      open(unit=nunit,file=filenameEB,status='old',err=903)
 37   read(nunit,*,end=38) dumi,dumr
         nflag=0
         do 39 j=1,npt
            if(dumi.eq.kid(j))then
               nflag=1
            endif
 39      continue
         do 40 j=1,nEB
            if(dumi.eq.kidEB(j))then
               nflag=1
            endif
 40      continue
         if(nflag.eq.0) then
            i=i+1
            nEB=nEB+1
            kidEB(nEB)=dumi
            PEB(nEB)=dumr
         endif
      goto 37
 38   continue
      close(nunit)
      write(0,*) "EBs: ",nEB,i

      j=0
      do 41 i=1,nEB
         if(PEB(i).gt.0.0)then
            j=j+1
            kidEB(j)=kidEB(i)
            PEB(j)=PEB(i)
         endif
 41   continue
      nEB=j
      write(0,*) "EB (per cut)",nEB

c      do 42 i=1,nEB
c         write(6,*) kidEB(i),PEB(i)
c 42   continue

      call pgopen('?') !open PGPlot device
      call pgask(.true.) !don't ask for new page.. just do it.
      call PGPAP ( 8.0 ,1.0) !paper size
      call pgsubp(2,6)

      header="------------------------------------------------------"
      write(6,501) header
      header="Comment                     1    2    3    4    5    6"
      write(6,501) header
      header="------------------------------------------------------"
      write(6,501) header
 501  format(A60)
CCCCCCCCCCCCCCCCCCCCCCCCCCc
C     count up number of planets in each KOI systems
      koimax=0
      do 13 i=1,nmax
         np(i)=0
         n(i)=0
 13   continue

      do 12 i=1,npt
         koinum=int(koi(i)+0.1)
         koimax=max(koinum,koimax) !largest KOI number
         np(koinum)=np(koinum)+1
 12   continue

      header="No Cuts:"
      do 14 i=1,koimax
         if(np(i).gt.0)then
            n(np(i))=n(np(i))+1
         endif
 14   continue
      write(6,500) header,(n(i),i=1,6)
 500  format(A25,6(I4,1X))

      j=0
      do 28 i=1,npt
         koinum=int(koi(i)+0.1)
         if(np(koinum).ge.2)then
            j=j+1
            dp(j)=b(i)
            dp2(j)=rp(i)
         endif
 28   continue
      nbin=30
      title='b (No Cuts)'
      call pgpanl(1,1)
      datamin2=0.0
      datamax2=1.5
      call histogram(j,realp,dp,work,np2,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs,1,datamin2,datamax2)
      title='Rp (No Cuts)'
      call pgpanl(2,1)
      datamin2=0.0
      datamax2=20.0
      call histogram(j,realp,dp2,work,np2,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs,1,datamin2,datamax2)

      j=0
      do 29 i=1,npt
         koinum=int(koi(i)+0.1)
         if((np(koinum).ge.2).and.(flag(i).eq.1))then
            j=j+1
            dp(j)=b(i)
            dp2(j)=rp(i)
         endif
 29   continue
      nbin=30
      title='b cut (FPs)'
      call pgpanl(1,2)
      datamin2=0.0
      datamax2=1.5
      call histogram(j,realp,dp,work,np2,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs,1,datamin2,datamax2)
      title='Rp cut (FPs)'
      call pgpanl(2,2)
      datamin2=0.0
      datamax2=20.0
      call histogram(j,realp,dp2,work,np2,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs,1,datamin2,datamax2)

      j=0
      do 30 i=1,npt
         koinum=int(koi(i)+0.1)
         if((np(koinum).ge.2).and.((flag(i).eq.1).or.(sn(i).le.10)))then
            j=j+1
            dp(j)=b(i)
            dp2(j)=rp(i)
         endif
 30   continue
      nbin=30
      title='b cut (FPs and S/N)'
      call pgpanl(1,3)
      datamin2=0.0
      datamax2=1.5
      call histogram(j,realp,dp,work,np2,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs,1,datamin2,datamax2)
      title='Rp cut (FPs and S/N)'
      call pgpanl(2,3)
      datamin2=0.0
      datamax2=20.0
      call histogram(j,realp,dp2,work,np2,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs,1,datamin2,datamax2)

      j=0
      do 31 i=1,npt
         koinum=int(koi(i)+0.1)
         if((np(koinum).ge.2).and.((flag(i).eq.1).or.(sn(i).le.10)
     .     .or.(per(i).lt.1.0)))then
            j=j+1
            dp(j)=b(i)
            dp2(j)=rp(i)
         endif
 31   continue
      nbin=30
      title='b cut (FPs, S/N and P)'
      call pgpanl(1,4)
      datamin2=0.0
      datamax2=1.5
      call histogram(j,realp,dp,work,np2,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs,1,datamin2,datamax2)
      title='Rp cut (FPs, S/N and P)'
      call pgpanl(2,4)
      datamin2=0.0
      datamax2=20.0
      call histogram(j,realp,dp2,work,np2,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs,1,datamin2,datamax2)

      j=0
      do 32 i=1,npt
         koinum=int(koi(i)+0.1)
         if((np(koinum).ge.2).and.((flag(i).eq.1).or.(sn(i).le.10)
     .     .or.(per(i).lt.1.0).or.(rdr(i)+b(i).ge.0.98)))then
            j=j+1
            dp(j)=b(i)
            dp2(j)=rp(i)
         endif
 32   continue
      nbin=30
      title='b cut (FPs, S/N, P and b)'
      call pgpanl(1,5)
      datamin2=0.0
      datamax2=1.5
      call histogram(j,realp,dp,work,np2,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs,1,datamin2,datamax2)
      title='Rp cut (FPs, S/N, P and b)'
      call pgpanl(2,5)
      datamin2=0.0
      datamax2=20.0
      call histogram(j,realp,dp2,work,np2,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs,1,datamin2,datamax2)

ccccccccccccccccccccccccccc
C     count up number of planets in each KOI systems
      koimax=0
      do 16 i=1,nmax
         np(i)=0
         n(i)=0
 16   continue

      do 17 i=1,npt
         if(flag(i).lt.1)then
            koinum=int(koi(i)+0.1)
            koimax=max(koinum,koimax) !largest KOI number
            np(koinum)=np(koinum)+1
         endif
c         write(6,*) koi(i),np(koinum),j
 17   continue

      header="Cut FPs:"
      do 18 i=1,koimax
         if(np(i).gt.0)then
            n(np(i))=n(np(i))+1
         endif
 18   continue
      write(6,500) header,(n(i),i=1,6)

ccccccccccccccccccccccccccc

ccccccccccccccccccccccccccc
C     count up number of planets in each KOI systems
      koimax=0
      do 19 i=1,nmax
         np(i)=0
         n(i)=0
 19   continue
      do 20 i=1,npt
         if((flag(i).lt.1).and.(sn(i).gt.10.0))then
            koinum=int(koi(i)+0.1)
            koimax=max(koinum,koimax) !largest KOI number
            np(koinum)=np(koinum)+1
         endif
c         write(6,*) koi(i),np(koinum),j
 20   continue

      header="Cut FPs and S/N:"
      do 21 i=1,koimax
         if(np(i).gt.0)then
            n(np(i))=n(np(i))+1
         endif
 21   continue
      write(6,500) header,(n(i),i=1,6)


CCCCCCCCCCCCCCCCCCCCCCCCCCcc
C     count up number of planets in each KOI systems
      koimax=0
      do 22 i=1,nmax
         np(i)=0
         n(i)=0
 22   continue
      do 23 i=1,npt
         if((flag(i).lt.1).and.(sn(i).gt.10.0).and.(per(i).gt.1.0))then
            koinum=int(koi(i)+0.1)
            koimax=max(koinum,koimax) !largest KOI number
            np(koinum)=np(koinum)+1
         endif
c         write(6,*) koi(i),np(koinum),j
 23   continue

      header="Cut FPs, S/N and P:"
      do 24 i=1,koimax
         if(np(i).gt.0)then
            n(np(i))=n(np(i))+1
         endif
 24   continue
      write(6,500) header,(n(i),i=1,6)

CCCCCCCCCCCCCCCCCCCCCCCCCCcc

      do 43 i=1,nmax
         do 44 j=1,nplanetmax
            koimat(i,j)=0
 44      continue
 43   continue

C     count up number of planets in each KOI systems
      koimax=0
      do 25 i=1,nmax
         np(i)=0
         n(i)=0
 25   continue
      do 26 i=1,npt
c         if(flag(i).lt.1)then
         if((flag(i).lt.1).and.(sn(i).gt.10.0).and.(per(i).gt.1.0)
     .      .and.(rdr(i)+b(i).lt.0.98))then
            koinum=int(koi(i)+0.1)
            koimax=max(koinum,koimax) !largest KOI number
            np(koinum)=np(koinum)+1
            koip=int(100*(koi(i)-koinum+0.001))
            if(koip.eq.0) koip=1
            koimat(koinum,koip)=1
c            write(6,*) koi(i),kid(i)
         endif
c         write(6,*) koi(i),np(koinum),j
 26   continue

      header="Cut FPs, S/N, P and b:"
      do 27 i=1,koimax
         if(np(i).gt.0)then
            n(np(i))=n(np(i))+1

C           Trying to output Good KOI systems for modelling
            do 42 j=1,npt
               if(i.eq.int(koi(j)+0.1))then
                  kidp=kid(j)
c            if((np(i).ge.1).and.flag(i).lt.1) write(6,*) kid(j),per(j)
               endif
 42         continue
c            if(np(i).gt.1)write(6,502) i,kidp,
c     .       (koimat(i,j),j=1,nplanetmax)
 502        format(I4,1X,I11,10(1X,I1))
         endif
 27   continue
      write(6,500) header,(n(i),i=1,6)

      header="------------------------------------------------------"
      write(6,501) header

      j=0
      do 33 i=1,npt
         koinum=int(koi(i)+0.1)
         if((np(koinum).ge.2).and.(flag(i).lt.1).and.(sn(i).gt.10)
     .     .and.(per(i).gt.1.0).and.(rdr(i)+b(i).lt.0.98))then
            j=j+1
            dp(j)=b(i)
            dp2(j)=rp(i)
         endif
 33   continue
      nbin=30
      title='b (All Cuts and still multi)'
      call pgpanl(1,6)
      datamin2=0.0
      datamax2=1.5
      call histogram(j,realp,dp,work,np2,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs,1,datamin2,datamax2)
      title='Rp (All Cuts and still multi)'
      call pgpanl(2,6)
      datamin2=0.0
      datamax2=20.0
      call histogram(j,realp,dp2,work,np2,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs,1,datamin2,datamax2)


      call pgclos()

      goto 999
901   write(0,*) "Cannot open: ",filename
      goto 999
902   write(0,*) "Error reading on line: ",i+2
      goto 999
903   write(0,*) "Cannot open: ",filenameEB
      goto 999
999   end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine histogram(npt,rp,dp,dp2,np,nbin,nbinmax,bdatax,bdatay,
     .  title,rmed,std,errs,icol,datamin2,datamax2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i,j,npt2,nbin,nbinmax,icol
      integer np(npt)
      real rp(npt),datamin,datamax,bdatax(nbinmax),bdatay(nbinmax),bmax,
     .  rave,errs(6),std,rmed,datamin2,datamax2
      double precision dp(npt),dp2(npt),ave,var,med
      character*80 title

C     First we remove the average (helps a lot with double -> real)
      call avevar(dp,npt,ave,var)
      rave=real(ave) !convert to real*4
      std=real(sqrt(var))

C     Now we convert dp to rp
      j=0 !because we have sigma-clipping, we need a counter.
      do 10 i=1,npt
        if((abs(dp(i)-ave).lt.4.0*sqrt(var)).and.(dp(i).lt.1500.0))then
            j=j+1
            dp2(j)=dp(i)-ave
            rp(j)=real(dp2(j))
        endif
 10   continue
      npt2=j

C     Find median
      call rqsort(npt2,dp2,np) !changed npt to k
      i=npt2/2
      if(i.le.0) i=1
      med=dp2(np(i))
      rmed=real(med)

cC     Find datarange
c      datamin=rp(1)
c      datamax=rp(1)
c      do 12 i=2,npt2
c        datamin=min(rp(i),datamin)
c        datamax=max(rp(i),datamax)
c 12   continue
cc      write(0,*)datamin,datamax
c      datamin=0.0-real(ave)
c      datamax=1.5-real(ave)
      datamin=datamin2-ave
      datamax=datamax2-ave


      call bindata(nbin,npt2,rp,bdatax,bdatay,datamin,datamax,bmax)

c      call pgpage() !fresh plotting surface
      call pgslw(1)
      call pgsch(2.0)
c         call windowsetup(xb1,xb2,yb1,yb2) !make a square plotting surface
c         call pgvport(xb1,xb2,yb1,yb2)
      call pgvport(0.2,1.0,0.2,0.9)
c      fc=2.*(datamax-datamin)/real(nbin) !center histograms
      call pgwindow(datamin,datamax,0.,bmax+0.1*bmax) !set size

C        Add axis labels
      call pglabel(title,"Relative Probability","")
      call pgsci(icol)
      call pgbin(nbin,bdatax,bdatay,.false.) !plot the histogram
      call pgsci(1)

C     Shift axis scale to account for average removal
      call pgwindow(datamin+rave,datamax+rave,0.,1.0+0.1*1.0)
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0) !add boarders

      rmed=rmed+rave !correct for average removal

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bindata(nbin,npt,pdata,bdatax,bdatay,datamin,datamax,
     .   bmax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nbin,npt,i,bin
      real pdata(npt),bdatax(nbin),bdatay(nbin),datamin,datamax,
     .   binsize,bmax,tsum

C     get size of bin.
      binsize=(datamax-datamin)/real(nbin-1)

C     initalize bdata
      do 20 i=1,nbin
         bdatay(i)=0.
 20   continue

C     loop over data and place data into proper bin.
      do 10 i=1,npt
C        calculate bin number
         bin=int((pdata(i)-datamin)/binsize)+1
         if((bin.gt.0).and.(bin.le.nbin)) bdatay(bin)=bdatay(bin)+1.0
 10   continue

C     get the max value in a bin and assign bin value.
      tsum=0.
      bmax=0.
      do 30 i=1,nbin
         bmax=max(bdatay(i),bmax)
         bdatax(i)=datamin+real(i-1)*binsize !added -1
         tsum=tsum+bdatay(i)
 30   continue
c      write(6,*) "Sum:",tsum

      return
      end

c**********************************************************************
      subroutine rqsortr(n,a,p)
c======================================================================
c     Return integer array p which indexes array a in increasing order.
c     Array a is not disturbed.  The Quicksort algorithm is used.
c
c     B. G. Knapp, 86/12/23
c
c     Reference: N. Wirth, Algorithms and Data Structures,
c     Prentice-Hall, 1986
c======================================================================
      implicit none

c     Input:
      integer   n
      real      a(n)

c     Output:
      integer   p(n)

c     Constants
      integer   LGN, Q
      parameter (LGN=32, Q=11)
c        (LGN = log base 2 of maximum n;
c         Q = smallest subfile to use quicksort on)

c     Local:
      real      x
      integer   stackl(LGN),stackr(LGN),s,t,l,m,r,i,j

c     Initialize the stack
      stackl(1)=1
      stackr(1)=n
      s=1

c     Initialize the pointer array
      do 1 i=1,n
         p(i)=i
    1 continue

    2 if (s.gt.0) then
         l=stackl(s)
         r=stackr(s)
         s=s-1

    3    if ((r-l).lt.Q) then

c           Use straight insertion
            do 6 i=l+1,r
               t = p(i)
               x = a(t)
               do 4 j=i-1,l,-1
                  if (a(p(j)).le.x) goto 5
                  p(j+1) = p(j)
    4          continue
               j=l-1
    5          p(j+1) = t
    6       continue
         else

c           Use quicksort, with pivot as median of a(l), a(m), a(r)
            m=(l+r)/2
            t=p(m)
            if (a(t).lt.a(p(l))) then
               p(m)=p(l)
               p(l)=t
               t=p(m)
            endif
            if (a(t).gt.a(p(r))) then
               p(m)=p(r)
               p(r)=t
               t=p(m)
               if (a(t).lt.a(p(l))) then
                  p(m)=p(l)
                  p(l)=t
                  t=p(m)
               endif
            endif

c           Partition
            x=a(t)
            i=l+1
            j=r-1
    7       if (i.le.j) then
    8          if (a(p(i)).lt.x) then
                  i=i+1
                  goto 8
               endif
    9          if (x.lt.a(p(j))) then
                  j=j-1
                  goto 9
               endif
               if (i.le.j) then
                  t=p(i)
                  p(i)=p(j)
                  p(j)=t
                  i=i+1
                  j=j-1
               endif
               goto 7
            endif

c           Stack the larger subfile
            s=s+1
            if ((j-l).gt.(r-i)) then
               stackl(s)=l
               stackr(s)=j
               l=i
            else
               stackl(s)=i
               stackr(s)=r
               r=j
            endif
            goto 3
         endif
         goto 2
      endif
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE avevar(data,n,ave,var)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n,m
      REAL*8 ave,var,data(n)
C Given array data(1:n), returns its mean as ave and its variance as var.
      INTEGER j
      REAL*8 s,ep
      ave=0.0
      m=min(1000,n)
      do 11 j=1,m
         ave=ave+data(j)
 11   continue
      ave=ave/m
      var=0.0
      ep=0.0
      do 12 j=1,n
         s=data(j)-ave
         ep=ep+s
         var=var+s*s
 12   continue
      var=(var-ep**2/n)/(n-1) !Corrected two-pass formula (14.1.8).
      return
      END

c**********************************************************************
      subroutine rqsort(n,a,p)
c======================================================================
c     Return integer array p which indexes array a in increasing order.
c     Array a is not disturbed.  The Quicksort algorithm is used.
c
c     B. G. Knapp, 86/12/23
c
c     Reference: N. Wirth, Algorithms and Data Structures,
c     Prentice-Hall, 1986
c======================================================================
      implicit none

c     Input:
      integer   n
      real*8      a(n)

c     Output:
      integer   p(n)

c     Constants
      integer   LGN, Q
      parameter (LGN=32, Q=11)
c        (LGN = log base 2 of maximum n;
c         Q = smallest subfile to use quicksort on)

c     Local:
      real*8      x
      integer   stackl(LGN),stackr(LGN),s,t,l,m,r,i,j

c     Initialize the stack
      stackl(1)=1
      stackr(1)=n
      s=1

c     Initialize the pointer array
      do 1 i=1,n
         p(i)=i
    1 continue

    2 if (s.gt.0) then
         l=stackl(s)
         r=stackr(s)
         s=s-1

    3    if ((r-l).lt.Q) then

c           Use straight insertion
            do 6 i=l+1,r
               t = p(i)
               x = a(t)
               do 4 j=i-1,l,-1
                  if (a(p(j)).le.x) goto 5
                  p(j+1) = p(j)
    4          continue
               j=l-1
    5          p(j+1) = t
    6       continue
         else

c           Use quicksort, with pivot as median of a(l), a(m), a(r)
            m=(l+r)/2
            t=p(m)
            if (a(t).lt.a(p(l))) then
               p(m)=p(l)
               p(l)=t
               t=p(m)
            endif
            if (a(t).gt.a(p(r))) then
               p(m)=p(r)
               p(r)=t
               t=p(m)
               if (a(t).lt.a(p(l))) then
                  p(m)=p(l)
                  p(l)=t
                  t=p(m)
               endif
            endif

c           Partition
            x=a(t)
            i=l+1
            j=r-1
    7       if (i.le.j) then
    8          if (a(p(i)).lt.x) then
                  i=i+1
                  goto 8
               endif
    9          if (x.lt.a(p(j))) then
                  j=j-1
                  goto 9
               endif
               if (i.le.j) then
                  t=p(i)
                  p(i)=p(j)
                  p(j)=t
                  i=i+1
                  j=j-1
               endif
               goto 7
            endif

c           Stack the larger subfile
            s=s+1
            if ((j-l).gt.(r-i)) then
               stackl(s)=l
               stackr(s)=j
               l=i
            else
               stackl(s)=i
               stackr(s)=r
               r=j
            endif
            goto 3
         endif
         goto 2
      endif
      return
      end
