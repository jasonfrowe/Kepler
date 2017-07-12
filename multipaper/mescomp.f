      program mescomp
      implicit none
      integer nmax,nunit,i,npt,j,k,nkoi,nbin,nbinmax,kidq,nqs,nsp
      parameter(nmax=500000,nbinmax=500)
      integer kid(nmax),flag,koikid(nmax),np(nmax),koi,nq(nmax),
     .   kidsp(nmax)
      real rp(nmax),bdatax(nbinmax),bdatay(nbinmax),errs(6),ave,std,
     .   bdatay2(nbinmax),dumr,teff(nmax),logg(nmax),teffsp(nmax),
     .   loggsp(nmax),loggcut
      double precision cdpp(nmax),work(nmax),koicdpp(nmax),med
      character*80 filename,filenameKOI,title,dumc,filenameQ,filenameSP

      filename='cdppq1q103h.dat'
      filenameKOI="koi_multi.20120926.dat"
      filenameQ="kidq1q8.txt"
      filenameSP="ktc_consol_v2.sp.dat"

      loggcut=3.5

      do 5 i=1,nmax
         logg(i)=0.0
         teff(i)=0.0
 5    continue

      nunit=10
      open(unit=nunit,file=filenameSP,status='old',err=904)
      read(nunit,*) dumc
      i=1
 21   read(nunit,*,end=22) kidsp(i),teffsp(i),dumr,loggsp(i)
         if(loggsp(i).eq.0.0) loggsp(i)=4.5
         i=i+1
      goto 21
 22   continue
      nsp=i-1
      close(nunit)
      write(0,*) "nsp: ",nsp


      open(unit=nunit,file=filename,status='old',err=901)
      i=1
 10   read(nunit,*,end=11) kid(i),cdpp(i)
c         write(0,*) kid(i)
         nq(i)=0

         do 23 j=1,nsp
            if(kid(i).eq.kidsp(j))then
               teff(i)=teffsp(j)
               logg(i)=loggsp(j)
            endif
 23      continue
c         write(0,*) kid(i),teff(i),logg(i)

         i=i+1
      goto 10
 11   continue
      npt=i-1  !is this correct? - depends on flag value
      close(nunit)
      write(0,*) "npt: ",npt

      open(unit=nunit,file=filenameQ,status='old',err=903)

 17   read(nunit,*,end=18) kidq,nqs
         flag=0
         i=1
         do while(flag.eq.0)
            if(kid(i).eq.kidq)then
               flag=1
               nq(i)=nqs
c               write(0,*) kid(i),kidq,nq(i)
            elseif(i.ge.nmax)then
               flag=1
               nq(i)=0
            else
               i=i+1
            endif
         enddo
      goto 17
 18   continue

      close(nunit)

      open(unit=nunit,file=filenameKOI,status='old',err=902)
      read(nunit,*,err=902) dumc
      read(nunit,*,err=902) dumc
      i=1
 12   read(nunit,*,end=13) koi,koikid(i),dumr,dumr,dumr,
     .   dumr,dumr,dumr,dumr,dumr,dumr,dumr,dumr,
     .   dumr,dumr,dumr,dumr,dumr,dumr,dumr,dumr,
     .   dumr,dumr,flag
         if(koi.gt.3059) goto 12
         if(flag.ge.1) goto 12
         if((koikid(i).ne.koikid(i-1)).and.(i.gt.1))then
c            write(6,*) i,koikid(i),koi
            i=i+1
         endif
         if(i.eq.1) i=i+1

      goto 12
 13   continue
      nkoi=i-1
      close(nunit)
      write(0,*) "nkoi: ",nkoi

      call pgopen('?') !open PGPlot device
      call pgask(.true.) !don't ask for new page.. just do it.
      call PGPAP ( 8.0 ,1.0) !paper size
      call pgsubp(1,1)


      nbin=30
      title='CDPP'
      call histogram(npt,rp,cdpp,work,np,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs,1)

      k=0
      do 14 i=1,npt
         do 15 j=1,nkoi
            if(koikid(j).eq.kid(i))then
               k=k+1
               koicdpp(k)=cdpp(i)
            endif
 15      continue
 14   continue
      write(0,*) "k:",k

C     Find median
      call rqsort(k,koicdpp,np) !changed npt to k
      i=k/2
      if(i.le.0) i=1
      med=koicdpp(np(i))
      write(6,*) "KOI Median CDPP: ",med,i

      j=0
      do 20 i=1,npt
         if((nq(i).ge.5).and.(cdpp(i).le.med))then
            j=j+1
          endif
 20   continue
      write(6,*) "KIDs with CDDP below median with at least 5Qs: ",j

      call histogram(k,rp,koicdpp,work,np,nbin,nbinmax,bdatax,bdatay2,
     .  title,ave,std,errs,2)


CCCCCC Working with logg cuts
      k=0
      do 24 i=1,npt
         do 25 j=1,nkoi
            if((koikid(j).eq.kid(i)).and.(logg(i).gt.loggcut))then
               k=k+1
               koicdpp(k)=cdpp(i)
            endif
 25      continue
 24   continue
      write(0,*) "k:",k

C     Find median
      call rqsort(k,koicdpp,np) !changed npt to k
      i=k/2
      if(i.le.0) i=1
      med=koicdpp(np(i))
      write(6,*) "KOI (logg) Median CDPP: ",med,i

      j=0
      do 26 i=1,npt
         if((nq(i).ge.5).and.(cdpp(i).le.med).and.
     .    (logg(i).gt.loggcut))then
            j=j+1
          endif
 26   continue
      write(6,*) "KIDs with CDDP below median, at least 5Qs, loggcut: "
     .   ,j




      call pgclos()

c      do 16 i=1,nbin
c         write(6,*) bdatay(i),bdatay2(i)
c 16   continue

      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 902  write(0,*) "Cannot open ",filenameKOI
      goto 999
 903  write(0,*) "Cannot open ",filenameQ
      goto 999
 904  write(0,*) "Cannot open ",filenameSP
      goto 999
 999  end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine histogram(npt,rp,dp,dp2,np,nbin,nbinmax,bdatax,bdatay,
     .  title,rmed,std,errs,icol)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i,j,npt2,nbin,nbinmax,icol
      integer np(npt)
      real rp(npt),datamin,datamax,bdatax(nbinmax),bdatay(nbinmax),bmax,
     .  rave,errs(6),std,rmed
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

C     Find datarange
      datamin=rp(1)
      datamax=rp(1)
      do 12 i=2,npt2
        datamin=min(rp(i),datamin)
        datamax=max(rp(i),datamax)
 12   continue
      write(0,*)datamin,datamax
c      datamin=-218.9510!-292.9580
c      datamax=300.0!199.042


      call bindata(nbin,npt2,rp,bdatax,bdatay,datamin,datamax,bmax)

c      call pgpage() !fresh plotting surface
      call pgslw(3)
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
