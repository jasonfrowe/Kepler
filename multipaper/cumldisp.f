      program cumldisp
      implicit none
      integer nmax,nunit,i,j,dumi,nflag,nEB,koimax,koinum,ncd
      parameter(nmax=10000)
      integer npt,kid(nmax),pflag(nmax),flag(nmax),kideb(nmax),np(nmax),
     .   ns(nmax),fp,nloop
      real xp(nmax),yp(nmax),pc
      double precision koi(nmax),kepmag(nmax),Rp(nmax),T0(nmax),
     .   Per(nmax),Teff(nmax),logg(nmax),Rstar(nmax),FeH(nmax),
     .   adrs(nmax),rdr(nmax),b(nmax),asemi(nmax),Teq(nmax),gmag(nmax),
     .   RA(nmax),DEC(nmax),tdepth(nmax),SN(nmax),Tdur(nmax),
     .   chisq(nmax),PEB(nmax),dumr,koifp,ckoi,maxper,rsn,koisn
      character*80 filename,dumc,filenameEB,fpfile,cfile,koicharsn

      pc=150.0 !at what value to extract value of CD

      nunit=10
c      filename="koi_characteristics.20130509.dat"
      filename="koi_multi.20120926.dat"

      maxper=335.0 !maximum period for distributions

      i=0
      open(unit=nunit,file=filename,status='old',err=901)
      read(nunit,*,err=902) dumc
      read(nunit,*,err=902) dumc
      i=1 !count number of lines read
 10   read(nunit,*,end=11,err=15) koi(i),kid(i),kepmag(i),Rp(i),T0(i),
     .   Per(i),Teff(i),logg(i),Rstar(i),FeH(i),pflag(i),adrs(i),rdr(i),
     .   b(i),asemi(i),Teq(i),gmag(i),RA(i),DEC(i),Tdepth(i),SN(i),
     .   Tdur(i),chisq(i),flag(i)!,cf(i),sflag(i)
         if(koi(i).eq.floor(koi(i))) koi(i)=koi(i)+0.01d0
         if(koi(i).gt.3059.0) goto 10
         if(per(i).gt.maxper) goto 10
c         if(cf(i).eq.2) goto 10
c         write(0,*) koi(i)
         i=i+1
      goto 10
 15   write(0,*) "Error on line: ",i+2
      goto 10
 11   continue
      close(nunit)
      npt=i-1
      write(0,*) "Number of lines read: ",npt

      koicharsn="koi_characteristics.20130509.dat"
      open(unit=nunit,file=koicharsn,status='old')
      read(nunit,*) dumc
      read(nunit,*) dumc
 61   read(nunit,*,end=62) koisn,(dumr,j=1,19),rsn
         nloop=0
         i=0
         do while (nloop.eq.0)
            i=i+1
            if(i.gt.npt) nloop=2
            if((koisn.eq.koi(i)).and.(nloop.eq.0))then
               nloop=1
               sn(i)=rsn
            endif
         enddo
      goto 61
 62   continue
      close(nunit)


      fpfile="Fergalfile.20130801.q1q8.awk.txt"
      open(unit=nunit,file=fpfile,status='old')
 53   read(nunit,*,end=54) koifp,fp
         nloop=0
         i=0
         do while (nloop.eq.0)
            i=i+1
            if(i.gt.npt) nloop=2
            if((koifp.eq.koi(i)).and.(nloop.eq.0))then
               nloop=1
               flag(i)=fp
            endif
         enddo
      goto 53
 54   continue
      close(nunit)

      cfile="PerT0collide.20130805.txt"
      open(unit=nunit,file=cfile,status='old')
      read(nunit,*) dumc
 55   read(nunit,*,end=56) ckoi
         nloop=0
         i=0
         do while (nloop.eq.0)
            i=i+1
            if(i.gt.npt) nloop=2
            if((ckoi.eq.koi(i)).and.(nloop.eq.0))then
               nloop=1
               flag(i)=2
c               write(0,*) ckoi,koi(i),flag(i)
            endif
         enddo
         goto 55
 56   continue
      close(nunit)


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
      call pgslw(4)
      call pgsch(2.0)
      call pgvport(0.2,0.9,0.2,0.9)
      call pgwindow(0.0,8.0,0.0,0.8) !set size
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel ("Period","Cumulative Distribution","")

      koimax=0
      do 13 i=1,nmax
         np(i)=0
 13   continue

C     Figure out with systems are multi's
      do 12 i=1,npt
         koinum=int(koi(i)+0.1)
         koimax=max(koinum,koimax) !largest KOI number
         np(koinum)=np(koinum)+1
 12   continue

      j=0
      do 43 i=1,npt
         koinum=int(koi(i)+0.1)
         if((np(koinum).gt.1).and.(flag(i).lt.1))then
            j=j+1
            xp(j)=real(per(i))
         endif
 43   continue

      call rqsortr(j,xp,ns)

      do 44 i=1,j
         yp(i)=xp(ns(i))
 44   continue
      ncd=0
      do 45 i=1,j
         xp(i)=real(i)/real(j)
         if((yp(i).ge.pc).and.(ncd.eq.0)) ncd=i
c         write(6,*) 1,yp(i),xp(i)
 45   continue
      write(6,500) "Multi CD(",pc,"): ",xp(ncd),"+/-",
     .   sqrt(xp(ncd)*(1-xp(ncd))*real(j))/real(j),j
 500  format(A9,F5.2,2(A3,F6.3),1X,I5)

C     Plot cummulative disp for 'good' multis.
      call pgline(j,yp,xp)

c     okay.. lets plot 'good' singles
      j=0
      do 46 i=1,npt
         koinum=int(koi(i)+0.1)
         if((np(koinum).eq.1).and.(flag(i).lt.1))then
            j=j+1
            xp(j)=real(per(i))
         endif
 46   continue

      call rqsortr(j,xp,ns)

      do 47 i=1,j
         yp(i)=xp(ns(i))
 47   continue
      ncd=0
      do 48 i=1,j
         xp(i)=real(i)/real(j)
         if((yp(i).ge.pc).and.(ncd.eq.0)) ncd=i
c         write(6,*) 2,yp(i),xp(i)
 48   continue
      write(6,500) "Singl CD(",pc,"): ",xp(ncd),"+/-",
     .   sqrt(xp(ncd)*(1-xp(ncd))*real(j))/real(j),j
      ncd=0

C     Plot cummulative disp for 'good' singles
      call pgsci(2)
      call pgline(j,yp,xp)
      call pgsci(1)

c     okay.. lets plot the poor KOIs
      j=0
      do 49 i=1,npt
         if((flag(i).eq.1).and.(sn(i).gt.7.1))then
            j=j+1
            xp(j)=real(per(i))
         endif
 49   continue

      call rqsortr(j,xp,ns)

      do 50 i=1,j
         yp(i)=xp(ns(i))
 50   continue
      ncd=0
      do 51 i=1,j
         xp(i)=real(i)/real(j)
         if((yp(i).ge.pc).and.(ncd.eq.0)) ncd=i
c         write(6,*) 3,yp(i),xp(i)
 51   continue
      write(6,500) "FP    CD(",pc,"): ",xp(ncd),"+/-",
     .   sqrt(xp(ncd)*(1.0-xp(ncd))*real(j))/real(j)
      ncd=0

C     Plot cummulative disp for poor KOIs
      call pgsci(3)
      call pgline(j,yp,xp)
      call pgsci(1)

C     now deal with the EBs
      call rqsort(nEB,PEB,ns)

      do 52 i=1,nEB
         xp(i)=real(i)/real(nEB)
         yp(i)=PEB(ns(i))
 52   continue

c     okay.. lets plot the P/T0 collisions
      j=0
      do 57 i=1,npt
         if((flag(i).eq.2).and.(sn(i).gt.7.1))then
            j=j+1
            xp(j)=real(per(i))
         endif
 57   continue

      call rqsortr(j,xp,ns)

      do 58 i=1,j
         yp(i)=xp(ns(i))
 58   continue
      ncd=0
      do 59 i=1,j
         xp(i)=real(i)/real(j)
         if((yp(i).ge.pc).and.(ncd.eq.0)) ncd=i
c         write(6,*) 4,yp(i),xp(i)
 59   continue
      write(6,500) "P/T0  CD(",pc,"): ",xp(ncd),"+/-",
     .   sqrt(xp(ncd)*(1-xp(ncd))*real(j))/real(j)
      ncd=0

C     Plot cummulative disp for poor KOIs
      call pgsci(5)
      call pgline(j,yp,xp)
      call pgsci(1)

C     now deal with the EBs
      call rqsort(nEB,PEB,ns)

      ncd=0
      do 60 i=1,nEB
         xp(i)=real(i)/real(nEB)
         yp(i)=PEB(ns(i))
         if((yp(i).ge.pc).and.(ncd.eq.0)) ncd=i
c         write(6,*) 5,yp(i),xp(i)
 60   continue
      write(6,500) "EB    CD(",pc,"): ",xp(ncd),"+/-",
     .   sqrt(xp(ncd)*(1-xp(ncd))*real(nEB))/real(nEB)
      ncd=0


C     Plot cummulative disp for EBs
      call pgsci(4)
      call pgline(nEB,yp,xp)
      call pgsci(1)

      call pgclos()

      goto 999
901   write(0,*) "Cannot open: ",filename
      goto 999
902   write(0,*) "Error reading on line: ",i+2
      goto 999
903   write(0,*) "Cannot open: ",filenameEB
      goto 999
999   end

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
