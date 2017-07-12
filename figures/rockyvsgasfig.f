      program rockyvsgas
C     Figure for Evidence of a Bi-model Planet Population from Kepler..
      implicit none
      integer nmax,i,nunit,kic,npt,j,k,nd,ii,jj
      parameter(nmax=4000)
      integer nRp(nmax)
      real xp(nmax),yp(nmax),dumr,koi,kepmag(nmax),T0,T0err,Per,Pererr,
     .   Rp(nmax),asemi,Teq,tdur,tdepth,adrs,adrserr,rdr,rdrerr,b,berr,
     .   rnpt,rmin,rmax,deriv(nmax),maxd,xp2(nmax),bx(4),by(4),snr,chi,
     .   Teff(nmax),logg,rstar,Teff2(nmax)
      character*80 filename,dumc

      rmin=1.2
      rmax=7.0

      nunit=10
      filename="appendixTable_2012Feb26.txt"
      open(unit=nunit,file=filename,status='old',err=901)

      do 10 i=1,57
         read(nunit,*) dumc
 10   continue

      i=1
 11   read(nunit,*,end=12) koi,kic,kepmag(i),T0,T0err,Per,Pererr,Rp(i),
     .   asemi,Teq,tdur,tdepth,adrs,adrserr,rdr,rdrerr,b,berr,snr,chi,
     .   Teff(i),logg,rstar
         i=i+1
      goto 11
 12   continue
      close(nunit)

      npt=i-1
      write(0,*) "Number of planets: ",npt

C     All the 'pg..' calls are for plotting purposes only.
      call pgopen('?')
      call PGPAP ( 6.0 ,1.0)
      call pgsubp(1,2)
      call pgpage()
      call pgvport(0.15,0.85,0.0,0.7)
      call pgwindow(log10(rmin),log10(rmax),0.0,1.1)
      call pgslw(3)
      call pgsch(2.5)
      CALL PGBOX('BCLTS',0.0,0,'BCNTS1',0.0,0)
      call pglabel(" ","Cumulative Fraction","")

      call rqsortr(npt,Rp,nRp)
      rnpt=real(npt)

      j=0
      do 13 i=1,npt
         if((rp(nRp(i)).gt.rmin).and.(Kepmag(i).gt.1.0))then
            j=j+1
            xp(j)=log10(rp(nRp(i)))
            yp(j)=real(i)/rnpt
            Teff2(j)=Teff(i)
         endif
 13   continue
      call pgline(j,xp,yp)
      write(0,*) "Number of planets used: ",j

      nd=200
      maxd=0.0
      k=0
      do 14 i=2,j-1!,nd
         ii=nd
         if(i+ii.gt.j) ii=j-i
         jj=nd
         if(i-jj.lt.1) jj=i-1
         ii=min(ii,jj)
         jj=ii
c         write(6,*) i,ii
         if(xp(i+ii).ne.xp(i-jj))then
            k=k+1
            xp2(k)=xp(i)
            deriv(k)=(yp(i+ii)-yp(i-jj))/(10.0**xp(i+ii)-10.0**xp(i-jj))
            maxd=max(maxd,deriv(k))
         endif
c         write(0,*) deriv(k)
 14   continue

      call pgpage()
      call pgvport(0.15,0.85,0.3,1.0)
      call pgwindow(log10(rmin),log10(rmax),0.0,maxd+0.1*maxd)
      call pgslw(3)
      call pgsch(2.5)

      bx(1)=log10(rmin)
      by(1)=0.0
      bx(2)=log10(1.5)
      by(2)=0.0
      bx(3)=log10(1.5)
      by(3)=maxd+0.1*maxd
      bx(4)=log10(rmin)
      by(4)=maxd+0.1*maxd

      call pgsfs(1)
      call pgsci(5)
      call pgpoly(4,bx,by)
      call pgsci(1)

      CALL PGBOX('BCNLTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Planet Radius (R\d\(2284)\u)",
     .   "Probability Density","")

      call pgline(k,xp2,deriv)

      call pgclos()


      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 999  end

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


