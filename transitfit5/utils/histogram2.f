CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine histogram(npt,rp,dp,dp2,np,nbin,nbinmax,bdatax,bdatay,
     .  title,rmed,std,errs,bfit)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i,j,npt2,nbin,nbinmax
      integer np(npt)
      real rp(npt),datamin,datamax,bdatax(nbinmax),bdatay(nbinmax),bmax,
     .  rave,errs(6),std,rmed
      double precision dp(npt),dp2(npt),ave,var,med,bfit
      character*80 title

      ave=bfit
      call avevar2(dp,dp2,npt,ave,var)
      rave=real(ave) !convert to real*4
      std=real(sqrt(var))
      
      j=0
      do 10 i=1,npt
            j=j+1
            dp2(j)=dp(i)-ave
            rp(j)=real(dp2(j))
 10   continue
      npt2=j

      call rqsort(npt2,dp2,np)
      med=bfit-ave
      rmed=real(med)

      datamin=rp(1)
      datamax=rp(1)
      do 12 i=2,npt2
        datamin=min(rp(i),datamin)
        datamax=max(rp(i),datamax)
 12   continue
 
      call bindata(nbin,npt2,rp,bdatax,bdatay,datamin,datamax,bmax)
      
      call pgpage() !fresh plotting surface
      call pgslw(1)
      call pgsch(2.0)
      call pgvport(0.2,1.0,0.2,0.9)
      call pgwindow(datamin,datamax,0.,bmax+0.1*bmax) !set size
         
C        Add axis labels
      call pglabel(title,"Relative Probability","")
      call pgbin(nbin,bdatax,bdatay,.false.) !plot the histogram
      
C     Shift axis scale to account for average removal
      call pgwindow(datamin+rave,datamax+rave,0.,1.0+0.1*1.0)
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0) !add boarders
      
      call errorest(npt2,dp2,np,nbin,bdatax,bdatay,bmax,rmed,errs)
      rmed=rmed+rave !correct for average removal

      return
      end      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer function count_range(npt, dp2, np, intl, inth)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt, np(npt), l, r, mid, ilow, ihigh
      double precision dp2(npt), intl, inth
      
      if (npt .le. 0 .or. inth .lt. intl) then
         count_range = 0
         return
      endif

      if (dp2(np(npt)) .lt. intl) then
         count_range = 0
         return
      endif
      if (dp2(np(1)) .ge. intl) then
         ilow = 1
      else
         l = 1
         r = npt
 10      if (r - l .gt. 1) then
            mid = (l + r) / 2
            if (dp2(np(mid)) .ge. intl) then
               r = mid
            else
               l = mid
            endif
            goto 10
         endif
         ilow = r
      endif

      if (dp2(np(1)) .gt. inth) then
         count_range = 0
         return
      endif
      if (dp2(np(npt)) .le. inth) then
         ihigh = npt
      else
         l = 1
         r = npt
 20      if (r - l .gt. 1) then
            mid = (l + r) / 2
            if (dp2(np(mid)) .le. inth) then
               l = mid
            else
               r = mid
            endif
            goto 20
         endif
         ihigh = l
      endif

      if (ihigh .ge. ilow) then
         count_range = ihigh - ilow + 1
      else
         count_range = 0
      endif
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine errorest(npt, dp2, np, nbin, bdatax, bdatay, bmax,
     .  frsol, errs)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt, nbin, np(npt), i, k, mdpnt, count_range, cnt
      integer bsy(500)
      real bdatax(nbin), bdatay(nbin), bmax, frsol, errs(6)
      double precision dp2(npt), dx, xc(500), shgt, intl, inth, per,
     .  perold, intlold, inthold, e1, e2, dy
      logical found

      do 4 i = 1, 6
         errs(i) = 0.0
 4    continue
      if (nbin .gt. 500 .or. nbin .le. 1) return

      dx = dble(bdatax(2) - bdatax(1))
      do 5 i = 1, nbin
         xc(i) = dble(bdatax(i)) + 0.5d0 * dx
 5    continue

      mdpnt = 1
      do 6 i = 1, nbin - 1
         if (dble(bdatax(i)) .le. dble(frsol) .and.
     .       dble(bdatax(i+1)) .gt. dble(frsol)) then
            mdpnt = i
         endif
 6    continue
      if (dble(frsol) .ge. dble(bdatax(nbin))) mdpnt = nbin

      call rqsortr(nbin, bdatay, bsy)

      perold = 0.0d0
      intlold = dble(frsol)
      inthold = dble(frsol)

      do 13 k = nbin, 1, -1
         shgt = dble(bdatay(bsy(k)))
         if (shgt .le. 0.0d0) goto 13

         found = .false.
         do 10 i = mdpnt, nbin - 1
            if (dble(bdatay(i)) .ge. shgt .and.
     .          dble(bdatay(i+1)) .lt. shgt) then
               dy = dble(bdatay(i+1) - bdatay(i))
               if (abs(dy) .gt. 1.0d-12) then
                  inth = xc(i) + (shgt - dble(bdatay(i))) / dy * dx
               else
                  inth = xc(i)
               endif
               found = .true.
               goto 11
            endif
 10      continue
 11      if (.not. found) then
            inth = dble(bdatax(nbin)) + dx
         endif

         found = .false.
         do 20 i = mdpnt, 2, -1
            if (dble(bdatay(i)) .ge. shgt .and.
     .          dble(bdatay(i-1)) .lt. shgt) then
               dy = dble(bdatay(i-1) - bdatay(i))
               if (abs(dy) .gt. 1.0d-12) then
                  intl = xc(i) + (shgt - dble(bdatay(i))) / dy * (-dx)
               else
                  intl = xc(i)
               endif
               found = .true.
               goto 21
            endif
 20      continue
 21      if (.not. found) then
            intl = dble(bdatax(1))
         endif

         cnt = count_range(npt, dp2, np, intl, inth)
         per = dble(cnt) / dble(npt)

         if (per .ge. 0.682689d0 .and. perold .lt. 0.682689d0) then
            if (per .gt. perold) then
               e1 = inthold + (0.682689d0 - perold) / (per - perold) *
     .              (inth - inthold)
               e2 = intlold + (0.682689d0 - perold) / (per - perold) *
     .              (intl - intlold)
            else
               e1 = inth
               e2 = intl
            endif
            errs(1) = real(e1 - dble(frsol))
            errs(2) = real(e2 - dble(frsol))
         endif

         if (per .ge. 0.954500d0 .and. perold .lt. 0.954500d0) then
            if (per .gt. perold) then
               e1 = inthold + (0.954500d0 - perold) / (per - perold) *
     .              (inth - inthold)
               e2 = intlold + (0.954500d0 - perold) / (per - perold) *
     .              (intl - intlold)
            else
               e1 = inth
               e2 = intl
            endif
            errs(3) = real(e1 - dble(frsol))
            errs(4) = real(e2 - dble(frsol))
         endif

         if (per .ge. 0.997300d0 .and. perold .lt. 0.997300d0) then
            if (per .gt. perold) then
               e1 = inthold + (0.997300d0 - perold) / (per - perold) *
     .              (inth - inthold)
               e2 = intlold + (0.997300d0 - perold) / (per - perold) *
     .              (intl - intlold)
            else
               e1 = inth
               e2 = intl
            endif
            errs(5) = real(e1 - dble(frsol))
            errs(6) = real(e2 - dble(frsol))
         endif

         perold = per
         inthold = inth
         intlold = intl
 13   continue

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

      binsize=(datamax-datamin)/real(nbin-1)
      
      do 20 i=1,nbin
         bdatay(i)=0.
 20   continue

      do 10 i=1,npt
         bin=int((pdata(i)-datamin)/binsize)+1
         if((bin.gt.0).and.(bin.le.nbin)) bdatay(bin)=bdatay(bin)+1.0
 10   continue

      tsum=0.
      bmax=0.
      do 30 i=1,nbin
         bmax=max(bdatay(i),bmax)
         bdatax(i)=datamin+real(i-1)*binsize
         tsum=tsum+bdatay(i)
 30   continue
 
      return
      end

c**********************************************************************
      subroutine rqsortr(n,a,p)
c======================================================================
      implicit none
      integer   n
      real      a(n)
      integer   p(n)
      integer   LGN, Q
      parameter (LGN=32, Q=11)
      real      x
      integer   stackl(LGN),stackr(LGN),s,t,l,m,r,i,j

      stackl(1)=1
      stackr(1)=n
      s=1

      do 1 i=1,n
         p(i)=i
    1 continue

    2 if (s.gt.0) then
         l=stackl(s)
         r=stackr(s)
         s=s-1

    3    if ((r-l).lt.Q) then
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
