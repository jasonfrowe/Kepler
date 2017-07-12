CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine histogram(npt,rp,dp,dp2,np,nbin,nbinmax,bdatax,bdatay,
     .  title,rmed,std,errs)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i,j,npt2,nbin,nbinmax,nmode
      integer np(npt)
      real rp(npt),datamin,datamax,bdatax(nbinmax),bdatay(nbinmax),bmax,
     .  rave,errs(6),std,rmed,mode,fc
      double precision dp(npt),dp2(npt),ave,var,med
      character*80 title

C     Find median
      call rqsort(npt,dp,np) !changed npt to k
      i=npt/2
      if(i.le.0) i=1
      ave=dp(np(i))
C     First we remove the average (helps a lot with double -> real)
      call avevar2(dp,dp2,npt,ave,var)
      rave=real(ave) !convert to real*4
      std=real(sqrt(var))
c      write(0,*) ave,std
      
C     Now we convert dp to rp
      j=0 !because we have sigma-clipping, we need a counter.
      do 10 i=1,npt
        if(abs(dp(i)-ave).lt.4.0*sqrt(var))then
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
 
      call bindata(nbin,npt2,rp,bdatax,bdatay,datamin,datamax,bmax)

!     calculating modes
      fc=(datamax-datamin)/real(nbin)/2.0d0 !center histograms
      mode=bdatax(1)+fc
      nmode=bdatay(1)
      do 13 i=2,nbin
         if(bdatay(i).gt.nmode)then
            nmode=bdatay(i)
            mode=bdatax(i)+fc
         endif
 13   continue
      rmed=mode

      call pgpage() !fresh plotting surface
      call pgslw(1)
      call pgsch(2.0)
c         call windowsetup(xb1,xb2,yb1,yb2) !make a square plotting surface
c         call pgvport(xb1,xb2,yb1,yb2)
      call pgvport(0.2,0.99,0.2,0.9)
c      fc=2.*(datamax-datamin)/real(nbin) !center histograms
      call pgwindow(datamin,datamax,0.,bmax+0.1*bmax) !set size
         
C        Add axis labels
      call pglabel(title,"Relative Probability","")
      call pgbin(nbin,bdatax,bdatay,.false.) !plot the histogram
      
C     Shift axis scale to account for average removal
      call pgwindow(datamin+rave,datamax+rave,0.,1.0+0.1*1.0)
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0) !add boarders
      
      call errorest(npt2,rp,nbin,bdatax,bdatay,bmax,rmed,errs)
c      call errorest2(nbin,bdatax,bdatay,bmax,rmed,errs)
      rmed=rmed+rave !correct for average removal
C     Need to recalulate standard deviation at this point!      

      return
      end      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine errorest2(nbin,bdatax,bdatay,bmax,frsol,errs)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nbin,i,j,k,nbinmax,mdpnt
      parameter(nbinmax=500)
      integer bsy(nbinmax)
      real bdatax(nbin),bdatay(nbin),bmax,frsol,errs(6),sumby,shgt,
     .  inth,intl,sumval,per,perold,inthold,intlold,e1,e2,bx
      logical loop
      
      sumby=0
      do 6 i=1,nbin
        sumby=sumby+bdatay(i)
 6    continue
c      write(0,*) "sumby:",sumby
      
      call rqsortr(nbin,bdatay,bsy)
      bx=bdatax(bsy(nbin))
      
C     find out which bin contains the middle.      
      do 5 i=1,nbin-1
c        write(0,*) frsol,bdatax(i),bdatax(i+1)
        if((bdatax(i).le.frsol).and.(bdatax(i+1).gt.frsol))then
            mdpnt=i
        endif
 5    continue
C     case when "bestfit" is outside histogram (i.e. eccentricity)
      if(frsol.le.bdatax(1)) mdpnt=1
      if(frsol.ge.bdatax(nbin)) mdpnt=nbin 
c      write(0,*) "mdpnt:",mdpnt,bdatax(mdpnt)

      perold=0 !initalization of variables
      intlold=0
      inthold=0
      do 14 i=1,6  !initialize variable to zero.
        errs(i)=0.0
 14   continue
      do 13 k=nbin,1,-1
        shgt=real(bdatay(bsy(k)))
c        write(0,*) shgt

        loop=.true. !loop until loop is false
        i=mdpnt  !removed -1
        do 10 while(loop)
            if((bdatay(i).ge.shgt).and.(bdatay(i+1).lt.shgt))then
                loop=.false.
                inth=real(i)!bdatax(i)
c                inth=(bdatay(i+1)-shgt)/(bdatay(i+1)-bdatay(i))
c     .              *(bdatax(i+1)-bdatax(i))+
c     .              bdatax(i)
            endif
            if(i.ge.nbin-1) then
                loop=.false.
                inth=real(nbin)!bdatax(nbin)
            endif
            i=i+1
c           write(6,*) "i:",i
 10     enddo
 
C       now we move in the reverse direction
        loop=.true. !loop until loop is false
        i=mdpnt+1
        do 11 while(loop)
            if((bdatay(i).ge.shgt).and.(bdatay(i-1).lt.shgt))then
                loop=.false.
                intl=real(i)!bdatax(i)
c                intl=bdatax(i)-
c     .              (bdatay(i)-shgt)/(bdatay(i)-bdatay(i-1))
c     .              *(bdatax(i)-bdatax(i-1))
            endif
            if(i.le.2) then
                loop=.false.
                intl=real(1)!bdatax(1)
            endif
            i=i-1
 11     enddo

        sumval=0.0
        do 15 i=intl,inth
            sumval=sumval+bdatay(i)
 15     continue
        per=sumval/sumby
        
        if((per.ge.0.683).and.(perold.lt.0.683))then
            e1=bdatax(inth)-(per-0.683)/(per-perold)*
     .          (bdatax(inth)-bdatax(inthold))
            e2=bdatax(intl)-(per-0.683)/(per-perold)*
     .          (bdatax(intl)-bdatax(intlold))
            errs(1)=e1-bx!-frsol
            errs(2)=e2-bx!-frsol
        endif
        if((per.ge.0.954).and.(perold.lt.0.954))then
            e1=bdatax(inth)-(per-0.954)/(per-perold)*
     .          (bdatax(inth)-bdatax(inthold))
            e2=bdatax(intl)-(per-0.954)/(per-perold)*
     .          (bdatax(intl)-bdatax(intlold))
            errs(3)=e1-bx!-frsol
            errs(4)=e2-bx!-frsol
        endif
        if((per.ge.0.9973).and.(perold.lt.0.9973))then
            e1=bdatax(inth)-(per-0.9973)/(per-perold)*
     .          (bdatax(inth)-bdatax(inthold))
            e2=bdatax(intl)-(per-0.9973)/(per-perold)*
     .          (bdatax(intl)-bdatax(intlold))
            errs(5)=e1-bx!-frsol
            errs(6)=e2-bx!-frsol
        endif
        
        perold=per
        intlold=intl
        inthold=inth

c      write(0,*) shgt,intl,inth,per
c      write(0,*) e1,e2
        
 13   continue
      
c      write(0,*) frsol
c      read(5,*)      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine errorest(npt,dd,nbin,bdatax,bdatay,bmax,frsol,errs)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nbin,i,j,mdpnt,k,nbinmax
      parameter(nbinmax=500)
      integer bsy(nbinmax)
      real dd(npt),bdatax(nbin),bdatay(nbin),bmax,frsol,errs(6),shgt,
     .  intl,inth,per,perold,e1,e2,inthold,intlold
      logical loop
      
      call rqsortr(nbin,bdatay,bsy)
      
C     find out which bin contains the middle.      
      do 5 i=1,nbin-1
c        write(0,*) frsol,bdatax(i),bdatax(i+1)
        if((bdatax(i).le.frsol).and.(bdatax(i+1).gt.frsol))then
            mdpnt=i
        endif
 5    continue
C     case when "bestfit" is outside histogram (i.e. eccentricity)
      if(frsol.le.bdatax(1)) mdpnt=1
      if(frsol.ge.bdatax(nbin)) mdpnt=nbin 
c      write(0,*) "mdpnt:",mdpnt
      
      
      perold=0 !initalization of variables
      intlold=0
      inthold=0
      do 14 i=1,6  !initialize variable to zero.
        errs(i)=0.0
 14   continue
C     bmax is the value of the largest histogram bin
c      do 13 k=1,bmax-1
c        shgt=real(bmax-k)
      do 13 k=nbin,1,-1
        shgt=real(bdatay(bsy(k)))
        
C       first we move in the forward direction      
        loop=.true. !loop until loop is false
        i=mdpnt  !removed -1
        do 10 while(loop)
            if((bdatay(i).ge.shgt).and.(bdatay(i+1).lt.shgt))then
                loop=.false.
c               inth=bdatax(i)
                inth=(bdatay(i+1)-shgt)/(bdatay(i+1)-bdatay(i))
     .              *(bdatax(i+1)-bdatax(i))+
     .              bdatax(i)
            endif
            if(i.ge.nbin-1) then
                loop=.false.
                inth=bdatax(nbin)
            endif
            i=i+1
c           write(6,*) "i:",i
 10     enddo
       
C       now we move in the reverse direction
        loop=.true. !loop until loop is false
        i=mdpnt+1
        do 11 while(loop)
            if((bdatay(i).ge.shgt).and.(bdatay(i-1).lt.shgt))then
                loop=.false.
c               intl=bdatax(i)
                intl=bdatax(i)-
     .              (bdatay(i)-shgt)/(bdatay(i)-bdatay(i-1))
     .              *(bdatax(i)-bdatax(i-1))
            endif
            if(i.le.2) then
                loop=.false.
                intl=bdatax(1)
            endif
            i=i-1
 11     enddo

        j=0
        do 12 i=1,npt
            if((dd(i).gt.intl).and.(dd(i).lt.inth))then
                j=j+1    
            endif
 12     continue
 
        per=real(j)/real(npt)
        if((per.ge.0.683).and.(perold.lt.0.683))then
            e1=inth-(per-0.683)/(per-perold)*(inth-inthold)
            e2=intl+(per-0.683)/(per-perold)*(intlold-intl)
            errs(1)=e1-frsol
            errs(2)=e2-frsol
c           write(6,*) "1 sigma:",frsol,errs(1),errs(2)
        endif
        if((per.ge.0.954).and.(perold.lt.0.954))then
            e1=inth-(per-0.954)/(per-perold)*(inth-inthold)
            e2=intl+(per-0.954)/(per-perold)*(intlold-intl)
            errs(3)=e1-frsol
            errs(4)=e2-frsol
c           write(6,*) "2 sigma:",frsol,errs(3),errs(4)
        endif 
        if((per.ge.0.9973).and.(perold.lt.0.9973))then
            e1=inth-(per-0.9973)/(per-perold)*(inth-inthold)
            e2=intl+(per-0.9973)/(per-perold)*(intlold-intl)
            errs(5)=e1-frsol
            errs(6)=e2-frsol
c           write(6,*) "3 sigma:",frsol,errs(5),errs(6)
        endif       
        perold=per
        intlold=intl
        inthold=inth
c       write(6,*) per,intl,inth,shgt
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
