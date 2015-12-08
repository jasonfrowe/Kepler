      program q1read
      implicit none
      integer nT,nstar,nunit,nbzi,i,j
      parameter(nT=1639,nstar=156097)
      integer kID(nstar),nf(nT)
      real*8 time(nT)
      real flux(nT,nstar),ff(nT),med,median
      character*80 filename
      
      nunit=10
      filename="/home/rowe/Kepler/q1/flux.dat"
      
C     Open the binary data file
      nbzi=nT*nstar
      open(unit=nunit,file=filename,status='old',form='unformatted',
     .  access='direct',recl=nbzi,err=901)
      read(nunit,rec=1) flux
      close(nunit) 
      
      filename="/home/rowe/Kepler/q1/kepId.txt"
      open(unit=nunit,file=filename,status='old',err=901)
      do 10 i=1,nstar
        read(nunit,*) kID(i)
 10   continue
      close(nunit)
      
      nbzi=2*nT
      filename="/home/rowe/Kepler/q1/cadenceTimes.dat"
      open(unit=nunit,file=filename,status='old',form='unformatted',
     .  access='direct',recl=nbzi,err=901)
      read(nunit,rec=1) time
      close(nunit)
      
      do 11 i=1,nstar
        if(kID(i).lt.10)then
            write(filename,507) "klc0000000",kID(i),".dat"
 507        format(A10,I1,A4)      
        elseif(kID(i).lt.100)then
            write(filename,506) "klc000000",kID(i),".dat"
 506        format(A9,I2,A4)      
        elseif(kID(i).lt.1000)then
            write(filename,505) "klc00000",kID(i),".dat"
 505        format(A8,I3,A4)      
        elseif(kID(i).lt.10000)then
            write(filename,504) "klc0000",kID(i),".dat"
 504        format(A7,I4,A4)      
        elseif(kID(i).lt.100000)then
            write(filename,503) "klc000",kID(i),".dat"
 503        format(A6,I5,A4)      
        elseif(kID(i).lt.1000000)then
            write(filename,502) "klc00",kID(i),".dat"
 502        format(A5,I6,A4)
        elseif(kID(i).lt.10000000)then
            write(filename,501) "klc0",kID(i),".dat"
 501        format(A4,I7,A4)
        else
            write(filename,500) "klc",kID(i),".dat"
 500        format(A3,I8,A4)
        endif
        
c        write(6,508) filename
c 508    format(A15)
        
        med=median(nT,nstar,flux,nf,ff,i)
        open(unit=nunit,file=filename)
        do 12 j=1,nT
            if(flux(j,i).gt.0.0)then
                write(nunit,*) time(j),(flux(j,i)-med)/med
            endif
 12     continue       
        close(nunit)
 
 11   continue
      
      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 999  end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function median(nT,nstar,flux,nf,ff,i)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nT,nstar,i,j,k
      integer nf(nT)
      real flux(nT,nstar),ff(nT)
      
      k=0
      do 10 j=1,nT
        if(flux(j,i).gt.0.0)then
            k=k+1
            ff(k)=flux(j,i)
        endif
 10   continue
 
      call rqsort(k,ff,nf)
      
      median=ff(nf(k/2))      
      
      return
      end
      
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