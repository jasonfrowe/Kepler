CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine boxcar(npt,time,mag,merr,boxbin)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Boxcar filter for unevenly sampled data.
C
C     This routine assumes that the data are sorted in time from 
C     earliest to latest. 
      implicit none
      integer npt,nmax,i,j,k,ii
      parameter(nmax=1650000)
      double precision time(nmax),mag(nmax),merr(nmax),boxbin,diffm,
     .  diffp,bbd2,avgm(nmax),avge(nmax),abin(nmax)
      
C     Initialize arrays to zero.
      do 5 i=1,npt
        avgm(i)=0.0d0
        avge(i)=0.0d0
        abin(i)=0.0d0
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
            avgm(k)=avgm(k)+mag(ii)/merr(ii)
            avge(k)=avge(k)+1.0
            abin(k)=abin(k)+1.0/merr(ii)
 11     continue
        avgm(k)=avgm(k)/abin(k)
        avge(k)=(avge(k)**0.5)/abin(k)
 10   continue
 
      i=1
      diffm=time(i)-time(1)
      do while(diffm.lt.bbd2)
        i=i+1
        diffm=time(i)-time(1)
      enddo
      j=npt
      diffp=time(npt)-time(j)
      do while(diffp.lt.bbd2)
        j=j-1
        diffp=time(npt)-time(j)
      enddo
 
      npt=0
      do 12 k=i,j
        npt=npt+1
        time(npt)=time(k)
        mag(npt)=mag(k)-avgm(k)
        merr(npt)=sqrt(merr(k)*merr(k)+avge(k)*avge(k))
 12   continue
 
      
      return
      end
        