CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine boxcar(npt,time,mag,merr,boxbin)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Robust Boxcar filter for unevenly sampled data.
C
C     This routine assumes that the data are sorted in time from 
C     earliest to latest. 
      implicit none
      integer npt,nmax,i,j,k,ii,dcsum,dcsumold,it,itmax
      parameter(nmax=650000)
      integer dc(nmax),nd(nmax)
      double precision time(nmax),mag(nmax),merr(nmax),boxbin,diffm,
     .  diffp,bbd2,avgm(nmax),avge(nmax),abin(nmax),std,adif(nmax),mean,
     .  stdev,sigcut,pts(nmax)
      
      sigcut=3.0
      itmax=10
      it=0
      
      do 3 i=1,npt
        dc(i)=0
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
        if(adif(i).gt.mean+sigcut*std) then
            dc(i)=1
c            write(11,*) time(i),mag(i)
        else
            dc(i)=0
        endif
c        write(0,*) "adif:",adif(i),sigcut*std
c        read(5,*)
        dcsum=dcsum+dc(i)
 15   continue
      write(0,*) dcsum,dcsumold
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

