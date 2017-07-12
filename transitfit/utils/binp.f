CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine binp(npt,phase,mag,merr,bins,flag)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,bins,N,i,j,flag(npt)
      parameter(N=1200000)
      double precision phase(npt),tt(N),mag(npt),tn(N),merr(npt),te(n),
     .  tp(N)
      
c      write(0,*) "Bins:",bins

      do 5 i=1,bins
         tn(i)=0.0
         tp(i)=0.0
         tt(i)=0.0
         te(i)=0.0
         flag(i)=0
 5    continue

      do 10 i=1,npt
         j=int(real(bins)*phase(i))+1
         if(j.le.bins)then
            tn(j)=tn(j)+1.0/merr(i)
            tp(j)=tp(j)+phase(i)/merr(i)
            tt(j)=tt(j)+mag(i)/merr(i)
            te(j)=te(j)+1.0
         endif
 10   continue

      do 20 i=1,bins
c         phase(i)=(dble(i)-0.5)/dble(bins)
         phase(i)=tp(i)/tn(i)
         mag(i)=tt(i)/tn(i)
         merr(i)=(te(i)**0.5)/tn(i)
         if(tn(i).eq.0.0) flag(i)=1
c         write(0,*) flag(i)
c         read(5,*)
c         write(6,*) phase(i),mag(i),merr(i)
 20   continue
      npt=bins

      return
      end