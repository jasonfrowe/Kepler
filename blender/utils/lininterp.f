      subroutine lininterp(tage,varin,ntrack,agein,varout)
      implicit none
      integer ntrack,i
      double precision tage(ntrack),varin(ntrack),agein,varout
      
C     Default is to set varout to non-sense to force new MC
      varout=-1.0d0

c      write(0,*) (rhoin(i),i=1,npt)
c      write(0,*) (rhoierr(i),i=1,npt)
      do 10 i=1,ntrack-1
        if((agein.gt.tage(i)).and.(agein.le.tage(i+1)))then
            varout=varin(i)+(agein-tage(i))/(tage(i+1)-tage(i))*
     .          (varin(i+1)-varin(i))
        endif
 10   continue
c      write(0,*) drho,dsig
c      read(5,*)
 
      return
      end