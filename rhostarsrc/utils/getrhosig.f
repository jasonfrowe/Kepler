      subroutine getrhosig(rhoierr,rhoin,npt,drho,dsig)
      implicit none
      integer npt,i
      double precision rhoierr(npt),rhoin(npt),drho,dsig
      
C     Default is to reject drho is not within bounds below
      dsig=100.0d0

c      write(0,*) (rhoin(i),i=1,npt)
c      write(0,*) (rhoierr(i),i=1,npt)
      do 10 i=1,npt-1
        if((drho.gt.rhoierr(i)).and.(drho.le.rhoierr(i+1)))then
            dsig=rhoin(i)+(drho-rhoierr(i))/(rhoierr(i+1)-rhoierr(i))*
     .          (rhoin(i+1)-rhoin(i))
        endif
 10   continue
c      write(0,*) drho,dsig
c      read(5,*)
 
      return
      end