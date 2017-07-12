CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mag2flux(npt,mag,flux)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i
      double precision mag(npt),flux(npt)
      
      do 10 i=1,npt
        flux(i)=10.0**(mag(i)/-2.5)
 10   continue
 
      return
      end