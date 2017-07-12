CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine surfacearea(nlat,theta,dlam,sarea,Pi)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nlat,i
      double precision theta(nlat),sarea(nlat),dlam,th2,th1,Pi
     
      do 10 i=1,nlat
        if(i+1.gt.nlat)then  !watch for indice problems
            th2=theta(1)+Pi
        else
            th2=theta(i+1)
        endif
        if(i-1.lt.1) then
            th1=theta(nlat)-Pi
        else
            th1=theta(i-1)
        endif
C       The relative surface element correction
        sarea(i)=dlam*(
     .      sin( (th2+theta(i))/2-Pi/2.0 )-
     .      sin( (theta(i)+th1)/2-Pi/2.0 )  )
        sarea(i)=abs(sarea(i))
 10   continue
      
      return
      end