CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function distance(asep,eccn,Tanom)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision asep,eccn,Tanom
      
      distance=asep*(1.0d0-eccn*eccn)/(1+eccn*cos(Tanom))
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function trueanomaly(eccn,Eanom)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision eccn,Eanom,temp(2)
      
      temp(1)=sqrt((1.0d0+eccn)/(1.0d0-eccn))
      temp(2)=tan(Eanom/2.)
      trueanomaly=2.0d0*atan(temp(1)*temp(2))
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine kepler(Manom,Eanom,eccn)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer i,itmax
      parameter(itmax=100)
      double precision Manom,Eanom,Eold,eccn,diff,thres
      thres=1.0d-8

      Eold=Eanom
      Eanom=Manom+eccn*sin(Eanom)
      diff=abs(1.0d0-Eanom/Eold)
      Eold=Eanom
      i=0
      do while ((diff.gt.thres).and.(i.lt.itmax))
        Eanom=Manom+eccn*sin(Eanom)
        diff=abs(1.0d0-Eanom/Eold)
        Eold=Eanom
        i=i+1
      enddo

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine invkepler(Eanom,Manom,eccn)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer i,itmax
      parameter(itmax=100)
      double precision Manom,Eanom,Mold,eccn,diff,thres
      thres=1.0d-6

      Mold=Manom
      Manom=Eanom-eccn*sin(Manom)
      diff=abs(1.0d0-Manom/Mold)
      Mold=Manom
      i=0
      do while ((diff.gt.thres).and.(i.lt.itmax))
        Manom=Eanom-eccn*sin(Manom)
        diff=abs(1.0d0-Manom/Mold)
        Mold=Manom
        i=i+1
      enddo
c      if(i.ge.itmax) write(0,*) "invkepler itmax"

      return
      end