CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine stellarrad(nlat,nlon,lambda,theta,dEn,L,dist,Ab,pi,fpi,
     .  sarea,rad,dtime,sangle)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nlat,nlon,i,j
      double precision dEn(nlat,nlon),L,dist,Ab,fpi,sarea(nlat),dtime,
     .  Fstar,a1,a2,pi,lambda(nlon),theta(nlat),sangle(2),alpha,rad
     
      Fstar=L/(fPi*dist*dist)*(1.0-Ab) ![J/s/m^2]
  
      do 10 j=1,nlon
        do 11 i=1,nlat
C     Calculate zenith angle of star in sky. 
            a1=theta(i)!real(i-1)*dth
            a2=lambda(j)!real(j-1)*dps
            if(a2.lt.0.0d0)a2=a2+2.0d0*Pi
            alpha=acos(sin(a1-Pi/2.0d0)*sin(sangle(1)-Pi/2.0d0)
     .      +cos(a1-Pi/2.0d0)*cos(sangle(1)-Pi/2.0d0)*cos(a2-sangle(2)))
            if(alpha.ge.Pi/2.)then
                dEn(i,j)=0.0
            else
                dEn(i,j)=Fstar*dtime*cos(alpha)*sarea(i)*rad*rad
            endif
 11     continue
 10   continue
      
      
      return
      end