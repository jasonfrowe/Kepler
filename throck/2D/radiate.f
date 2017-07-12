CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine radiate(nlat,nlon,sbolt,dtime,sarea,rad,T,dEn)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nlat,nlon,i,j
      double precision sbolt,dtime,sarea(nlat),T(nlat,nlon),
     .  dEn(nlat,nlon),rad
      
      
      do 10 j=1,nlon
        do 11 i=1,nlat
           dEn(i,j)=dEn(i,j)-dtime*sarea(i)*rad*rad*sbolt*T(i,j)**4.0d0
 11     continue
 10   continue
 
      return
      end
            