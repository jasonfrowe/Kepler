CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine radiate(nlat,nlon,nR,sbolt,dtime,sarea,rad,T,dEn)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nlat,nlon,nR,i,j,k
      double precision sbolt,dtime,sarea(nlat),T(nlat,nlon,nR),
     .  dEn(nlat,nlon,nR),rad(nR)
      
      k=1
      do 10 j=1,nlon
        do 11 i=1,nlat
            dEn(i,j,k)=dEn(i,j,k)-dtime*sarea(i)*rad(k)*rad(k)*sbolt*
     .          T(i,j,k)**4.0d0
 11     continue
 10   continue
 
      return
      end
            