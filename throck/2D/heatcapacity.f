CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine heatcapacity(nlat,nlon,Cp,T,dEn,M)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nlat,nlon,i,j
      double precision Cp,T(nlat,nlon),dEn(nlat,nlon),M(nlat)
      
      do 10 j=1,nlon
        do 11 i=1,nlat
            T(i,j)=T(i,j)+dEn(i,j)/Cp/M(i)
 11     continue
 10   continue
 
      return
      end