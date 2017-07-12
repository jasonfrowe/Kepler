CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine heatcapacity(nlat,nlon,nR,Cp,T,dEn,M)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nlat,nlon,nR,i,j,k
      double precision Cp(nlat,nlon,nR),T(nlat,nlon,nR),
     .  dEn(nlat,nlon,nR),M(nlat,nR)
      
      do 12 k=1,nR
        do 10 j=1,nlon
            do 11 i=1,nlat
                T(i,j,k)=T(i,j,k)+dEn(i,j,k)/Cp(i,j,k)/M(i,k)
 11         continue
 10     continue
 12   continue
 
      return
      end