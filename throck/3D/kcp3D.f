CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine kcp(nlat,nlon,nR,tcond,Cp,T,rho)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nlat,nlon,nR,i,j,k
      double precision tcond(nlat,nlon,nR),Cp(nlat,nlon,nR),
     .  T(nlat,nlon,nR),rho,kappa
     
      do 11 k=1,nR
        do 12 j=1,nlon
            do 13 i=1,nlat
C               calculate diffusivity and specific heat
C                         [mm^2 s^-1]     [J mole^-1 K^-1]
                if(T(i,j,k).lt.846.0)then
                    kappa=567.3/T(i,j,k)-0.062
                    Cp(i,j,k)=199.50+0.0857*T(i,j,k)-5.0d6/T(i,j,k)/
     .                  T(i,j,k)
                else
                    kappa=0.732-0.000135*T(i,j,k)
                    Cp(i,j,k)=229.32+0.0323*T(i,j,k)-47.9d-6/T(i,j,k)/
     .                  T(i,j,k)
                endif
C               average molar mass : 221.78 g mol^-1
                Cp(i,j,k)=Cp(i,j,k)*1000.0d0/221.78d0
                tcond(i,j,k)=kappa*1.0d-6*rho*Cp(i,j,k) ![W m-1 K-1]
 13         continue
 12     continue
 11   continue
 
      return
      end
                
                
