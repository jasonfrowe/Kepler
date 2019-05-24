subroutine calcimpact(nbodies,y,sol,b_cur)
use precision
implicit none
!import vars
integer :: nbodies
real(double), dimension(:) :: y,sol,b_cur
!local vars
integer :: i
real(double) :: fourpid3,rhostar,Mstar,rstar,xp,yp,zp,xs,ys,zs,LU2

!check if we have any transits and initialize b_old
fourpid3=4.0d0*Pi/3.0d0
LU2=LU*LU
rhostar=abs(sol(1)*1000.0) !mean stellar density (kg/m^3)
Mstar=abs(sol(12)) !mass of central object (MEarth)
rstar=(Mearth*Mstar/(fourpid3*rhostar))**(1.0d0/3.0d0)

xs=y(1)
ys=y(2)
zs=y(3)
do i=2,nbodies
	xp=y(6*i-5) !X
	yp=y(6*i-4) !Y
	zp=y(6*i-3) !Z
	if(yp-ys.gt.0.0d0)then !condition for a transit
		b_cur(i)=(zp-zs)*(zp-zs)+(xp-xs)*(xp-xs)
		b_cur(i)=sqrt(LU2*b_cur(i))/rstar
	else
		b_cur(i)=1.0e30 !no transit, wrong side of star
	endif
enddo

return
end