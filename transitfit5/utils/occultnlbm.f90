subroutine occultnlbm(rprs,c1,c2,c3,c4,b0,flux,nb)
use precision
implicit none
integer :: nb
real(double) :: rprs, c1, c2, c3, c4, b0(nb), flux(nb)
!local vars
integer :: i
real(double) d, x_in, xout, delta, x, dx, A_i, A_f, fac, M_PI, norm, &
        x_out, intens

fac = 5.0d-6
M_PI = acos(-1.d0) 
norm = (-c1/10.0 - c2/6.0 - 3.0*c3/14 - c4/4 + 0.5)*2.0*M_PI

do i = 1, nb
        
        d = b0(i)

        x_in = max(d - rprs, 0.0)
        x_out = min(d + rprs, 1.0)

        if (x_in .ge. 1.0) then
                flux(i) = 1.0d0
        elseif (x_out - x_in .lt. 1.0d-7) then
                flux(i) = 1.0d0
        else
                delta = 0.0d0
                x = x_in
                dx = fac*acos(x)
                
                x = x + dx

                A_i = 0.0d0

                do while (x < x_out)
                        call area(d, x, rprs, A_f)
                        call intensity(x - dx/2.0d0, c1, c2, c3, c4, norm, intens)
                        delta = delta + (A_f - A_i) * intens
                        dx = fac * acos(x)
                        x = x + dx
                        A_i = A_f
                enddo

                dx = x_out - x + dx
                x = x_out
                call area(d, x, rprs, A_f)
                call intensity(x - dx/2.0, c1, c2, c3, c4, norm, intens)
                delta = delta + (A_f - A_i)*intens

                flux(i) = 1.0 - delta

        endif
enddo

return
end

subroutine intensity(x, c1, c2, c3, c4, norm, intens)
use precision
implicit none
!import args
real(double) :: x, c1, c2, c3, c4, norm, intens
!local args
real(double) :: sqrtmu

if (x > 0.99995) x = 0.99995
sqrtmu = (1.0 - x*x)**(0.25d0)
intens = (1. - c1*(1. - sqrtmu) - c2*(1. - (sqrtmu)**2.0) - c3*(1. - (sqrtmu)**3.0d0) - c4*(1. - (sqrtmu)**4.0))/norm

return
end

subroutine area(d, x, R, ans)
use precision
implicit none
!import args
real(double) :: d, x, r, ans
!local args
real(double) :: arg1, arg2, arg3, M_PI


M_PI=acos(-1.d0)

arg1 = (d*d + x*x - R*R)/(2.*d*x)
arg2 = (d*d + R*R - x*x)/(2.*d*R)
arg3 = (-d + x + R)*(d + x - R)*(d - x + R)*(d + x + R)
arg3 = max(arg3, 0.0d0)

if (x .le. R - d) then
        ans = M_PI*x*x
elseif (x .ge. R + d) then
        ans = M_PI*R*R
else
        ans = x*x*acos(arg1) + R*R*acos(arg2) - 0.5*sqrt(arg3)
endif

return
end

