!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine qfunc(f,Mstar,Rstar,q)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
real(double) :: q,f,Mstar,Rstar,fac,G,Rsun,Msun,fsec

G=6.674d-11 !N m^2 kg^-2  Gravitation constant
fac=1.083852140278d0 !(2pi)^(2/3)/pi
Rsun=696265.0d0*1000.0d0 !m  radius of Sun
Msun=1.9891d30 !kg  mass of Sun

fsec=f/86400.0d0 !convert from c/d to c/sec

q=fac*Rstar*Rsun/(G*Mstar*Msun)**(1.0/3.0)*fsec**(2.0/3.0)

return
end
