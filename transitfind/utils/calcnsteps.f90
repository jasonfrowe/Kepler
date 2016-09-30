!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine calcnsteps(nstep,Mstar,Rstar,freq1,freq2,nyq,ofac,n)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer nstep,n
real(double) :: Mstar,Rstar,freq1,freq2,nyq,ofac,df0,f,q,df,steps

interface
   subroutine qfunc(f,Mstar,Rstar,q)
      use precision
      real(double), intent(in) :: f,Mstar,Rstar
      real(double), intent(out) :: q
   end subroutine qfunc
end interface

!write(0,*) ofac,n,nyq
steps=ofac*(freq2-freq1)*n/nyq
df0=(freq2-freq1)/dble(steps)
f=freq1
nstep=0
!write(0,*) "steps: ",steps,freq1,freq2
do while(f.lt.freq2)
!   write(0,*) "nstep: ",nstep,f,df0
   nstep=nstep+1 !increase counter to count number of steps
   call qfunc(f,Mstar,Rstar,q)
   df=q*df0 !estimate next frequency step based on optimum q
   f=f+df !next frequency to scan
enddo

return
end
