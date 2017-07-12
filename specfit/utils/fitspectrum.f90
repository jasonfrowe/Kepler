subroutine fitspectrum(npt,wv,flux,nmodel,wmod,fmod,fmodbin,nfit)
use precision
use fittingmod
implicit none
integer :: npt,n,lwa,info,i,j,nfit,k,maxiter
integer, target :: nmodel
integer, allocatable, dimension(:) :: iwa,q
real(double) :: tollm,norm,rv,med,chi,chiold,dchi
real(double), dimension(:) :: fmodbin
real(double), dimension(:), target :: wv,flux,wmod,fmod
real(double), allocatable, dimension(:), target :: serr
real(double), allocatable, dimension(:) :: sol,fvec,wa,ratio
external fcn

interface
   subroutine binmodel(npt,wv,nmodel,wmod,fmod,fmodbin,rv)
      use precision
      implicit none
      integer, intent(inout) :: npt,nmodel
      real(double), intent(inout) :: rv
      real(double), dimension(:), intent(inout) :: wv,wmod,fmod
      real(double), dimension(:), intent(inout) :: fmodbin
   end subroutine
end interface

wv2 => wv
flux2 => flux
nmodel2 => nmodel
wmod2 => wmod
fmod2 => fmod

!order of fit
n=nfit+2 !1=R/V (m/s), 2=zero-point
allocate(sol(n))
sol=0.0d0

rv=sol(1)
call binmodel(npt,wv,nmodel,wmod,fmod,fmodbin,rv)
allocate(ratio(npt),q(npt))
ratio=flux(1:npt)/fmodbin !ratio of observed and model
call rqsort(npt,ratio,q) !sort to get median
sol(2)=ratio(q(npt/2)) !estimate of spectrum normalization
deallocate(ratio,q)

lwa=npt*n+5*npt*n
allocate(fvec(npt),iwa(n),wa(lwa))
tollm=1.0d-12

allocate(serr(npt)) !allocate array for data weights.
serr=1.0d0
serr2 => serr

call lmdif1(fcn,npt,n,sol,fvec,tollm,info,iwa,wa,lwa) !chi-sq min
!write(0,*) sol
rv=sol(1)
call binmodel(npt,wv,nmodel,wmod,fmod,fmodbin,rv) !bin-data
do i=1,npt
   norm=0.0d0
   do j=3,n
      norm=norm+sol(j)*(wv(i)-wv(1))**(j-3)
   enddo
   fmodbin(i)=fmodbin(i)*norm+sol(2) !scale data according to weights.
enddo

maxiter=1000 !maximum iterations allowed
chiold=100.0
dchi=0.0
k=0 !count interactions
allocate(q(npt))
do while((dchi.gt.1.0e-8).and.(k.lt.maxiter))
   do i=1,npt
      serr(i)=abs(flux(i)-fmodbin(i))**2.0d0
   enddo
   call rqsort(npt,serr,q)
   med=serr(q(npt/2)) !get median difference
   serr=serr+med/10.0 !no errors=0
   call lmdif1(fcn,npt,n,sol,fvec,tollm,info,iwa,wa,lwa) !chi-sq min
!   write(0,*) sol
   rv=sol(1)
   call binmodel(npt,wv,nmodel,wmod,fmod,fmodbin,rv) !bin-data
   do i=1,npt
      norm=0.0d0
      do j=3,n
         norm=norm+sol(j)*(wv(i)-wv(1))**(j-3)
      enddo
      fmodbin(i)=fmodbin(i)*norm+sol(2) !scale data according to weights.
   enddo
!   write(0,*) "Chi: ",(sum(flux(1:npt)-fmodbin))**2.0
   chi=(sum(flux(1:npt)-fmodbin))**2.0
   dchi=abs(chi-chiold)/chiold
!   write(0,*) chi,chiold,dchi
   chiold=chi
   k=k+1 !counter to break infinite loops
enddo
deallocate(q)
if(k.ge.maxiter) write(0,*) "warning, max iterations reached"
write(0,*) "sol:",sol


return
end subroutine fitspectrum

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine fcn(npt,n,x,fvec,iflag)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
use fittingmod
implicit none
integer :: npt,n,iflag,i,j
real(double) :: norm,rv
real(double), dimension(n) :: x
real(double), dimension(npt) :: fvec
real(double), allocatable, dimension(:) :: fmodbin

interface
   subroutine binmodel(npt,wv,nmodel,wmod,fmod,fmodbin,rv)
      use precision
      implicit none
      integer, intent(inout) :: npt,nmodel
      real(double), intent(inout) :: rv
      real(double), dimension(:), intent(inout) :: wv,wmod,fmod
      real(double), dimension(:), intent(inout) :: fmodbin
   end subroutine
end interface

!write(0,*) x

allocate(fmodbin(npt))
fmodbin=0.0d0
rv=x(1)
call binmodel(npt,wv2,nmodel2,wmod2,fmod2,fmodbin,rv)

do i=1,npt
   norm=0.0d0
   do j=3,n
      norm=norm+x(j)*(wv2(i)-wv2(1))**(j-3)
   enddo
   fvec(i)=(flux2(i)-fmodbin(i)*norm-x(2))/serr2(i)
enddo

return
end subroutine
