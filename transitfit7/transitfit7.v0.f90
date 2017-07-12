program transitfit7
use precision
implicit none
integer :: iargc,i,j,npt,nmax,k,nplanet
integer, target :: nbodies
real(double), target :: tol
real(double), allocatable, dimension(:), target :: y, sol,time,flux,ferr,itime
real(double), allocatable, dimension(:) :: ans,err,yerr
real(double), allocatable, dimension(:,:), target :: serr,yserr
character(80) :: inputsol,photfile
external f,jac

interface
   subroutine readkeplc(photfile,npt,time,flux,ferr,itime)
      use precision
      implicit none
      character(80), intent(in) :: photfile
      integer, intent(out) :: npt
      real(double), dimension(:), intent(inout) :: time,flux,ferr,itime
   end subroutine readkeplc
end interface

interface
   subroutine calcnbodies(inputsol,nbodies)
      use precision
      implicit none
      character(80), intent(in) :: inputsol
      integer, intent(out) :: nbodies
   end subroutine calcnbodies
end interface

interface
   subroutine readinputsol(inputsol,nbodies,y,yserr,yerr,sol,serr,err)
      use precision
      implicit none
      character(80), intent(in) :: inputsol
      integer, intent(in) :: nbodies
      real(double), dimension(:), intent(out) :: y,yerr,sol,err
      real(double), dimension(:,:), intent(out) :: yserr,serr
   end subroutine readinputsol
end interface

interface
   subroutine lcmodel(nbodies,npt,tol,y,sol,time,itime,ans)
      use precision
      implicit none
      integer, intent(in) :: nbodies,npt
      real(double), intent(in) :: tol
      real(double), dimension(:), intent(inout) :: y,ans
      real(double), dimension(:), intent(in) :: sol,time,itime
   end subroutine lcmodel
end interface

interface
   subroutine fittransitmodel(nbodies,npt,tol,y,yserr,sol,serr,time,flux,ferr,itime)
      use precision
      implicit none
      integer, intent(in) :: nbodies,npt
      real(double), intent(in) :: tol
      real(double), dimension(:), intent(in) :: time,flux,ferr,itime
      real(double), dimension(:), intent(inout) :: y,sol
      real(double), dimension(:,:), intent(in) :: yserr,serr
   end subroutine fittransitmodel
end interface

interface
   subroutine exportfit(nplanet,sol,serr,err,y,yerr,yserr)
      use precision, only: double
      implicit none
      integer, intent(in) :: nplanet
      real(double), dimension(:), intent(in) :: sol,err,y,yerr
      real(double), dimension(:,:), intent(in) :: serr,yserr
   end subroutine exportfit
end interface

tol=1.0d-12  !default tolerance level

if(iargc().lt.2) then
   write(0,*) "Usage: transitfit7 <photometry> <inputsol>"
   stop
endif

call getarg(1,photfile)
nmax=100000
allocate(time(nmax),flux(nmax),ferr(nmax),itime(nmax))
call readkeplc(photfile,npt,time,flux,ferr,itime)
write(0,*) "Number of observations read: ",npt

call getarg(2,inputsol)
call calcnbodies(inputsol,nbodies)
allocate(m(nbodies),mserr(nbodies,2),merr(nbodies),y(nbodies*6), &
   yserr(nbodies*6,2),yerr(nbodies*6),sol(8+nbodies), &
   serr(8+nbodies,2),err(8+nbodies))
write(0,*) "nbodies: ",nbodies
call readinputsol(inputsol,nbodies,y,yserr,yerr,sol,serr,err)
do i=1,nbodies
   write(0,500) m(i),(y(j),j=6*i-5,6*i)
enddo
500  format(28(1PE10.3,1X))

write(0,*) "Starting Best Fit"
call fittransitmodel(nbodies,npt,tol,y,yserr,sol,serr,time,flux,ferr,itime)
write(0,*) "Done fitting"

write(0,*) "Exporting fit"
nplanet=nbodies-1
call exportfit(nplanet,sol,serr,err,y,yerr,yserr)

write(0,*) "Generating output.."
allocate(ans(npt))
call lcmodel(nbodies,npt,tol,y,sol,time,itime,ans)
do i=1,npt
   write(6,501) time(i),flux(i),ans(i)
enddo
501 format(5(1X,1PE17.10))

end program transitfit7
