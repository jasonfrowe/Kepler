program transitfit7
use precision
use ocmod
implicit none
integer :: iargc,nmax,npt,nbodies,i,np
real(double) :: tol
real(double), allocatable, dimension(:) :: time,flux,ferr,itime,sol,ans, &
 Pers
real(double), allocatable, dimension(:,:) :: serr
character(80) :: inputsol,photfile

!below are all the interfaces to allow dynamic arrays.
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
   subroutine readinputsol(inputsol,nbodies,sol,serr)
      use precision
      implicit none
      character(80), intent(in) :: inputsol
      integer, intent(inout) :: nbodies
      real(double), dimension(:), intent(inout) :: sol
      real(double), dimension(:,:), intent(inout) :: serr
   end subroutine readinputsol
end interface

interface
   subroutine lcmodel(nbodies,npt,tol,sol,time,itime,ans)
      use precision
      implicit none
      integer, intent(inout) :: nbodies,npt
      real(double), intent(inout) :: tol
      real(double), dimension(:), intent(inout) :: sol,time,itime
      real(double), dimension(:), intent(inout) :: ans
   end subroutine lcmodel
end interface

interface
   subroutine fittransitmodel(nbodies,npt,tol,sol,serr,time,flux,ferr,itime)
      use precision
      implicit none
      integer, intent(inout) :: nbodies,npt
      real(double), intent(inout) :: tol
      real(double), dimension(:), intent(inout) :: sol,time,flux,ferr,itime
      real(double), dimension(:,:), intent(inout) :: serr
   end subroutine fittransitmodel
end interface

interface
   subroutine exportfit(nbodies,sol,serr)
      use precision
      implicit none
      integer nbodies
      real(double), dimension(:), intent(inout) :: sol
      real(double), dimension(:,:), intent(inout) :: serr
   end subroutine exportfit
end interface
!end of interfaces

tol=1.0d-12  !default tolerance level

if(iargc().lt.2) then  !check that we have sufficient commandline arguments
   write(0,*) "Usage: transitfit7 <photometry> <inputsol>"
   stop
endif

call getarg(1,photfile)  !get filename for photometry
nmax=100000  !hard-coded maximum number of points to read -- will update
allocate(time(nmax),flux(nmax),ferr(nmax),itime(nmax)) !allocate photometry
call readkeplc(photfile,npt,time,flux,ferr,itime) !read in the Kepler phot
write(0,*) "Number of observations read: ",npt !report number of obs

call getarg(2,inputsol) !get name of input solution
call calcnbodies(inputsol,nbodies) !first pass to get number of bodies
write(0,*) "nbodies: ",nbodies
allocate(sol(7+nbodies*7),serr(7+nbodies*7,2)) !use nbodies to allocate
allocate(m(nbodies))
call readinputsol(inputsol,nbodies,sol,serr) !read in initial solution

!The following is for enabling determination of model TTVs
allocate(Pers(nbodies-1))
do i=2,nbodies
   np=7+7*(i-1)
   Pers(i-1)=sol(np+2)
enddo
ntmidmax=(maxval(time)-minval(time))/minval(Pers)*2 !need to add checks for overflows
!write(0,*) "ntmidmax:",ntmidmax
allocate(ntmid(nbodies),tmid(nbodies,ntmidmax))
ntmid=0
deallocate(Pers)

!call fittransitmodel(nbodies,npt,tol,sol,serr,time,flux,ferr,itime)

!write(0,*) "Exporting fit"
!call exportfit(nbodies,sol,serr)

allocate(ans(npt)) !ans contains the model to match the data.
call lcmodel(nbodies,npt,tol,sol,time,itime,ans) !generate a LC model.
do i=1,npt
   write(6,501) time(i),flux(i),ans(i)
enddo
501 format(5(1X,1PE17.10))

end program transitfit7
