program transitfit7
use precision
implicit none
integer :: iargc,nmax,npt,nbodies,i,np,itprint,itmodel
real(double) :: tol
real(double), allocatable, dimension(:) :: time,flux,ferr,itime,sol,ans, &
 Pers,percor
real(double), allocatable, dimension(:,:) :: serr
character(80) :: inputsol,photfile
!percor vars
integer :: npt2
real(double) :: tsamp,ts,te
real(double), allocatable, dimension(:) :: time2,itime2,ans2
!ocmod vars
integer :: ntmidmax
integer, allocatable, dimension(:) :: ntmid
real(double), allocatable, dimension(:,:) :: tmid

!below are all the interfaces to allow dynamic arrays.
interface
   subroutine readkeplc(photfile,npt,time,flux,ferr,itime)
      use precision
      implicit none
      character(80), intent(in) :: photfile
      integer, intent(out) :: npt
      real(double), dimension(:), intent(inout) :: time,flux,ferr,itime
   end subroutine readkeplc
   subroutine calcnbodies(inputsol,nbodies)
      use precision
      implicit none
      character(80), intent(in) :: inputsol
      integer, intent(out) :: nbodies
   end subroutine calcnbodies
   subroutine readinputsol(inputsol,nbodies,sol,serr)
      use precision
      implicit none
      character(80), intent(in) :: inputsol
      integer, intent(inout) :: nbodies
      real(double), dimension(:), intent(inout) :: sol
      real(double), dimension(:,:), intent(inout) :: serr
   end subroutine readinputsol
   subroutine lcmodel(nbodies,npt,tol,sol,time,itime,ntmid,tmid,percor,ans,itprint,itmodel)
      use precision
      implicit none
      integer, intent(inout) :: nbodies
      integer, intent(inout) :: npt,itprint,itmodel
      integer, dimension(:), intent(inout) :: ntmid
      real(double), intent(inout) :: tol
      real(double), dimension(:), intent(inout) :: sol,time,itime,percor
      real(double), dimension(:), intent(inout) :: ans
      real(double), dimension(:,:), intent(inout) :: tmid
   end subroutine lcmodel
   subroutine fittransitmodel(nbodies,npt,tol,sol,serr,time,flux,ferr,itime,ntmidmax)
      use precision
      implicit none
      integer, intent(inout) :: nbodies,npt,ntmidmax
      real(double), intent(inout) :: tol
      real(double), dimension(:), intent(inout) :: sol,time,flux,ferr,itime
      real(double), dimension(:,:), intent(inout) :: serr
   end subroutine fittransitmodel
   subroutine exportfit(nbodies,sol,serr)
      use precision
      implicit none
      integer nbodies
      real(double), dimension(:), intent(inout) :: sol
      real(double), dimension(:,:), intent(inout) :: serr
   end subroutine exportfit
   subroutine percorcalc(nbodies,sol,ntmidmax,ntmid,tmid,percor)
      use precision
      implicit none
      !import vars
      integer, intent(in) :: nbodies,ntmidmax
      integer, dimension(:), intent(in) :: ntmid
      real(double), dimension(:), intent(in) :: sol
      real(double), dimension(:,:), intent(in) :: tmid
      real(double), dimension(:), intent(inout):: percor
   end subroutine percorcalc
end interface
!end of interfaces

tol=1.0d-8  !default tolerance level (can probably be d-12 -- needs testing)

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
!allocate(m(nbodies))
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

call fittransitmodel(nbodies,npt,tol,sol,serr,time,flux,ferr,itime,ntmidmax)
write(0,*) "Exporting fit"
call exportfit(nbodies,sol,serr)

!for percorcalc
allocate(percor(nbodies))

allocate(ans(npt)) !ans contains the model to match the data.

!setting up percorcalc
tsamp=maxintg/86400.0 !sampling [days]  !1-5 min seems to be fine for Kepler.
ts=minval(time(1:npt))
te=maxval(time(1:npt))
npt2=int((te-ts)/tsamp)+1
allocate(time2(npt2),itime2(npt2),ans2(npt2))
do i=1,npt2
   time2(i)=ts+dble(i)*tsamp
   itime2(i)=tsamp
enddo
itprint=0 !no output of timing measurements
itmodel=0 !do not need a transit model
percor=0.0d0 !initialize percor to zero.
call lcmodel(nbodies,npt2,tol,sol,time2,itime2,ntmid,tmid,percor,ans2,itprint,itmodel) !generate a LC model.
call percorcalc(nbodies,sol,ntmidmax,ntmid,tmid,percor)
itprint=1 !create output of timing measurements
itmodel=0 !do not calculate a transit model
call lcmodel(nbodies,npt2,tol,sol,time2,itime2,ntmid,tmid,percor,ans2,itprint,itmodel)
itprint=0 !do not create output of timing measurements
itmodel=1 !calculate a transit model
call lcmodel(nbodies,npt,tol,sol,time,itime,ntmid,tmid,percor,ans,itprint,itmodel)
!call percorcalc(nbodies,sol,percor)
do i=1,npt
   write(6,501) time(i),flux(i),ans(i)
enddo

!itprint=0
!call lcmodel(nbodies,npt,tol,sol,time,itime,percor,ans,itprint) !generate a LC model.
!call percorcalc(nbodies,sol,percor)
!itprint=1
!call lcmodel(nbodies,npt,tol,sol,time,itime,percor,ans,itprint) !generate a LC model.
!do i=1,npt
!   write(6,501) time(i),ans(i)
!enddo


!do i=1,npt
!   write(6,501) time(i),flux(i),ans(i)
!enddo
501 format(5(1X,1PE17.10))

end program transitfit7
