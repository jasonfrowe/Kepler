program transitfit7
use precision
implicit none
integer :: iargc,nmax,npt,nbodies,i,np,itprint,itmodel,colflag
real(double) :: tol
real(double), allocatable, dimension(:) :: time,flux,ferr,itime,sol,ans, &
 Pers,percor
real(double), allocatable, dimension(:,:) :: serr
character(80) :: inputsol,photfile
!percor vars
!integer :: npt2
!real(double) :: tsamp,ts,te
!real(double), allocatable, dimension(:) :: time2,itime2,ans2
!ocmod vars
integer :: ntmidmax
integer, allocatable, dimension(:) :: ntmid
real(double), allocatable, dimension(:,:) :: tmid

!below are all the interfaces to allow dynamic arrays.
interface
   subroutine readkeplc(photfile,npt,time,flux,ferr,itime)
      use precision
      implicit none
      character(80) :: photfile
      integer :: npt
      real(double), dimension(:) :: time,flux,ferr,itime
   end subroutine readkeplc
   subroutine calcnbodies(inputsol,nbodies)
      use precision
      implicit none
      character(80) :: inputsol
      integer :: nbodies
   end subroutine calcnbodies
   subroutine readinputsol(inputsol,sol,serr)
      use precision
      implicit none
      character(80) :: inputsol
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: serr
   end subroutine readinputsol
   subroutine lcmodel(nbodies,npt,tol,sol,time,itime,ntmidmax,ntmid,tmid,percor,ans,colflag,itprint,itmodel)
      use precision
      implicit none
      integer :: nbodies,npt,itprint,itmodel,colflag,ntmidmax
      integer, dimension(:) :: ntmid
      real(double) :: tol
      real(double), dimension(:) :: sol,time,itime,percor
      real(double), dimension(:) :: ans
      real(double), dimension(:,:) :: tmid
   end subroutine lcmodel
   subroutine lcmodel_pc(nbodies,npt,tol,sol,time,ntmidmax,ntmid,tmid,percor,colflag,itprint)
      use precision
      implicit none
      integer :: nbodies,colflag,itprint,npt,ntmidmax
      real(double) :: tol
      real(double), dimension(:) :: sol,time,percor
      integer, dimension(:) :: ntmid !used with octiming 
      real(double), dimension(:,:) :: tmid !used with octiming
   end subroutine lcmodel_pc
   subroutine fittransitmodel(nbodies,npt,tol,sol,serr,time,flux,ferr,itime,ntmidmax)
      use precision
      implicit none
      integer :: nbodies,npt,ntmidmax
      real(double) :: tol
      real(double), dimension(:) :: sol,time,flux,ferr,itime
      real(double), dimension(:,:) :: serr
   end subroutine fittransitmodel
   subroutine exportfit(nbodies,sol,serr)
      use precision
      implicit none
      integer nbodies
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: serr
   end subroutine exportfit
   subroutine percorcalc(nbodies,sol,ntmidmax,ntmid,tmid,percor)
      use precision
      implicit none
      !import vars
      integer :: nbodies,ntmidmax
      integer, dimension(:) :: ntmid
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: tmid
      real(double), dimension(:) :: percor
   end subroutine percorcalc
end interface
!end of interfaces

tol=1.0d-8  !default tolerance level (can probably be d-12 -- needs testing)

if(iargc().lt.2) then  !check that we have sufficient commandline arguments
   write(0,*) "Usage: transitfit7 <photometry> <inputsol>"
   stop
endif

call getarg(1,photfile)  !get filename for photometry
nmax=200000  !hard-coded maximum number of points to read -- will update
allocate(time(nmax),flux(nmax),ferr(nmax),itime(nmax)) !allocate photometry
call readkeplc(photfile,npt,time,flux,ferr,itime) !read in the Kepler phot
write(0,*) "Number of observations read: ",npt !report number of obs
if(npt.gt.nmax)then
   write(0,*) "Increase nmax to at least: ",npt
   stop
endif

call getarg(2,inputsol) !get name of input solution
call calcnbodies(inputsol,nbodies) !first pass to get number of bodies
write(0,*) "nbodies: ",nbodies
allocate(sol(7+nbodies*7),serr(7+nbodies*7,2)) !use nbodies to allocate
call readinputsol(inputsol,sol,serr) !read in initial solution

!The following is for enabling determination of model TTVs
allocate(Pers(nbodies-1))
do i=2,nbodies
   np=7+7*(i-1)
   Pers(i-1)=sol(np+2)
enddo
ntmidmax=int((maxval(time(1:npt))-minval(time(1:npt)))/minval(Pers(1:nbodies-1)))*2 !need to add checks for overflows
!write(0,*) "ntmidmax:",ntmidmax,nbodies
allocate(ntmid(nbodies),tmid(nbodies,ntmidmax))
ntmid=0
deallocate(Pers)

write(0,*) "Fitting data"
call fittransitmodel(nbodies,npt,tol,sol,serr,time,flux,ferr,itime,ntmidmax)
write(0,*) "Exporting fit"
call exportfit(nbodies,sol,serr)

!for percorcalc
allocate(percor(nbodies))

allocate(ans(npt)) !ans contains the model to match the data.

itprint=1 !no output of timing measurements, create percor
percor=0.0d0 !initialize percor to zero.
call lcmodel_pc(nbodies,npt,tol,sol,time,ntmidmax,ntmid,tmid,percor,colflag,itprint) !generate a LC model.
if (colflag.eq.0) call percorcalc(nbodies,sol,ntmidmax,ntmid,tmid,percor)
itprint=1 !create output of timing measurements, no percor
call lcmodel_pc(nbodies,npt,tol,sol,time,ntmidmax,ntmid,tmid,percor,colflag,itprint) !generate a LC model.
itprint=0 !do not create output of timing measurements
itmodel=1 !calculate a transit model
call lcmodel(nbodies,npt,tol,sol,time,itime,ntmidmax,ntmid,tmid,percor,ans,colflag,itprint,itmodel)
do i=1,npt
   write(6,501) time(i),flux(i),ans(i)
enddo

501 format(5(1X,1PE17.10))

end program transitfit7
