subroutine fittransitmodel(nbodies,npt,tol,sol,serr,time,flux,ferr,itime,ntmidmax)
use precision
use fittingmod
implicit none
!import vars
integer :: npt
integer, target :: nbodies,ntmidmax
real(double), target :: tol
real(double), dimension(:), target :: sol,time,flux,ferr,itime
!local vars
integer n,i,info,lwa,j
integer, allocatable, dimension(:) :: iwa
real(double), target :: tollm
real(double), dimension(:,:), target :: serr
real(double), allocatable, dimension(:) :: solfit,wa,fvec
!ecosw,esinw convertion vars
integer :: np
real(double) :: ecosw,esinw,sqecosw,sqesinw,ecc
external fcn

!convert sqecosw, sqesinw to ecosw, esinw for LSQ
do i=1,nbodies
   np=7+7*(i-1)
   sqecosw=sol(np+6)  !e^1/2 cos(w)
   sqesinw=sol(np+7)  !e^1/2 sin(w)
   ecc=(sqecosw*sqecosw+sqesinw*sqesinw) !eccentricity
   sol(np+6)=sqrt(ecc)*sqecosw
   sol(np+7)=sqrt(ecc)*sqesinw
!   write(0,*) "sqecosw,sqesinw: ",sqecosw,sqesinw
!   write(0,*) "ecosw,esinw: ",sol(np+6),sol(np+7)
enddo

allocate(solfit(7+nbodies*7))

!how many variables are we fitting?
n=0
do i=1,7+nbodies*7
   if(serr(i,2).ne.0.0d0)then
      n=n+1
      solfit(n)=sol(i)
   endif
enddo
write(0,*) "n:",n

!assign pointers from module that is shared with FCN
tol2 => tol
nbodies2 => nbodies
sol2 => sol
serr2 => serr
time2 => time
flux2 => flux
ferr2 => ferr
itime2 => itime
!pointers for octiming calculations
ntmidmax2 => ntmidmax

lwa=npt*n+5*npt*n
allocate(fvec(npt),iwa(n),wa(lwa))
tollm=1.0d-8
write(0,*) "Calling lmdif1"
call lmdif1(fcn,npt,n,solfit,fvec,tollm,info,iwa,wa,lwa)
write(0,*) "info: ",info

n=0
!write(0,*) "out:"
do i=1,7+nbodies*7
   if(serr(i,2).ne.0.0d0)then
      n=n+1
      !write(0,501) sol(i),solfit(n),sol(i)-solfit(n)
      sol(i)=solfit(n)
   endif
enddo
501 format(5(1X,1PE17.10))

!convert ecosw, esinw to sqecosw, sqesinw for storage
do i=1,nbodies
   np=7+7*(i-1)
   ecosw=sol(np+6)  !ecos(w)
   esinw=sol(np+7)  !esin(w)
   ecc=sqrt(ecosw*ecosw+esinw*esinw) !eccentricity
   if(ecc.gt.0.0)then
      sol(np+6)=ecosw/sqrt(ecc)
      sol(np+7)=esinw/sqrt(ecc)
   else
      sol(np+6)=0.0d0
      sol(np+7)=0.0d0
   endif
enddo

return
end subroutine fittransitmodel

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine fcn(npt,n,x,fvec,iflag)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
use fittingmod
implicit none
!import vars
integer :: npt,n,iflag
real(double), dimension(n) :: x
real(double), dimension(npt) :: fvec
!local vars
integer :: i,j,nunit,itprint,itmodel,np,colflag
real(double) :: yp,ecosw,esinw,sqecosw,sqesinw,ecc
real(double), allocatable, dimension(:) :: sol3,m3,y3,percor
!percorcalc
integer :: npt3
integer, allocatable, dimension(:) :: ntmid2
real(double) :: tsamp,ts,te
real(double), allocatable, dimension(:) :: time3,itime3,ans3
real(double), allocatable, dimension(:,:) :: tmid2

interface
   subroutine lcmodel(nbodies,npt,tol,sol,time,itime,ntmid,tmid,percor,ans,colflag,itprint,itmodel)
      use precision
      implicit none
      integer, intent(inout) :: nbodies
      integer, intent(inout) :: npt,itprint,itmodel,colflag
      integer, dimension(:), intent(inout) :: ntmid
      real(double), intent(inout) :: tol
      real(double), dimension(:), intent(inout) :: sol,time,itime,percor
      real(double), dimension(:), intent(inout) :: ans
      real(double), dimension(:,:), intent(inout) :: tmid
   end subroutine lcmodel
   subroutine lcmodel_pc(nbodies,npt,tol,sol,time,ntmid,tmid,percor,colflag,itprint)
      use precision
      implicit none
      integer :: nbodies,colflag,itprint,npt
      real(double) :: tol
      real(double), dimension(:) :: sol,time,percor
      integer, dimension(:) :: ntmid !used with octiming 
      real(double), dimension(:,:) :: tmid !used with octiming
   end subroutine lcmodel_pc
   subroutine percorcalc(nbodies,sol,ntmidmax,ntmid,tmid,percor)
      use precision
      implicit none
      integer, intent(in) :: nbodies,ntmidmax
      integer, dimension(:), intent(in) :: ntmid
      real(double), dimension(:), intent(in) :: sol
      real(double), dimension(:,:), intent(in) :: tmid
      real(double), dimension(:), intent(inout):: percor
   end subroutine percorcalc
end interface

!write(0,*) "FCN1"

!allocate for percorcal
allocate(ntmid2(nbodies2),tmid2(nbodies2,ntmidmax2))

allocate(sol3(7+nbodies2*7))

j=0
do i=1,7+nbodies2*7
   sol3(i)=sol2(i)
enddo
do i=1,7+nbodies2*7
   if(serr2(i,2).ne.0.0d0)then
      j=j+1
      sol3(i)=x(j)
   endif
enddo


!convert ecosw, esinw to sqecosw, sqesinw for model
do i=1,nbodies2
   np=7+7*(i-1)
   ecosw=sol3(np+6)  !ecos(w)
   esinw=sol3(np+7)  !esin(w)
   ecc=sqrt(ecosw*ecosw+esinw*esinw) !eccentricity
   if(ecc.gt.0.0)then
      sol3(np+6)=ecosw/sqrt(ecc)
      sol3(np+7)=esinw/sqrt(ecc)
   else
      sol3(np+6)=0.0d0
      sol3(np+7)=0.0d0
   endif
!   write(0,*) "ecosw,esinw: ",ecosw,esinw
!   write(0,*) "sqecosw,sqesinw: ",sol3(np+6),sol3(np+7)
enddo

!!make sure masses are positive
!do i=1,nbodies2
!   np=7+7*(i-1)
!   sol3(np+5)=abs(sol3(np+5))
!enddo

!do i=1,7+nbodies2*7
!   write(0,501) sol3(i),sol2(i),sol3(i)-sol2(i),serr2(i,2)
!enddo

!for percorcalc
allocate(percor(nbodies2))

!!setting up percorcalc
!tsamp=maxintg/86400.0 !sampling [days]  !1-5 min seems to be fine for Kepler.
!ts=minval(time2(1:npt))
!te=maxval(time2(1:npt))
!npt3=int((te-ts)/tsamp)+1
!allocate(time3(npt3),itime3(npt3),ans3(npt3))
!do i=1,npt3
!   time3(i)=ts+dble(i)*tsamp
!   itime3(i)=tsamp
!enddo
itprint=0 !no output of timing measurements
itmodel=0 !do not need a transit model
percor=0.0d0 !initialize percor to zero.
!write(0,*) "Calling lcmodel1",npt3
call lcmodel_pc(nbodies2,npt,tol2,sol3,time2,ntmid2,tmid2,percor,colflag,itprint) !generate a LC model.
!call lcmodel(nbodies2,npt3,tol2,sol3,time3,itime3,ntmid2,tmid2,percor,ans3,colflag,itprint,itmodel) !generate a LC model.
!write(0,*) "Calling percorcal.."
if (colflag.eq.0) call percorcalc(nbodies2,sol3,ntmidmax2,ntmid2,tmid2,percor)
itprint=0 !no output of timing measurements
itmodel=1 !calculate a transit model
!write(0,*) "Calling lcmodel2"
call lcmodel(nbodies2,npt,tol2,sol3,time2,itime2,ntmid2,tmid2,percor,fvec,colflag,itprint,itmodel)

!nunit=10
!open(unit=nunit,file="junk.dat")
!do i=1,npt
!   write(nunit,501) time2(i),flux2(i),fvec(i)
!enddo
!close(nunit)
501 format(5(1X,1PE17.10))
!write(0,*) "lcmodel in FCN done"
!read(5,*)

yp=1.0d0
!$OMP PARALLEL DO
do i=1,npt
   fvec(i)=(fvec(i)-flux2(i))/ferr2(i)*yp  !yp is not needed. 
enddo
!$OMP END PARALLEL DO


!write(0,*) "FCN1.. done"

return
end


