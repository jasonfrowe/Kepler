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
integer n,i,info,lwa
integer, allocatable, dimension(:) :: iwa
real(double), target :: tollm
real(double), dimension(:,:), target :: serr
real(double), allocatable, dimension(:) :: solfit,wa,fvec
!ecosw,esinw convertion vars
!integer :: np
!real(double) :: ecosw,esinw,sqecosw,sqesinw,ecc
external fcn2

allocate(solfit(7+nbodies*7))

!how many variables are we fitting?
n=0
do i=1,7+nbodies*7
   if(serr(i,2).ne.0.0d0)then
      n=n+1
      solfit(n)=sol(i)
   endif
enddo
!write(0,*) "n:",n

!assign pointers from module that is shared with FCN
tol2f => tol
nbodies2f => nbodies
sol2f => sol
serr2f => serr
time2f => time
flux2f => flux
ferr2f => ferr
itime2f => itime
!pointers for octiming calculations
ntmidmax2f => ntmidmax

lwa=npt*n+5*npt+n
allocate(fvec(npt),iwa(n),wa(lwa))
tollm=1.0d-8
write(0,*) "Calling lmdif1"
call lmdif1(fcn2,npt,n,solfit,fvec,tollm,info,iwa,wa,lwa)
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

deallocate(solfit,fvec,iwa,wa)
nullify(tol2f,nbodies2f,sol2f,serr2f,time2f,flux2f,ferr2f,itime2f,ntmidmax2f)

return
end subroutine fittransitmodel

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine fcn2(npt,n,x,fvec,iflag)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
use fittingmod
implicit none
!import vars
integer :: npt,n,iflag
real(double), dimension(n) :: x
real(double), dimension(npt) :: fvec
!local vars
integer :: i,j,itprint,itmodel,np,colflag
real(double) :: chisq
!real(double) :: ecosw,esinw,sqecosw,sqesinw,ecc
real(double), allocatable, dimension(:) :: sol3,percor
!percorcalc
!integer :: npt3
integer, allocatable, dimension(:) :: ntmid2f
!real(double) :: tsamp,ts,te
!real(double), allocatable, dimension(:) :: time3,itime3,ans3
real(double), allocatable, dimension(:,:) :: tmid2f

interface
   subroutine lcmodel(nbodies,npt,tol,sol,time,itime,ntmid,tmid,percor,ans,colflag,itprint,itmodel)
      use precision
      implicit none
      integer :: nbodies
      integer :: npt,itprint,itmodel,colflag
      integer, dimension(:) :: ntmid
      real(double) :: tol
      real(double), dimension(:) :: sol,time,itime,percor
      real(double), dimension(:) :: ans
      real(double), dimension(:,:) :: tmid
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
      integer :: nbodies,ntmidmax
      integer, dimension(:) :: ntmid
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: tmid
      real(double), dimension(:) :: percor
   end subroutine percorcalc
end interface

!write(0,*) "FCN2"

!allocate for percorcal
allocate(ntmid2f(nbodies2f),tmid2f(nbodies2f,ntmidmax2f))

allocate(sol3(7+nbodies2f*7))

j=0
do i=1,7+nbodies2f*7
   sol3(i)=sol2f(i)
enddo
do i=1,7+nbodies2f*7
   if(serr2f(i,2).ne.0.0d0)then
      j=j+1
      sol3(i)=x(j)
   endif
enddo

!for percorcalc
allocate(percor(nbodies2f))

itprint=0 !no output of timing measurements
percor=0.0d0 !initialize percor to zero.
!write(0,*) "Calling lcmodel1",npt3
call lcmodel_pc(nbodies2f,npt,tol2f,sol3,time2f,ntmid2f,tmid2f,percor,colflag,itprint) !generate a LC model.
if (colflag.eq.0) call percorcalc(nbodies2f,sol3,ntmidmax2f,ntmid2f,tmid2f,percor)
itprint=0 !no output of timing measurements
itmodel=1 !calculate a transit model
!write(0,*) "Calling lcmodel2"
call lcmodel(nbodies2f,npt,tol2f,sol3,time2f,itime2f,ntmid2f,tmid2f,percor,fvec,colflag,itprint,itmodel)

!write(0,*) "done with lcmodel"

chisq=0.0d0
!do i=1,npt
!   chisq=chisq+(fvec(i)-flux2f(i))*(fvec(i)-flux2f(i))/(ferr2f(i)*ferr2f(i))
!enddo

!$OMP PARALLEL DO
do i=1,npt
   fvec(i)=(fvec(i)-flux2f(i))/ferr2f(i)   
enddo
!$OMP END PARALLEL DO

deallocate(ntmid2f,tmid2f,sol3,percor)

!write(0,*) "FCN2.. done",chisq/dble(npt)

return
end


