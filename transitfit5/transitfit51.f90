program transitfit51
!F90 version.
use precision
implicit none
integer :: nfrho,iargc,nunit,nmax,npt,nptv,nfit,nplanet,npta,i,         &
 nplanetmax
integer, allocatable, dimension(:) :: dtype,ntt
real(double) :: ztime,chisq,rhoi,rhoierr(9)
real(double), allocatable, dimension(:) :: time,flux,ferr,itime,vtime,  &
 vel,verr,vetime,sol,err,aT,aM,aE,aIT,tmodel
real(double), allocatable, dimension(:,:) :: serr,tobs,omc
character(80) :: obsfile,rvfile,inputsol,ttfile
character(3) :: titles(18)

interface
   subroutine readdata2(obsfile,nmax,npt,time,flux,ferr,itime,ztime)
      use precision
      implicit none
      integer :: nmax,npt
      real(double) :: ztime
      real(double), dimension(:) :: time,flux,ferr,itime
      character(80) :: obsfile
   end subroutine readdata2
   subroutine fittransitmodel3(nfit,sol,serr,nplanet,npt,aT,aM,aE,aIT,     &
 dtype,ntt,tobs,omc,nfrho,rhoi,rhoierr)
      use precision
      implicit none
      integer, target :: nfit,npt,nplanet,nfrho
      integer, dimension(:), target :: dtype,ntt
      real(double), target :: rhoi
      real(double), dimension(:), target :: sol,aT,aM,aE,aIT,rhoierr
      real(double), dimension(:,:), target :: serr,tobs,omc
   end subroutine fittransitmodel3
end interface

nfrho=1 !flag for folding in rho-star priors

if(iargc().lt.3)then
   write(6,*) "Usage: transitfit5 <photfile> <rvfile> <fitpars> [ttfiles]"
   stop
endif

!maximum number of data points
nmax=2000000

call getarg(1,obsfile)
allocate(time(nmax),flux(nmax),ferr(nmax),itime(nmax))
ztime=54900.0
call readdata2(obsfile,nmax,npt,time,flux,ferr,itime,ztime)
write(0,*) "Number of data points read: ",npt

call getarg(2,rvfile) !get filename for RV data
allocate(vtime(nmax),vel(nmax),verr(nmax),vetime(nmax))
if(rvfile.eq.'null')then
   nptv=0
else
   call readdata2(rvfile,nmax,nptv,vtime,vel,verr,vetime,ztime)
endif

call getarg(3,inputsol) !get filename for input solution
nfit=108
allocate(sol(nfit),serr(nfit,2),err(nfit))
nunit=10 !unit number used for file input
open(unit=nunit,file=inputsol,status='old',err=902)
write(0,*) "reading in input solution"
call getfitpars(nunit,nfit,nplanet,sol,serr,err)
write(0,*) "done reading input solution"
close(nunit) !release unit number as we are done with file
goto 903 !goofy goto to use F77 code
902 write(0,*) "Cannot open ",inputsol
stop
903 continue

!total number of data points read
npta=npt+nptv
allocate(dtype(npta),aT(npta),aM(npta),aE(npta),aIT(npta))
if(npt.gt.0)then
   dtype(1:npt)=0 !mark as photometry
   aT(1:npt)=time(1:npt)
   aM(1:npt)=flux(1:npt)
   aE(1:npt)=ferr(1:npt)
   aIT(1:npt)=itime(1:npt)
endif
!deallocate initial arrays for photometry
deallocate(time,flux,ferr,itime)
if(nptv.gt.0)then
   dtype(npt+1:npt+nptv)=1 !mark as RVs
   aT(npt+1:npt+nptv)=vtime(1:nptv)
   aM(npt+1:npt+nptv)=vel(1:nptv)
   aE(npt+1:npt+nptv)=verr(1:nptv)
   aIT(npt+1:npt+nptv)=vetime(1:nptv)
endif
deallocate(vtime,vel,verr,vetime)
!deallocate initial RV arrays

allocate(ntt(nplanet),tobs(nplanet,npta),omc(nplanet,npta))
do i=1,nplanet
   if(iargc().ge.3+i)then
      call getarg(3+i,ttfile)
      if(ttfile.eq.'null')then
         ntt(i)=0
      else
         nunit=10
         open(unit=nunit,file=ttfile,status='old',err=905)
         goto 906
          905 write(0,*) "Cannot open ", ttfile
          stop
         906 continue
         call readttfile(nunit,nplanet,npta,i,ntt,tobs,omc)
         close(nunit)
      endif
   else
      ntt(i)=0
   endif
enddo

call fittransitmodel3(nfit,sol,serr,nplanet,npt,aT,aM,aE,aIT,     &
 dtype,ntt,tobs,omc,nfrho,rhoi,rhoierr)

call exportfit(nfit,nplanet,sol,serr,err,titles)

nmax=npta
nplanetmax=nplanet
allocate(tmodel(npta))
call transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npta,aT,aIT,         &
 ntt,tobs,omc,tmodel,dtype)

chisq=0.0d0
!$OMP PARALLEL DO REDUCTION(+:chisq)
do i=1,npta
   chisq=chisq+((aM(i)-tmodel(i))/aE(i))**2.0d0
enddo
!$OMP END PARALLEL DO
write(0,*) "Chisq: ",chisq

end program transitfit51

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine exportfit(nfit,nplanet,sol,serr,err,titles)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: nfit,i,nplanet
real(double) :: sol(nfit),serr(nfit,2),err(nfit)
character(3) :: tit
character(3) :: titles(18)
titles(1)='RHO'
titles(2)='NL1'
titles(3)='NL2'
titles(4)='NL3'
titles(5)='NL4'
titles(6)='DIL'
titles(7)="VOF"
titles(8)="ZPT"
titles(9)="EP1"
titles(10)="PE1"
titles(11)='BB1'
titles(12)='RD1'
titles(13)='EC1'
titles(14)='ES1'
titles(15)='KR1'
titles(16)='TE1'
titles(17)='EL1'
titles(18)='AL1'

open(unit=11,file="newfit.dat")

sol(1)=abs(sol(1))

do i=1,nplanet
   sol(11+10*(i-1))=abs(sol(11+10*(i-1)))
enddo

do i=1,nplanet*10+8
   if(i.le.8)then
      tit=titles(i)
   else
      tit=titles(i-10*((i-9)/10))
!      write(0,*) titles(i-10*((i-9)/10))
      write(tit(3:3),501) ((i-9)/10)+1
 501  format(I1)
!      write(0,*) i,i-10*((i-9)/10),tit
   endif
   write(11,500) tit,sol(i),serr(i,1),serr(i,2),err(i)
enddo
 500  format(A3,5(1X,1PE17.10))

      close(11)

      return
      end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine readdata2(obsfile,nmax,npt,time,flux,ferr,itime,ztime)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: nmax,npt
real(double) :: ztime
real(double), dimension(:) :: time,flux,ferr,itime
character(80) :: obsfile
!local vars
integer :: nunit,filestatus,i
real(double) :: t,f,e,sec2day,it

sec2day=86400.0d0

nunit=10
open(unit=nunit,file=obsfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",obsfile
   stop
endif

i=0
do
   if(i.gt.nmax)then
      write(0,*) "Increase nmax to match data points"
      write(0,*) "nmax: ",nmax
      stop
   endif
   read(nunit,*,iostat=filestatus) t,f,e,it
   it=0.0 !Forces the use of Long Cadence
   if(filestatus == 0) then
      i=i+1
      time(i)=t-ztime+0.5d0
      flux(i)=f+1.0
      ferr(i)=e
      if (it.lt.0.5) then
         itime(i)=1765.5/sec2day
      else
         itime(i)=58.85/sec2day !short cadence
      endif
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i+1
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit) !close file
npt=i

return
end subroutine readdata2
