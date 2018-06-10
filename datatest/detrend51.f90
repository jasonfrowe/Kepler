program detrend51
use precision
implicit none
integer :: iargc,nmax,npt,nfit,nunit,nplanet,i,nfitp,ftype,filestatus
integer, allocatable, dimension(:) :: ntt,tflag
real(double) :: ztime,boxbin
real(double), allocatable, dimension(:) :: time,flux,ferr,itime,sol,err
real(double), allocatable, dimension(:,:) :: serr,tobs,omc
character(80) :: obsfile,inputsol,ttfile,cline

interface
   subroutine readdata2(obsfile,nmax,npt,time,flux,ferr,itime,ztime)
     	use precision
     	implicit none
     	integer :: nmax,npt
     	real(double) :: ztime
     	real(double), dimension(:) :: time,flux,ferr,itime
     	character(80) :: obsfile
   end subroutine readdata2
   subroutine marktransit(np,nplanet,npt,time,tflag,nfit,sol,ntt,tobs,omc)
		use precision
		implicit none
		integer :: np,nplanet,npt,nfit
		integer, dimension(:) :: tflag,ntt
		real(double), dimension(:) :: time,sol
		real(double), dimension(:,:) :: tobs,omc
	end subroutine marktransit
	subroutine polyfilter(npt,time,flux,ferr,tflag,boxbin,nfitp)
		use precision
		implicit none
		integer :: npt,nfitp
		integer, dimension(:) :: tflag
		real(double) :: boxbin
		real(double), dimension(:) :: time,flux,ferr
	end subroutine polyfilter
end interface

if(iargc().lt.3)then
   write(0,*) "Usage: detrend51 <ftype> <filein> <boxbin> [n1.dat]"
   write(0,*) "<ftype> 0=polyfilter 1=boxcar"
   write(0,*) "<filein> Kepler photometry"
   write(0,*) "<boxbin> width of filter in days"
   write(0,*) "[n1.dat] transit solution (optional)"
   stop
endif

!maximum number of data points
nmax=2000000

call getarg(1,cline)
read(cline,*,iostat=filestatus) ftype
if((filestatus == 0).and.((ftype.eq.1).or.(ftype.eq.0))) then
	continue
else
    write(0,'(A20,A80)') "**ERROR** <ftype> : ",cline
	write(0,*) "Usage: detrend51 <ftype> <filein> <boxbin> [n1.dat]"
    write(0,*) "<ftype> 0=polyfilter 1=boxcar"
    write(0,*) "<filein> Kepler photometry"
    write(0,*) "<boxbin> width of filter in days"
    write(0,*) "[n1.dat] transit solution (optional)"
    stop
endif

call getarg(3,cline)
read(cline,*,iostat=filestatus) boxbin
if((filestatus == 0).and.(boxbin.gt.0.0)) then
	continue
else
    write(0,'(A20,A80)') "**ERROR** <boxbin> : ",cline
	write(0,*) "Usage: detrend51 <ftype> <filein> <boxbin> [n1.dat]"
    write(0,*) "<ftype> 0=polyfilter 1=boxcar"
    write(0,*) "<filein> Kepler photometry"
    write(0,*) "<boxbin> width of filter in days"
    write(0,*) "[n1.dat] transit solution (optional)"
    stop
endif

call getarg(2,obsfile)
allocate(time(nmax),flux(nmax),ferr(nmax),itime(nmax))
ztime=54900.0
call readdata2(obsfile,nmax,npt,time,flux,ferr,itime,ztime)
write(0,*) "Number of data points read: ",npt

nplanet=0
if(iargc().ge.4) then	
	call getarg(4,inputsol) !get filename for input solution
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

	allocate(ntt(nplanet),tobs(nplanet,npt),omc(nplanet,npt))
	do i=1,nplanet
   		if(iargc().ge.4+i)then
      		call getarg(4+i,ttfile)
      		if(ttfile.eq.'null')then
         		ntt(i)=0
      		else
         		nunit=10
         		open(unit=nunit,file=ttfile,status='old',err=905)
         		goto 906
          		905 write(0,*) "Cannot open ", ttfile
          		stop
         		906 continue
         		call readttfile(nunit,nplanet,npt,i,ntt,tobs,omc)
         		close(nunit)
      		endif
   		else
      		ntt(i)=0
   		endif
	enddo
endif

allocate(tflag(npt))
tflag=0
do i=1,nplanet
	call marktransit(i,nplanet,npt,time,tflag,nfit,sol,ntt,tobs,omc)  
enddo

nfitp=4
      
if(ftype.eq.0)then
	call polyfilter(npt,time,flux,ferr,tflag,boxbin,nfitp)
! 9      do 10 i=1,npt
!c            if(tflag(i).eq.0) write(6,*) time(i)-0.5d0+54900.0d0,
!c     .          flux(i),ferr(i)
!            write(6,*) time(i)-0.5d0+54900.0d0,
!     .          flux(i),ferr(i)
! 10     continue
!      else
!        call boxfilter(npt,time,flux,ferr,ts,ngap,gaps,tflag,x,y,z,
!     .      boxbin,nx)
endif

do i=1,npt
	write(6,500) time(i)-0.5d0+54900.0d0,flux(i),ferr(i),itime(i)
enddo
500 format(3(F17.11,1X),I1)


end