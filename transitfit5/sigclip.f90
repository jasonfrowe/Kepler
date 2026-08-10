program sigclip_new
use precision
implicit none
integer :: nmax,npt,iostatus,nfit,nunit,nplanet,iargc,i,j,imodel,npt2,nplanetmax
integer, allocatable, dimension(:) :: ntt,tflag,cflag,dtype
real(double) :: ztime,sigclip,std,my,stdev
real(double), allocatable, dimension(:) :: time,flux,ferr,itime,sol,err,phase,&
  deriv,y,pcut,tmodel
real(double), allocatable, dimension(:,:) :: serr,tobs,omc
character(80) :: obsfile,cline,inputsol,ttfile

interface
   subroutine readdata2(obsfile,nmax,npt,time,flux,ferr,itime,ztime)
      use precision
      implicit none
      integer :: nmax,npt
      real(double) :: ztime
      real(double), dimension(:) :: time,flux,ferr,itime
      character(80) :: obsfile
   end subroutine readdata2
end interface

!maximum number of data points
nmax=2000000

!Usage Message
if(iargc().lt.2)then
   write(6,*) "Usage: sigclip <photfile> <sigma> [n1.dat] [ttfiles]"
   write(6,*) "  <photfile> : 4 columns: time, flux, ferr, itime"
   write(6,*) "  <sigma> : sigma-clipping value.  Standard deviations"
   write(6,*) "  [n1.dat]: optional transit solution file to exclude in-transit data"
   write(6,*) "  [ttfiles]: optional .tt files for transit timing variations"
   stop
endif

!read in photometry
call getarg(1,obsfile)
allocate(time(nmax),flux(nmax),ferr(nmax),itime(nmax))
ztime=54900.0
call readdata2(obsfile,nmax,npt,time,flux,ferr,itime,ztime)
write(0,*) "Number of data points read: ",npt

!read in sigma-clip value from commandline

call getarg(2,cline)
read(cline,*,iostat=iostatus) sigclip
if(iostatus>0)then 
   write(0,*) "*** Invalid sigclip entry, must be a number ***"
   stop
endif
if(sigclip.le.0)then
   write(0,*) "*** sigclip must be greater than zero ***"
   stop
endif

!Read in optional transit model solution

imodel=0 !default is no model
if(iargc().ge.3)then
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
   imodel=1 !set flag to note that we have a model
endif

allocate(ntt(nplanet),tobs(nplanet,npt),omc(nplanet,npt))
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
         call readttfile(nunit,nplanet,npt,i,ntt,tobs,omc)
         close(nunit)
      endif
   else
      ntt(i)=0
   endif
enddo

allocate(tflag(npt))
tflag=0 !initiate tflag to zero. tflag=1, means in-transit

if (imodel.eq.1) then
   allocate(phase(npt))
   do i=1,nplanet
      call marktransit(i,npt,phase,time,tflag,nfit,sol,ntt,tobs,omc)
   enddo
endif

allocate(cflag(npt),deriv(npt))
cflag=0
call outlier(npt,time,flux,cflag,deriv,sigclip,tflag)

allocate(y(npt),pcut(npt))

if(imodel.eq.1)then
   allocate(dtype(npt),tmodel(npt))
   dtype=0
   nplanetmax=nplanet
   call transitmodel(nfit,nplanet,nplanetmax,sol,nmax,npt,time,itime, &
      ntt,tobs,omc,tmodel,dtype)
   pcut=flux(1:npt)-tmodel(1:npt) !if we have a model use residuals
else
   pcut=flux
endif

j=0
do i=1,npt
   if(tflag(i).eq.0)then
      j=j+1
      y(j)=pcut(i)
   endif
enddo
npt2=j
      
std=stdev(npt2,y,my)

do i=1,npt
   if(abs(pcut(i)-my).lt.sigclip*std)then
      if(cflag(i).eq.0)then
         write(6,500) time(i)-0.5d0+ztime,flux(i)-1.0,ferr(i),itime(i)
      endif
   endif
enddo
500  format(4(F17.11,1X))

end program sigclip_new

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine outlier(npt,time,flux,cflag,deriv,sigcut,tflag)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: npt,cflag(npt),i,tflag(npt)
real(double) :: time(npt),flux(npt),deriv(npt),std,mean,stdev,sigcut

do i=1,npt-1
   deriv(i)=flux(i+1)-flux(i)
enddo

std=stdev(npt-1,deriv,mean)

do i=1,npt-1
  if((abs(deriv(i)-mean).gt.std*sigcut).and.  &
    (abs(deriv(i+1)-mean).gt.std*sigcut).and. &
    (deriv(i)/deriv(i+1).lt.0.0d0))then
      if(tflag(i+1).eq.0)then 
         cflag(i+1)=1 !only mark-of-transit data
      endif
  endif
enddo 

return
end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine marktransit(np,npt,phase,time,tflag,nfit,sol,ntt,tobs,omc)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: npt,nplanetmax,nmax,nfit
parameter(nmax=2000000,nplanetmax=10)
integer :: tflag(npt),i,np,ntt(nplanetmax),col
real(double) :: time(npt),sol(nfit),ph1,ph2,transitdur,toff,tdur, &
  phase(npt),period,tobs(nplanetmax,nmax),omc(nplanetmax,nmax),   &
  tcor(nmax),ttcor,epo,tdurfac

tdur=2.0d0*transitdur(np,nfit,sol)/86400.0d0+0.03
tdur=max(tdur,2.55d0/24.0) !make sure duration is at least 2.5 hours
!write(0,*) 'tdur: ',tdur

col=10*(np-1)
epo=sol(9+col)
period=sol(10+col)

!tdurfac - 0.5: remove exactly 1 transit duration
!tdurfac - 1.0: remove +/- 1 transit duration centred on transit
tdurfac=1.0d0
ph1=0.75-tdurfac*tdur/period
if(ph1.lt.0.5)ph1=0.5
ph2=0.75+tdurfac*tdur/period
if(ph2.gt.1.0)ph2=1.0
!write(0,*) "ph1,ph2",ph1,ph2

toff=0.75-(epo/period-int(epo/period))
if(toff.lt.0.0)toff=toff+1.0
!write(0,*) "Toff:",toff

do i=1,npt
   call lininterp(tobs,omc,nplanetmax,nmax,np,ntt,time(i),ttcor)
   tcor(i)=time(i)-ttcor
!   write(0,*) tcor(i),ttcor
enddo

call phasept(npt,tcor,phase,period,toff)

do i=1,npt
!  write(0,*) phase(i),ph1,ph2
  if((phase(i).ge.ph1).and.(phase(i).le.ph2))then
      tflag(i)=1
!      write(0,*) time(i)
!      read(5,*)
  endif
enddo

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function transitdur(np,nfit,sol)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer np,nfit
real(double) :: sol(nfit),b,Psec,G,aConst,Pi,adrs,cincl,temp(4),bb,rdr,&
  transitdur

Pi=acos(-1.d0)   !Pi
G=6.674d-11 !N m^2 kg^-2  Gravitation constant
aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)

b=sol(11+10*(np-1))
bb=b*b
Psec=sol(10+10*(np-1))*8.64d4 !sec ; period of planet
adrs=1000.0*sol(1)*G*(Psec)**2/(3.0d0*Pi)
adrs=adrs**(1.0d0/3.0d0) !a/R*
cincl=b/adrs !cos(i)
rdr=sol(12+10*(np-1))
!write(0,*) bb,adrs,rdr

temp(1)=Psec/Pi
temp(2)=1.0d0/adrs
temp(3)=(1+rdr)**2.0-bb
temp(4)=1-cincl*cincl
!Transit duration in days
transitdur=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine phasept(npt,time,phase,period,toff)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     npt - number of data points
!     time - times of observations
!     mag - the observations
!     phase - returned phase of observation
!     period - fixed period for data
use precision
implicit none
integer npt
double precision time(npt),phase(npt),period,toff

integer i
double precision temp

do i=1,npt
   temp=time(i)
!   Get the phase
   phase(i)=temp/period-int(temp/period)
!   apply optional phase offset to make plot pretty
   phase(i)=phase(i)+toff
!   make sure phase is between 0 and 1
   if(phase(i).lt.0.0) phase(i)=phase(i)+1.0
   if(phase(i).gt.1.0) phase(i)=phase(i)-1.0
enddo
return
end
