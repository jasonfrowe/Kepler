program pdmsearch
use precision
implicit none
integer iargc,nmax,npt,nstep,nplot
real(double) :: ztime,nyq,mint,maxt,freq1,freq2,fread,ofac,Mstar,Rstar, &
 bper,bpow,epo,depth,pmean,std,qtran
real(double), allocatable, dimension(:) :: time,flux,ferr,itime,t1,t2,  &
 t3,p,pow
character(80) :: obsfile,cline

interface
   subroutine readdata2(obsfile,nmax,npt,time,flux,ferr,itime,ztime)
      use precision
      implicit none
      integer :: nmax,npt
      real(double) :: ztime
      real(double), dimension(:) :: time,flux,ferr,itime
      character(80) :: obsfile
   end subroutine readdata2
   subroutine nyquest(npt,time,nyq)
      use precision
      integer, intent(in) :: npt
      real(double), dimension(:), intent(in) :: time
      real(double), intent(out) :: nyq
   end subroutine nyquest
   subroutine calcnsteps(nstep,Mstar,Rstar,freq1,freq2,nyq,ofac,n)
      use precision
      integer, intent(out) :: nstep
      integer, intent(in) :: n
      real(double), intent(in) :: Mstar,Rstar,freq1,freq2,nyq,ofac
   end subroutine calcnsteps
   subroutine pdm(nstep,p,pow,npt,time,flux,ferr,freq1,freq2,nyq,ofac,  &
    Mstar,Rstar,bper,bpow)
      use precision
      implicit none
      integer, intent(in) :: npt,nstep
      real(double), intent(in) :: ofac,nyq,freq1,freq2,Mstar,Rstar
      real(double), intent(inout) :: bper,bpow
      real(double), dimension(:), intent(inout) :: p,pow
      real(double), dimension(:), intent(in) :: time,flux,ferr
   end subroutine pdm
   subroutine makeplot(nplot,nstep,freqs,p,obsfile,epo,bper,bpow,depth, &
    pmean,std,qtran,npt,time,flux)
      use precision
      integer, intent(in) :: nplot,nstep,npt
      real(double), intent(in) :: bper,epo,bpow,depth,pmean,std,qtran
      real(double), dimension(:), intent(in) :: freqs,p,time,flux
      character(80), intent(in) :: obsfile
   end subroutine makeplot
end interface

if(iargc().lt.1) then !check number of command line arguments
   write(0,*) "Usage: pdmsearch <filename> [nplot] [freq1] [freq2] [Msun] [Rsun]"
   write(0,*) "  <filename> - photometry (time,flux)"
   write(0,*) "  [nplot] - 0: no plot, 1:xwindow 2:ps 3:png"
   write(0,*) "  [Per1] - low Period (days) to scan, -1 for default"
   write(0,*) "  [Per2] - high Period (days) to scan, -1 for default"
   write(0,*) "  [Msun]  - mass of host star (Msun) "
   write(0,*) "  [Rsun]  - radius of host star (Rsun) "
   stop
endif

call getarg(1,obsfile) !get filename for input data

!maximum number of data points - will be reduced after reading in data
nmax=2000000

call getarg(1,obsfile)  !get filename from commandline
allocate(time(nmax),flux(nmax),ferr(nmax),itime(nmax)) !allocate space
ztime=54900.0  !zero-point offset (unique to Kepler)
call readdata2(obsfile,nmax,npt,time,flux,ferr,itime,ztime) !read in Kepler-LC
allocate(t1(npt),t2(npt),t3(npt))  !allocate temp-array
t1(1:npt)=time(1:npt)  !transfer to compact array
t2(1:npt)=flux(1:npt)
t3(1:npt)=ferr(1:npt)
nmax=npt
deallocate(time,flux,ferr,itime)  !re-allocate main arrays
allocate(time(nmax),flux(nmax),ferr(nmax))
time=t1 !fill with data
flux=t2
ferr=t3
deallocate(t1,t2,t3) !de-allocate temp-arrays

write(0,*) "Number of data points read: ",npt

!check that we have enough data to continue
if(npt.lt.3)then
   write(0,*) "Need at least 3 data points"
   stop
endif

!estimate Nyquest freq base on median value of dT
call nyquest(npt,time,nyq)

!set low-frequency search bounds
mint=minval(time(1:npt)) !get length of observations (days)
maxt=maxval(time(1:npt))
freq1=2.0/(maxt-mint) !set lowest frequency
if(iargc().ge.4)then   !check if commandline wants a different value
   call getarg(4,cline)
   read(cline,*) fread
   if(fread.gt.0.0d0) freq1=1.0/fread
endif

freq2=0.5 !Default value for high-frequency (c/d)
if(iargc().ge.3)then
   call getarg(3,cline)
   read(cline,*) fread
   if(fread.gt.0.0d0) freq2=1.0/fread
endif

write(0,*) "freqs: ",freq1,freq2

ofac=1.0 !over sampling

Mstar=1.0  !Mass and radius of star - to be read in
if(iargc().ge.5)then   !check if commandline wants a different value
   call getarg(5,cline)
   read(cline,*) fread
   if(fread.gt.0.0d0) Mstar=fread
endif
Rstar=1.0
if(iargc().ge.6)then   !check if commandline wants a different value
   call getarg(6,cline)
   read(cline,*) fread
   if(fread.gt.0.0d0) Rstar=fread
endif
write(0,*) "Mass,Radius of host: ",Mstar,Rstar

write(0,*) "Calc number of steps"
call calcnsteps(nstep,Mstar,Rstar,freq1,freq2,nyq,ofac,npt) !number of f's
write(0,*) "nstep: ",nstep

!allocate array space for periods (p) and logl (log-likelihood)
allocate(p(nstep),pow(nstep))

call pdm(nstep,p,pow,npt,time,flux,ferr,freq1,freq2,nyq,ofac,Mstar,     &
 Rstar,bper,bpow)

epo=0.0
depth=0.0
pmean=0.0
std=1.0
qtran=1.0

nplot=0 !disable plotting by default
if(iargc().ge.2)then
   call getarg(2,cline)
   read(cline,*) nplot
endif
if(nplot.ge.1)then
   call makeplot(nplot,nstep,p,pow,obsfile,epo,bper,bpow,depth,pmean, &
      std,qtran,npt,time,flux)
endif

end program pdmsearch
