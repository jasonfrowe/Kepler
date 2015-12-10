program rhostariso
use precision
implicit none
integer i,nyear
real(double) :: teff,teffsig,logg,loggsig,feh,fehsig,rhoi
real(double), allocatable, dimension(:) :: rhoierr,years
character(80) filename

interface
   subroutine getinitmodel(teff,logg,feh,rhoi)
      use precision
      implicit none
      real(double), intent(in) :: teff,logg,feh,rhoi
   end subroutine getinitmodel
end interface

teff=3841.0
teffsig=50.0
logg=0.0
loggsig=0.0 !if an error is zero or below, that variable is not fit.
feh=-0.18
fehsig=0.10
rhoi=10.6 !g/cc
allocate(rhoierr(9))
rhoierr(1)=-rhoi
rhoierr(2)=-6.7
rhoierr(3)=-5.1
rhoierr(4)=-3.1
rhoierr(5)=0.0d0
rhoierr(6)= 2.4
rhoierr(7)= 6.2
rhoierr(8)=10.1
rhoierr(9)=rhoierr(8)*10.0d0

!okay, we start by scanning through all the isochrones to find the best
!match to teff, log(g) and [Fe/H]
call getinitmodel(teff,logg,feh,rhoi)

!There are 38 ages in the Dartmouth isochrones
nyear=38
allocate(years(nyear))


end program rhostariso

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine getinitmodel(teff,logg,feh,rhoi)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
use isochrones
implicit none
integer :: i,nunit,filestatus
real(double) :: teff,logg,feh,rhoi

do i = 1,niso

   nunit=10
   open(unit=nunit,file=isofile(i),iostat=filestatus,status='old')

   if(filestatus>0)then
   write(0,*) "Cannot open ",photfile
   stop
endif


return
end subroutine getinitmodel

