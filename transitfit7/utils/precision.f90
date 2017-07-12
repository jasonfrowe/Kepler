module precision
   implicit none
!   type fitlookuptable
!      integer :: fittype
!      integer :: fitparm
!   end type fitlookuptable
   integer, parameter :: double  = 8 !double precision
   real(double), parameter :: Msun=1.9891d30, Yr=31557600.0d0,G=1.0d0, &
    Rsun=696265000.0,AU=1.4959787d11,Mearth=5.97219e24,day=86400.0
   real(double), parameter :: MU=Msun,TU=day,LU=AU
   integer nbod
   real(double),allocatable, dimension(:) :: m!,merr
!   real(double),allocatable, dimension(:,:) :: mserr
end module precision
