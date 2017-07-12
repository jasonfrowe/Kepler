module precision
   implicit none
!   type fitlookuptable
!      integer :: fittype
!      integer :: fitparm
!   end type fitlookuptable
   integer, parameter :: double  = 8 !double precision
   real(double), parameter :: Msun=1.98911d30, Yr=31557600.0d0
   real(double), parameter :: MU=Msun,TU=Yr, &
      LU=(MU*TU*TU*6.674d-11)**(1.0/3.0), G=1.0d0, Rsun=696265000.0, &
      AU=149597871000.0,Mearth=5.97219e24
   real(double) nbod
   real(double),allocatable, dimension(:) :: m!,merr
!   real(double),allocatable, dimension(:,:) :: mserr
end module precision
