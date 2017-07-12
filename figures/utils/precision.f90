module precision
   implicit none
!   type fitlookuptable
!      integer :: fittype
!      integer :: fitparm
!   end type fitlookuptable
   integer, parameter :: double  = 8 !double precision
   real(double), parameter :: Msun=1.9891d30, Yr=31557600.0d0,          &
    Rsun=696265000.0,AU=1.4959787d11,Mearth=5.97219e24,day=86400.0,     &
    Tsun=5781.6d0,G=6.674d-11
end module precision
