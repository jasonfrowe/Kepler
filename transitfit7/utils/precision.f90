module precision
   implicit none
   type fitlookuptable
      integer :: fittype
      integer :: fitparm
   end type fitlookuptable
   integer, parameter :: double  = 8 
   real(double), parameter :: Yr=31557600.0
   real(double), parameter :: Rsun=696265000.0
   real(double), parameter :: Mearth=5.97219e24
   real(double), parameter :: day=86400.0
   real(double), parameter :: PI = 3.141592653589793
   !max integration time for n-body steps (seconds)
   real(double), parameter :: maxintg = 2000.0
   real(double), parameter :: maxintg_nt = 10.0
   
   !Mercury
   real(double), parameter :: Gr=6.6719842229637d-11 
   real(double), parameter :: K2=2.959122082855911d-04
   real(double), parameter :: Msun=1.9891e30
   real(double), parameter :: AU=1.4959787d11
   real(double), parameter :: MU=Msun,TU=day,LU=AU


   !BS-Step (modern)
   !real(double), parameter :: Gr=6.67408d-11 !standard value
   !real(double), parameter :: K2=1.0d0 !not used for BS-Step calcs, so set to 1.
   !real(double), parameter :: Msun=1.9891d30
   !real(double), parameter :: AU=1.4959787d11 
   !real(double), parameter :: MU=Msun,TU=Yr,LU=(MU*TU*TU*Gr)**(1.0/3.0)
   !real(double), parameter :: G=1.0d0

   !shared vars for BS-Step routines via pointers
   !integer, pointer :: nbodies2
   !real(double), dimension(:), pointer :: m2
end module precision
