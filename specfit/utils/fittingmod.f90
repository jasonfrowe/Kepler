module fittingmod
   use precision, only: double
   implicit none
   integer, pointer :: nmodel2
   real(double), dimension(:), pointer :: wv2,flux2,wmod2,fmod2,serr2
end module fittingmod
