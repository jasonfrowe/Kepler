module fittermod
   use precision, only: double
   implicit none
   integer, pointer :: ngsol2
   integer, dimension(:), pointer :: naxes2,numnei2
   integer, dimension(:,:), pointer :: starmap2,nei2
   real(double), pointer :: sat2
   real(double), dimension(:,:), pointer :: Image2
end module fittermod
