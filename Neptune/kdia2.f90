program kdia2
use precision
implicit none
integer :: iargc
character(80) :: Imagename,Refname,Outname

if(iargc().lt.3)then
   write(0,*) "Usage: kdia <Image> <Reference> <Subtracted> [Flat]"
   write(0,*) " <Image> : FITS image for subtraction"
   write(0,*) " <Reference> : FITS reference image to match"
   write(0,*) " <Subtracted>: FITS output of subtracted image"
   stop
endif
call getarg(1,Imagename) !get Image filename
call getarg(2,Refname) !get Reference filename
call getarg(3,Outname) !get output filename (subtracted image)

end program kdia2
