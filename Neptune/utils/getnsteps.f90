subroutine getnsteps(fitsfile,nsteps,npixmax)
use precision
implicit none
integer :: status,readwrite,unitfits,blocksize,nhuds,hud,hudtype,       &
   nrecmax,nkeys,j,nrows,ncols,datacode,repeat,width,nsteps,npixmax
character(80) :: fitsfile
character(80), allocatable, dimension(:) :: header
character(8) :: imagetype,extname

nrecmax=700 !maximum number of header entries.
allocate(header(nrecmax)) !readheader assumes nrecmax=700

!first we open the fits file
!cfitsio wants this initialized
status=0

!setting to zero makes fits file readwrite
readwrite=0

!gets an unused unit number to open fits file
call ftgiou(unitfits,status)
!open the fits file
call ftopen(unitfits,fitsfile,readwrite,blocksize,status)
!get the number of the extensions (huds) in the image
call ftthdu(unitfits,nhuds,status)

!write(6,*) "NHUDS: ",nhuds,status
do hud=1,nhuds
!This command moves us to the next HUD
   call ftmahd(unitfits,hud,hudtype,status)
   if(status.eq.0) then
      call readheader(unitfits,status,header,nkeys)
!      write(6,*) "HUD: ",hud,nkeys
      do j=1,nkeys
!         write(6,*) j,header(j)(1:8)
         if(header(j)(1:8).eq."XTENSION") then
            read(header(j)(12:19),503) imagetype
 503        format(A8)
!            write(6,*) header(j)(1:8),header(j)(12:19)
         endif
         if(header(j)(1:8).eq."EXTNAME") then
            read(header(j)(12:19),503) extname
!            write(6,*) header(j)(1:8),header(j)(12:19)
!            write(6,*) "extname:",extname
         endif
      enddo
      if((imagetype.eq."BINTABLE").and.(extname.eq."TARGETTA")) then
         call FTGNRW(unitfits,nrows, status)
         call FTGNCL(unitfits,ncols, status)
!         write(6,*) "R,C:",nrows,ncols
         nsteps=nrows
         npixmax=0
         do j=1,ncols
            call ftgtcl(unitfits,j,datacode,repeat,width,status)
            npixmax=max(npixmax,repeat)
         enddo
      endif
   endif
enddo

!close the fits file
call ftclos(unitfits,status)
!return unit number to unused list.
call ftfiou(unitfits,status)

return
end subroutine getnsteps
