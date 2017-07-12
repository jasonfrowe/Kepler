subroutine fitsread(fitsfile,nstep,npix,tpix,pix,crval,nkeysmax,        &
 nkeysout,headerout,ix,iy)
use precision
implicit none
integer :: status,readwrite,unitfits,blocksize,nhuds,hud,hudtype,       &
   nrecmax,nkeys,j,nrows,ncols,datacode,repeat,width,nstep,npix,k,      &
   nread,nkeysout,nkeysmax,ix,iy
integer, dimension(2) :: crval
real(double) :: ts,tpix
real(double), dimension(:) :: pix
logical anyf
character(80) :: fitsfile
character(80), allocatable, dimension(:) :: header
character(80), dimension(:) :: headerout
character(8) :: imagetype,extname

nrecmax=nkeysmax !maximum number of header entries.
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
!         write(6,*) j,header(j)(1:30)
!         if(header(j)(1:8).eq."NREADOUT") then
!!            write(6,*) header(j)(1:8),header(j)(12:19)
!            read(header(j)(12:30),*) nread
!            write(0,*) "NREAD:",nread
!         endif
         if(header(j)(1:8).eq."CRVAL1P") then
!            write(6,*) header(j)(1:8),header(j)(12:19)
            read(header(j)(12:30),*) crval(1)
         endif
         if(header(j)(1:8).eq."CRVAL2P") then
!            write(6,*) header(j)(1:8),header(j)(12:19)
            read(header(j)(12:30),*) crval(2)
         endif
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
         if(header(j)(1:8).eq."TDIM5") then
!            write(6,*) header(j)(1:8),header(j)(12:19)
            read(header(j)(13:30),*) ix,iy
!            write(0,*) ix,iy
!            write(0,'(A18)') header(j)(12:30)
         endif
      enddo
!      write(6,*) "extname:",extname,imagetype

!      if((imagetype.eq."IMAGE").and.(extname.eq."APERTURE")) then
!         write(0,*) "hello!"
!         read(5,*)
!      endif

      if(extname.eq."TARGETTA")then
         nkeysout=nkeys  !copy header info for export.
         headerout=header
      endif
      if((imagetype.eq."BINTABLE").and.(extname.eq."TARGETTA")) then
         call FTGNRW(unitfits,nrows, status)
         call FTGNCL(unitfits,ncols, status)
!         write(6,*) "R,C:",nrows,ncols
         do j=1,ncols
            call ftgtcl(unitfits,j,datacode,repeat,width,status)
!            write(6,*) "cc:",j,datacode,repeat,width,status
            if(j.eq.1)then
               call ftgcvd(unitfits,j,nstep,1,1," ",ts,anyf,status)
               tpix=ts
!               write(6,*) nstep,ts
            endif
            ! 5=corrected, 4=raw
            if(j.eq.5)then
               npix=repeat
!               write(0,*) "npix: ",npix
               call ftgcvd(unitfits,j,nstep,1,repeat," ",pix,anyf,status)
            endif
         enddo
      endif
   endif
enddo

!close the fits file
call ftclos(unitfits,status)
!return unit number to unused list.
call ftfiou(unitfits,status)

return
end subroutine fitsread
