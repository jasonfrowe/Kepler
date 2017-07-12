subroutine writefits(nxmax,nymax,parray,bpix,tavg,nkeys,header,nstep)
use precision
implicit none
integer :: nxmax,nymax,nkeys,nstep,status,blocksize,bitpix,naxis,funit, &
   npixels,group,firstpix,nbuf,i,j,nbuffer,npt,ii,jj
integer, dimension(2) :: naxes
integer, dimension(4) :: nr
real(double), allocatable, dimension(:) :: buffer
real(double) :: bpix,tavg
real(double), dimension(:,:) :: parray
real(double), allocatable, dimension(:,:) :: oarray
character(80) :: fileout,record
character(80), dimension(:) :: header
logical simple,extend

!create filename for output
if(nstep.lt.10)then
   write(fileout,501) 'n00000',nstep,'.fits'
   501 format(A6,I1,A5)
elseif(nstep.lt.100)then
   write(fileout,502) 'n0000',nstep,'.fits'
   502 format(A5,I2,A5)
elseif(nstep.lt.1000)then
   write(fileout,503) 'n000',nstep,'.fits'
   503 format(A4,I3,A5)
elseif(nstep.lt.10000)then
   write(fileout,504) 'n00',nstep,'.fits'
   504 format(A3,I4,A5)
elseif(nstep.lt.100000)then
   write(fileout,505) 'n0',nstep,'.fits'
   505 format(A2,I5,A5)
else
   write(fileout,506) 'n',nstep,'.fits'
   506 format(A1,I6,A5)
endif
write(0,'(A15,A12)') "Writing file.. ", fileout

!assemble pixels into compact square.
npt=0
do i=1,nxmax
   do j=1,nymax
      if(parray(i,j).lt.bpix)then
         npt=npt+1
         if(npt.eq.1)then
            nr(1)=i
            nr(2)=i
            nr(3)=j
            nr(4)=j
         else
            nr(1)=min(nr(1),i)
            nr(2)=max(nr(2),i)
            nr(3)=min(nr(3),j)
            nr(4)=max(nr(4),j)
         endif
      endif
   enddo
enddo

write(0,*) "Offset: ",nr(1)-1,nr(3)-1

naxes(1)=nr(2)-nr(1)+1
naxes(2)=nr(4)-nr(3)+1
allocate(oarray(naxes(1),naxes(2)))
ii=0
do i=nr(1),nr(2)
   ii=ii+1
   jj=0
   do j=nr(3),nr(4)
      jj=jj+1
      oarray(ii,jj)=parray(i,j)
   enddo
enddo

status=0
!if file already exists.. delete it.
call deletefile(fileout,status)
!get a unit number
call ftgiou(funit,status)
!Create the new empty FITS file.  The blocksize parameter is a
!historical artifact and the value is ignored by FITSIO.
blocksize=1
call ftinit(funit,fileout,blocksize,status)

!Initialize parameters about the FITS image.
!BITPIX = 16 means that the image pixels will consist of 16-bit
!integers.  The size of the image is given by the NAXES values.
!The EXTEND = TRUE parameter indicates that the FITS file
!may contain extensions following the primary array.
simple=.true.
bitpix=-32
naxis=2
extend=.true.

!Write the required header keywords to the file
call ftphpr(funit,simple,bitpix,naxis,naxes,0,1,extend,status)

!Write the array to the FITS file.
npixels=naxes(1)*naxes(2)
group=1
firstpix=1
nbuf=naxes(1)
j=0

allocate(buffer(nbuf))
do while (npixels.gt.0)
!read in 1 column at a time
   nbuffer=min(nbuf,npixels)

   j=j+1
!find max and min values
   do i=1,nbuffer
      buffer(i)=oarray(i,j)
   enddo

   call ftpprd(funit,group,firstpix,nbuffer,buffer,status)

!update pointers and counters

   npixels=npixels-nbuffer
   firstpix=firstpix+nbuffer

enddo

do i=6,nkeys
   record=header(i)
   if(record(1:8).eq.'TELESCOP')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'INSTRUME')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'OBJECT')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'KEPLERID')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'RA_OBJ')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'DEC_OBJ')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'EQUINOX')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'EXPOSURE')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'TIMEREF')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'TASSIGN')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'TIMESYS')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'BJDREFI')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'BJDREFF')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'TIMEUNIT')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'TELAPSE')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'LIVETIME')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'TSTART')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'TSTOP')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'LC_START')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'LC_END')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'DEADC')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'TIMEPIXR')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'TIERRELA')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'TIERABSO')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'READTIME')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'FRAMETIM')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'NUM_FRM')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'TIMEDEL')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'DATE-OBS')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'DATE-END')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'BACKAPP')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'DEADAPP')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'VIGNAPP')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'GAIN')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'READNOIS')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'NREADOUT')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'TIMSLICE')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'MEANBLCK')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'LCFXDOFF')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'SCFXDOFF')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'CROWDSAP')then
      call ftprec(funit,record,status)
   elseif(record(1:8).eq.'FLFRCSAP')then
      call ftprec(funit,record,status)
   endif
!   write(6,*) "ftprec:",status
enddo

!write(6,*) "ftprec:",status
write(record,'(A8,A3,F10.1)') 'DATAMAX ','=  ',bpix
!write(6,'(a80)') record
call ftprec(funit,record,status)
write(record,'(A8,A3,F10.1)') 'DATAMIN ','=  ',-10000.0
!write(6,'(a80)') record
call ftprec(funit,record,status)
write(record,'(A8,A3,F14.10)') 'TIMEOBS ','=  ',tavg
!write(6,'(a80)') record
call ftprec(funit,record,status)
!write(6,*) "ftprec:",status


!close fits file
call ftclos(funit,status)
call ftfiou(funit,status)

return
end subroutine writefits
