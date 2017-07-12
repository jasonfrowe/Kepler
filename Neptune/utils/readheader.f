CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readheader(unitfits,status,header,nkeys)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     reads in the fits header
C     unitsfits is the unit number assigned to the fits file
C     status is from openfits and is used by cfitsio
      implicit none
      integer nrecmax
      parameter(nrecmax=700)
      integer status,unitfits,nkeys,nspace,i
      character*80 record,header(nrecmax)

      nkeys=0

C     get number of headers in image
      call ftghsp(unitfits,nkeys,nspace,status)

      if (nkeys.gt.nrecmax) then
         write(6,*) "WARNING: nkeys is more than nrecmax!!!!"
         write(6,*) "nkeys, nrecmax", nkeys, nrecmax
        pause
      endif

C     read in each header and move it into the master list.
      do 10 i=1,nkeys
         call ftgrec(unitfits,i,record,status)
         header(i)=record
 10   continue

      return
      end
