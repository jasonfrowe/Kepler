      program q1q12nlist
      implicit none
      integer nmax,nbits,nunit,kid,ndtype,npt,i,j,flag
      parameter(nmax=10000,nbits=15)
      integer lineout(nmax,nbits)
      real koi
      character*80 filename

C     bits
C 1   KOI
C 2   KID
C 3   detrend type
C 4   detrend length
C 5   nplanet
C6-15 planet exist?

      do 5 i=1,nmax
         do 6 j=1,nbits
            lineout(i,j)=0
 6       continue
 5    continue

      filename="koiprops.list"
c      filename="test.list"
      nunit=10
      npt=0
      open(unit=nunit,file=filename,status='old',err=901)
 11   read(nunit,*,end=12) koi,kid,ndtype
c         write(0,*) int(100*(koi-int(koi))+0.1)
         flag=0
         do 10 i=1,npt
            if(int(koi).eq.lineout(i,1))then
               lineout(i,5)=lineout(i,5)+1
               lineout(i,5+int(100*(koi-int(koi))+0.1))=1
               flag=1
            endif
 10      continue
         if(flag.eq.0)then
            npt=npt+1
            lineout(npt,1)=int(koi)
            lineout(npt,2)=kid
            lineout(npt,3)=ndtype
            lineout(npt,4)=2
            lineout(npt,5)=1
            lineout(npt,5+int(100*(koi-int(koi))+0.1))=1
            lineout(npt,6)=1
         endif
      goto 11
 12   continue
      close(nunit)

      do 13 i=1,npt
         write(6,501) (lineout(i,j),j=1,nbits)
 501     format(I4,1X,I9,1X,I1,1X,I2,1X,11(I1,1X))
 13   continue

      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 999  end
