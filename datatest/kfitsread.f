      program kfitsread
C     Reads in FITS files with Kepler data from MAST
C     Jason Rowe - jasonfrowe@gmail.com
      implicit none
      integer unitfits,status,iargc,readwrite,blocksize,nhuds,hudtype,
     .  hud,nrecmax,j,nkeys,nmax,npt,k,i,nfile,qflagcomb
      parameter(nrecmax=700,nmax=2000000)
      integer qflag(nmax),ns(nmax),qfbit(21)
      double precision time(nmax),flux(nmax),ferr(nmax),median
      character*160 filename
      character*80 header(nrecmax)
      character*8 imagetype,extname
      
      if(iargc().lt.1) goto 901
      nfile=iargc()
      
      do 14 i=1,nfile
      
        call getarg(i,filename)
        write(0,*) i,filename
      
C       first we open the fits file
C       cfitsio wants this initialized
        status=0
      
C       setting to zero makes fits file readwrite
        readwrite=0

C       gets an unused unit number to open fits file      
        call ftgiou(unitfits,status)
C       open the fits file
        call ftopen(unitfits,filename,readwrite,blocksize,status)      
C       get the number of the extensions (huds) in the image
        call ftthdu(unitfits,nhuds,status)
      
        if(status.ne.0) then
            write(6,*) "whoops...",status
        endif
      
        do 10 hud=1,nhuds
      
c           This command moves us to the next HUD
            call ftmahd(unitfits,hud,hudtype,status)
        
            if(status.eq.0) then
C               Read in all the header information
                call readheader(unitfits,status,header,nkeys)
                do 11 j=1,nkeys
                    if(header(j)(1:8).eq."XTENSION") then
                        read(header(j)(12:19),503) imagetype
 503                    format(A8)
c                    write(6,*) header(j)(1:8),header(j)(12:19)
c                     write(0,*) imagetype
                    endif
                    if(header(j)(1:8).eq."EXTNAME") then
                        read(header(j)(12:19),503) extname
c                    write(6,*) header(j)(1:8),header(j)(12:19)
c                  write(6,*) "extname:",extname
                    endif
11              continue

                if((imagetype.eq."BINTABLE").and.
     .            ((extname.eq."LIGHTCUR").or.
     .            (extname.eq."INJECTED"))) then
c                if((imagetype.eq."BINTABLE").and.
c     .            (extname.eq."TARGETTA")) then
                    call rdfitbtab(unitfits,status,nmax,npt,time,flux,
     .                ferr,qflag)
                endif

            endif
        
c        write(0,*) imagetype,extname

C       check for any FITSIO error codes.
            if(status.ne.0) then
                write(6,*) "whoops check FITSIO error code:",status
                write(6,*) "filename:",filename,"[",HUD,"]"
                pause
                status=0
            endif
      
 10     continue     

C     Exclude bad data
        k=0
        do 13 j=1,npt
            call getbit(qflag(j),qfbit,qflagcomb)
            if((qflagcomb.eq.0).and.(isnan(flux(j)).eqv..false.))then
                  k=k+1
                  time(k)=time(j)
                  flux(k)=flux(j)
                  ferr(k)=ferr(j)
                  qflag(k)=qflag(j)
            endif
 13     continue
        npt=k
 
C       Calculate the median flux
        call rqsort(npt,flux,ns)
        median=flux(ns(npt/2+1))
c        median=1.0
c        write(0,*) "Median: ",median
  
C     Output the data 
        do 12 j=1,npt
            if((qflagcomb.eq.0).and.(flux(j)/median.gt.-0.9999))then
                write(6,*) time(j)+54833.0-0.5,flux(j)/median-1.0d0,
     .              ferr(j)/median
c               write(6,*) time(j)+54833.0-0.5,flux(j),ferr(j)
            endif
 12     continue
 500    format(F17.11,1X,F17.11,1X,F17.11)

C       close the fits file
        call ftclos(unitfits,status)
C       return unit number to unused list.
        call ftfiou(unitfits,status)
      
c       write(0,*) unitfits,status,nhuds
      
 14   continue
      
      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getbit(qflag,qfbit,qflagcomb)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer qflag,qfbit(21),i,qflagcomb
      character(len=21) cbin

      !write(0,*) "qflag: ",qflag
      write(cbin,500) qflag
!      write(0,*) cbin(1:21)
 500  format(B21.21)
      do 10 i=1,21
!         write(0,*) 22-i,cbin(i:i),i
         read(cbin(i:i),*) qfbit(22-i)
 10   continue

      qflagcomb=0
      if (qfbit(1).eq.1) qflagcomb=1
      if (qfbit(2).eq.1) qflagcomb=1
      if (qfbit(3).eq.1) qflagcomb=1
      if (qfbit(4).eq.1) qflagcomb=1
      if (qfbit(6).eq.1) qflagcomb=1
      if (qfbit(7).eq.1) qflagcomb=1
      if (qfbit(9).eq.1) qflagcomb=1
      if (qfbit(13).eq.1) qflagcomb=1
      if (qfbit(15).eq.1) qflagcomb=1
      if (qfbit(16).eq.1) qflagcomb=1
      if (qfbit(17).eq.1) qflagcomb=1

!      write(6,501) (qfbit(i),i=1,21)
 501  format(21I1)

!      read(5,*)

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rdfitbtab(unitfits,status,nmax,npt,time,flux,ferr,
     .  qflag)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     reads binary tables in FITS file and get the lightcurve
      implicit none
      integer unitfits,status,nrows,ncols,i,datacode,repeat,width,j,
     .  nmax,npt,qflag(nmax)
      double precision time(nmax),ts,flux(nmax),ferr(nmax)
      logical anyf
      
C     get the number of rows and columns in the table.
      call FTGNRW(unitfits,nrows, status)
      call FTGNCL(unitfits,ncols, status)
      
c      write(0,*) "nn:",nrows,ncols
      npt=nrows
      
C     loop over all the columns (have to do it this way with FITSIO
      do 10 i=1,ncols
C     see what kind of data string we are dealing with.
C     in the secondary target table, it should only be characters and
C     integers.
         call ftgtcl(unitfits,i,datacode,repeat,width,status)
c         write(6,*) "cc:",i,datacode,repeat,width,status
         if(i.eq.1)then
            do 11 j=1,nrows
C     read in the character string from the table.
                call ftgcvd(unitfits,i,j,1,1," ",ts,anyf,status)
                time(j)=ts
c                write(6,*) j,time(j)
 11         continue
         endif
C        8-PDC-Map,4-SAP_FLUX
         if(i.eq.8)then
            do 12 j=1,nrows
C     read in the character string from the table.
                call ftgcvd(unitfits,i,j,1,1," ",ts,anyf,status)
                flux(j)=ts
 12         continue
         endif
         if(i.eq.9)then
            do 13 j=1,nrows
C     read in the character string from the table.
                call ftgcvd(unitfits,i,j,1,1," ",ts,anyf,status)
                ferr(j)=ts
 13         continue
         endif
         if(i.eq.10)then
            do 14 j=1,nrows
C     read in the character string from the table.
                call ftgcvd(unitfits,i,j,1,1," ",ts,anyf,status)
                qflag(j)=int(ts)
 14         continue
         endif

 10   continue
      
      return
      end
 
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
      
c**********************************************************************
      subroutine rqsort(n,a,p)
c======================================================================
c     Return integer array p which indexes array a in increasing order.
c     Array a is not disturbed.  The Quicksort algorithm is used.
c
c     B. G. Knapp, 86/12/23
c
c     Reference: N. Wirth, Algorithms and Data Structures,
c     Prentice-Hall, 1986
c======================================================================
      implicit none

c     Input:
      integer   n
      real*8      a(n)

c     Output:
      integer   p(n)

c     Constants
      integer   LGN, Q
      parameter (LGN=32, Q=11)
c        (LGN = log base 2 of maximum n;
c         Q = smallest subfile to use quicksort on)

c     Local:
      real*8      x
      integer   stackl(LGN),stackr(LGN),s,t,l,m,r,i,j

c     Initialize the stack
      stackl(1)=1
      stackr(1)=n
      s=1

c     Initialize the pointer array
      do 1 i=1,n
         p(i)=i
    1 continue

    2 if (s.gt.0) then
         l=stackl(s)
         r=stackr(s)
         s=s-1

    3    if ((r-l).lt.Q) then

c           Use straight insertion
            do 6 i=l+1,r
               t = p(i)
               x = a(t)
               do 4 j=i-1,l,-1
                  if (a(p(j)).le.x) goto 5
                  p(j+1) = p(j)
    4          continue
               j=l-1
    5          p(j+1) = t
    6       continue
         else

c           Use quicksort, with pivot as median of a(l), a(m), a(r)
            m=(l+r)/2
            t=p(m)
            if (a(t).lt.a(p(l))) then
               p(m)=p(l)
               p(l)=t
               t=p(m)
            endif
            if (a(t).gt.a(p(r))) then
               p(m)=p(r)
               p(r)=t
               t=p(m)
               if (a(t).lt.a(p(l))) then
                  p(m)=p(l)
                  p(l)=t
                  t=p(m)
               endif
            endif

c           Partition
            x=a(t)
            i=l+1
            j=r-1
    7       if (i.le.j) then
    8          if (a(p(i)).lt.x) then
                  i=i+1
                  goto 8
               endif
    9          if (x.lt.a(p(j))) then
                  j=j-1
                  goto 9
               endif
               if (i.le.j) then
                  t=p(i)
                  p(i)=p(j)
                  p(j)=t
                  i=i+1
                  j=j-1
               endif
               goto 7
            endif

c           Stack the larger subfile
            s=s+1
            if ((j-l).gt.(r-i)) then
               stackl(s)=l
               stackr(s)=j
               l=i
            else
               stackl(s)=i
               stackr(s)=r
               r=j
            endif
            goto 3
         endif
         goto 2
      endif
      return
      end
