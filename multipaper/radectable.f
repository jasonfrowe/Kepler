      program radectable
      implicit none
      integer nmax,nunit,i,nkid,kidcddp,nloop,kidspars,kidchar,nkoimax,
     .   nkoi,nloop2,j,nfp,nc,nsn
      parameter(nmax=6570711,nkoimax=6000)
      integer flag(nmax),kid(nmax),kflag(nmax),intkoi(nmax),fp(nkoimax),
     .   rfp,radeckoi,numkoi(nmax),tapo(nmax),napo
      real qobs(nmax),cddp,logg(nmax),dumr,readlogg,loggcut,rper,rrprs,
     .   rb,rsn,readkoi,koi(nkoimax),per(nkoimax),b(nkoimax),
     .   rprs(nkoimax),mfp(nkoimax),fpkoi(nkoimax),ckoi(nkoimax),
     .   ra(nmax),dec(nmax),rra,rdec,akoi(nkoimax),pdone,sn(nkoimax),
     .   koisn(nkoimax)
      character*80 kidtable,cddptable,spars,dumc,koichar,mpars,fpfile,
     .   cfile,kicradec,apofile,cout,koicharsn

      loggcut=3.5

      nunit=10

      kidtable="kidq1q8.txt"
      open(unit=nunit,file=kidtable,status='old')
      i=1
 10   read(nunit,*,end=11) kid(i),qobs(i)
        i=i+1
      goto 10
 11   continue
      nkid=i-1
      close(nunit)

      do 12 i=1,nkid
         flag(i)=0
         logg(i)=0.0
         kflag(i)=0
         intkoi(i)=0.0
         numkoi(i)=0
         tapo(i)=0
 12   continue

      spars="ktc_consol_v2.sp.dat"
      open(unit=nunit,file=spars,status='old')
      read(nunit,*) dumc
 16   read(nunit,*,end=17) kidspars,dumr,dumr,readlogg
         nloop=0
         i=0
         do while(nloop.eq.0)
            i=i+1
            if(i.gt.nkid) nloop=2
            if((kidspars.eq.kid(i)).and.(nloop.eq.0))then
               nloop=1
               logg(i)=readlogg
            endif
         enddo
      goto 16
 17   continue
      close(nunit)

      cddptable="cdppq1q103h.dat"
      open(unit=nunit,file=cddptable,status='old')
 13   read(nunit,*,end=14) kidcddp,cddp
         nloop=0
         i=0
         do while(nloop.eq.0)
            i=i+1
            if(i.gt.nkid) nloop=2
            if((kidcddp.eq.kid(i)).and.(nloop.eq.0))then
               nloop=1
               if((cddp.lt.160.0).and.(qobs(i).gt.4).and.
     .              (logg(i).gt.3.5))then
                  flag(i)=1
               endif
            endif
         enddo
      goto 13
 14   continue
      close(nunit)

      mpars="multiprops.dat"
      open(unit=nunit,file=mpars,status='old')
      read(nunit,*) dumc
      i=1
 20   read(nunit,*,end=21) koi(i),dumr,dumr,per(i),dumr,b(i),dumr,
     .   rprs(i),dumr,dumr,dumr,dumr,dumr,dumr,dumr,dumr,dumr,dumr,
     .   dumr,dumr,dumr,dumr,dumr,mfp(i)
         i=i+1
      goto 20
 21   close(nunit)
      nkoi=i-1

      fpfile="Fergalfile.20130801.q1q8.awk.txt"
      open(unit=nunit,file=fpfile,status='old')
      i=1
 22   read(nunit,*,end=23) fpkoi(i),fp(i)
         i=i+1
         goto 22
 23   continue
      close(nunit)
      nfp=i-1

      cfile="PerT0collide.20130805.txt"
      open(unit=nunit,file=cfile,status='old')
      read(nunit,*) dumc
      i=1
 24   read(nunit,*,end=25) ckoi(i)
         i=i+1
         goto 24
 25   continue
      close(nunit)
      nc=i-1

      apofile="apolist.20130829.txt"
      open(unit=nunit,file=apofile,status='old')
      i=1
 28   read(nunit,*,end=29) akoi(i)
         i=i+1
         goto 28
 29   continue
      close(nunit)
      napo=i-1
      write(0,*) "napo: ",napo

      koicharsn="koi_characteristics.20130509.dat"
      open(unit=nunit,file=koicharsn,status='old')
      read(nunit,*) dumc
      read(nunit,*) dumc
      i=1
 30   read(nunit,*,end=31) koisn(i),(dumr,j=1,19),sn(i)
         i=i+1
      goto 30
 31   continue
      close(nunit)
      nsn=i-1

c      koichar="koi_characteristics.20130509.dat"
      koichar="koi_multi.20120926.dat"
      open(unit=nunit,file=koichar,status='old')
      read(nunit,*) dumc
      read(nunit,*) dumc
 18   read(nunit,*,end=19) readkoi,kidchar,dumr,dumr,dumr,rper,dumr,
     .     dumr,dumr,dumr,dumr,dumr,rrprs,rb,dumr,dumr,dumr,dumr,dumr,
     .     dumr,rsn
         if(readkoi.eq.floor(readkoi)) readkoi=readkoi+0.01
         if(readkoi.ge.3059.0) goto 18

         nloop2=0
         j=0
         do while(nloop2.eq.0)
            j=j+1
            if(j.gt.nsn) nloop2=2
            if((koisn(j).eq.readkoi).and.(nloop2.eq.0))then
               nloop2=1
               rsn=sn(j)
            endif
         enddo
         if(rsn.lt.7.1) goto 18 !ignore FAs.

         nloop=0
         i=0
         do while(nloop.eq.0)
            i=i+1
            if(i.gt.nkid) nloop=2
            if((kidchar.eq.kid(i)).and.(nloop.eq.0))then
               numkoi(i)=numkoi(i)+1

               nloop=1
               intkoi(i)=int(readkoi)

               nloop2=0
               j=0
               do while(nloop2.eq.0)
                  j=j+1
                  if(j.gt.nkoi) nloop2=2
                  if((koi(j).eq.readkoi).and.(nloop2.eq.0))then
                     nloop2=1
                     rper=per(j)
                     rb=b(j)
                     rrprs=rprs(j)
                  endif
               enddo
c               write(6,*) readkoi,rper,rb,rrprs,rsn,j
c               read(5,*)

               nloop2=0
               j=0
               do while(nloop2.eq.0)
                  j=j+1
                  if(j.gt.nfp) nloop2=2
                  if((fpkoi(j).eq.readkoi).and.(nloop2.eq.0))then
                     nloop2=1
                     rfp=fp(j)
                  endif
               enddo

               nloop2=0
               j=0
               do while(nloop2.eq.0)
                  j=j+1
                  if(j.gt.nc) nloop2=2
                  if((ckoi(j).eq.readkoi).and.(nloop2.eq.0))then
                     nloop2=1
                     rfp=2
                  endif
               enddo

               nloop2=0
               j=0
               do while(nloop2.eq.0)
                  j=j+1
                  if(j.gt.napo) nloop2=2
                  if((akoi(j).eq.readkoi).and.(nloop2.eq.0))then
                     nloop2=1
                     tapo(i)=tapo(i)+1
c                     write(0,*) akoi(j),tapo(i)
                  endif
               enddo

               if((rper.gt.1.6).and.(rb+rrprs.lt.0.98).and.
     .              (rsn.gt.10.0).and.(rfp.eq.0))then
                  kflag(i)=kflag(i)+1.0
               endif
            endif
         enddo
      goto 18
 19   continue
      close(nunit)

      write(0,*) "Starting RA/DEC Match"
      kicradec="kicradec.awk.txt"
      open(unit=nunit,file=kicradec,status='old')
      j=0
 26   read(nunit,*,end=27) radeckoi,rra,rdec
         j=j+1
         pdone=real(j)/6570711.0
         if(mod(int(100.0*pdone),1).eq.0)then
            write(cout,501) "radec:",pdone
 501        format(A6,1X,F6.3)
            call ovrwrt(cout,2)
         endif
         nloop=0
         i=0
         do while(nloop.eq.0)
            i=i+1
            if(i.gt.nkid) nloop=2
            if((radeckoi.eq.kid(i)).and.(nloop.eq.0))then
               nloop=1
               ra(i)=rra
               dec(i)=rdec
c               write(6,*) kid(i),ra(i),dec(i)
c               read(5,*)
            endif
         enddo
      goto 26
 27   continue

      do 15 i=1,nkid
         write(6,500) kid(i),ra(i),dec(i),flag(i),numkoi(i),kflag(i),
     .      tapo(i),intkoi(i)
 15   continue
 500  format(I8,2(1X,F9.6),4(1X,I1),1X,I4)

 999  end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ovrwrt (line, iwhich)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Taken from the DAOPhot code by Stetson.
C     This is a cool routine for writing lots of info to the screen
C     and keeping everything on the same line.
C
      character*(*) line
      character*79 output
      integer len
      if (iwhich .eq. 1) then
         write (0,1) line
    1    format (a)
      else if (iwhich .eq. 2) then
         if (len(line) .lt. 79) then
            output = ' '
            output = line
            write (0,2) output, char(13), char(13)
            write (0,2) output, char(13), char(13)
            write (0,2) output, char(13), char(13)
    2       format (a, 2a1, $)
         else
            write (0,2) line, char(13), char(13)
         end if
      else if (iwhich .eq. 3) then
         write (0,3) line
    3    format (a)
      else
         write (0,4) line, char(13), char(13)
    4    format (/a, 2a1, $)
         write (0,2) line, char(13), char(13)
      end if
      return
      end
