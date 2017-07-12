      program multistats
      implicit none
      integer nmax,i,nunit,n(6),j,k,nall,nTeq
      parameter(nmax=10000)
      integer npt,kid(nmax),pflag(nmax),flag(nmax)
      real koi(nmax),kepmag(nmax),Rp(nmax),T0(nmax),Per(nmax),
     .   Teff(nmax),logg(nmax),Rstar(nmax),FeH(nmax),adrs(nmax),
     .   rdr(nmax),b(nmax),asemi(nmax),Teq(nmax),gmag(nmax),RA(nmax),
     .   DEC(nmax),tdepth(nmax),SN(nmax),Tdur(nmax),chisq(nmax)
      character*80 filename,dumr

      nunit=10
      filename="koi_characteristics.20120703.csv"

      i=0
      open(unit=nunit,file=filename,status='old',err=901)
      read(nunit,*,err=902) dumr
      read(nunit,*,err=902) dumr
      i=1 !count number of lines read
 10   read(nunit,*,end=11,err=15) koi(i),kid(i),kepmag(i),Rp(i),T0(i),
     .   Per(i),Teff(i),logg(i),Rstar(i),FeH(i),pflag(i),adrs(i),rdr(i),
     .   b(i),asemi(i),Teq(i),gmag(i),RA(i),DEC(i),Tdepth(i),SN(i),
     .   Tdur(i),chisq(i),flag(i)
c         flag(i)=0
         if(SN(i).lt.10.0) flag(i)=1
         if(rdr(i)+b(i).gt.1.0) flag(i)=1
         if(koi(i).gt.3150.0) flag(i)=1
         if(flag(i).gt.1) write(0,*) koi(i)
c         if(teq(i).gt.300.0) flag(i)=1
c         if(kepmag(i).lt.14.2) flag(i)=1
c         write(0,*) i+2,koi(i),flag(i)
         i=i+1
      goto 10
 15   write(0,*) "Error on line: ",i+2
      goto 10
 11   continue
      close(nunit)
      npt=i-1
      write(0,*) "Number of lines read: ",npt

      do 13 i=1,6
         n(i)=0
 13   continue
      i=1
      do 12 while(i.lt.npt)
         nTeq=0
         if(flag(i).eq.0) then
            n(1)=n(1)+1 !count systems with at least 1
            if(Teq(i).lt.300.0) nTeq=1
         endif
         j=1
         k=1
         do 14 while(j.le.5) !cycle through 2 - 6 planet systems
            if(i+j.le.npt)then
c               if(int(koi(i)).eq.int(koi(i+j)))then
               if(koi(i+j)-koi(i).lt.0.5)then
                  if(flag(i+j).eq.0) then !check FP flag
                     if(Teq(i+j).lt.300.0) nTeq=nTeq+1
                     k=k+1
                     n(k)=n(k)+1
                  endif
                  j=j+1
               else
                  k=j
                  j=6 !break from loop, no more matches to consider
               endif
            else
               k=j
               j=6 !break from loop, no more KOI left
            endif
c            write(0,*) koi(i),k,j

 14      continue
         if(nTeq.gt.0) write(6,501) koi(i),flag(i),(n(j),j=1,6)
         i=i+k
 501  format(F7.2,1X,10(I5,1X))
c         read(5,*)
 12   continue

      nall=n(1)+n(2)+n(3)+n(4)+n(5)+n(6)

      write(6,500) nall,(n(i),i=1,6)
 500  format(10(I5,1X))

      goto 999
901   write(0,*) "Cannot open: ",filename
      goto 999
902   write(0,*) "Error reading on line: ",i+2
      goto 999
999   end
