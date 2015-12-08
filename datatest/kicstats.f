      program kicstats
      implicit none
      integer nmax,kid,Teff,nunit,i,nQ,nt,j,nflag
      parameter(nmax=300000)
      integer kidobs(nmax)
      real kepmag,logg,FeH
      character*80 filename,kicdatabase,line,dumc

      kicdatabase="kic_colour.format.txt"

      nQ=10 !number of quarters

      nt=0 !count number of targets observed.
      nunit=10
      do 10 i=1,nQ
         if(i.lt.10)then
            write(filename,500) "targetq",i,".list"
 500        format(A7,I1,A5)
         else
            write(filename,501) "targetq",i,".list"
 501        format(A7,I2,A5)
         endif

         write(0,*) filename
         open(unit=nunit,file=filename,status='old',err=901)
 11      read(nunit,*,end=12) line
            read(line,502,err=11) dumc,kid
502         format(A3,I8)
            nflag=0
            do 13 j=1,nt
               if(kidobs(j).eq.kid) nflag=1
 13         continue
            if(nflag.eq.0)then
               nt=nt+1
               kidobs(nt)=kid
            endif
         goto 11
 12      close(nunit)

 10   continue

      write(0,*) "number of targets ever observed: ",nt

      do 14 i=1,nt
         write(6,*) kidobs(i)
 14   continue



      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 902  write(0,*) "error reading: ",line
      goto 999
 999  end
