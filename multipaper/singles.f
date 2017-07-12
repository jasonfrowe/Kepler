      program singles
      implicit none
      integer n,nmax,nunit,i,nEB,j,nex,k,nfp,nblarge
      parameter(nmax=10000,nex=17)
      integer kid(nmax),flag(nmax),pflag(nmax),nkoi(nmax),np(nmax),
     .   nkoi1(nmax),ebkid(nmax),nEBc,ebflag(nmax),ncol(nmax)
      real koi(nmax),kepmag(nmax),Rp(nmax),T0(nmax),Per(nmax),
     .   Teff(nmax),logg(nmax),Rstar(nmax),FeH(nmax),koi2,
     .   adrs(nmax),rdr(nmax),b(nmax),asemi(nmax),Teq(nmax),gmag(nmax),
     .   RA(nmax),DEC(nmax),Tdepth(nmax),SN(nmax),Tdur(nmax),SN2,
     .   chisq(nmax),ebper(nmax),dumr,ebpdepth(nmax),ebmorph(nmax),
     .   pars2(20),pexcl(nex),koifp,koicc
      character*80 filename,filenameEB,dumc,filename2,filenameFP,
     .   filenameCC
      data pexcl/375.01,422.01,435.02,490.02,771.01,1032.01,1096.01,
     .   1192.01,1206.01,1274.01,1356.01,1421.01,1463.01,1174.01,
     .   1208.01,1008.01,99.01/


      nunit=10
      filename="koi_multi.20120926.dat"
      filename2="koi_characteristics.20130509.dat"
      filenameEB="eb.20130517.list"
      filenameFP="Fergalfile.20130801.q1q8.awk.txt"
      filenameCC="PerT0collide.20130805.txt"

      open(unit=nunit,file=filename,status='old',err=901)
      read(nunit,*,err=902) dumc
      read(nunit,*,err=902) dumc
      i=1 !count number of lines read
 10   read(nunit,*,end=11,err=15) koi(i),kid(i),kepmag(i),Rp(i),T0(i),
     .   Per(i),Teff(i),logg(i),Rstar(i),FeH(i),pflag(i),adrs(i),rdr(i),
     .   b(i),asemi(i),Teq(i),gmag(i),RA(i),DEC(i),Tdepth(i),SN(i),
     .   Tdur(i),chisq(i),flag(i)!,cf(i),sflag(i)
         if(koi(i).gt.3058.01) goto 10
         if(floor(koi(i)).eq.koi(i)) koi(i)=koi(i)+0.01
         ncol(i)=0

c         if(cf(i).eq.2) goto 10
         i=i+1
      goto 10
 15   write(0,*) "Error on line: ",i+2
      goto 10
 11   continue
      close(nunit)
      n=i-1
      write(0,*) "Number of lines read: ",n

      open(unit=nunit,file=filenameFP,status='old',err=906)
 60   read(nunit,*,end=61) koifp,nfp
         do 62 i=1,n
            if(koifp.eq.koi(i))then
c               if(nfp.ne.flag(i)) write(6,*) koi(i),flag(i),nfp
               flag(i)=nfp
c               write(0,*) koi(i),flag(i)
            endif
 62      continue
      goto 60
 61   continue
      close(nunit)

      open(unit=nunit,file=filename2,status='old',err=905)
      read(nunit,*,err=904) dumc
      read(nunit,*,err=904) dumc
      i=1 !count number of lines read
 51   read(nunit,*,end=52,err=53) koi2,(pars2(j),j=1,20)
         if(koi2.gt.3058.01) goto 51
         if(floor(koi2).eq.koi2) koi2=koi2+0.01
         do 54 j=1,n
            if(koi(j).eq.koi2) then
               SN(j)=pars2(20) !SN replacement
c               write(0,*) koi(j),per(i),pars2(5)
               per(j)=pars2(5) !period

               do 55 k=1,nex
                  if(pexcl(k).eq.koi(j))then
c                     write(0,*) pexcl(k),koi(j)
                     per(j)=-1.0
                  endif
 55            continue

            endif
 54      continue


         i=i+1
      goto 51
 53   write(0,*) "Error on line: ",i+2
      goto 51
 52   continue
      close(nunit)
      n=i-1
      write(0,*) "Number of lines read: ",n

C     Read in Per/T0 collisions and mark them as FAs
      open(unit=nunit,file=filenameCC,status='old',err=907)
      read(nunit,*) dumc !first line is a comment in file
 63   read(nunit,*,end=64) koicc
         do 65 i=1,n
            if(koicc.eq.koi(i))then
c               write(0,*) "cc: ",koi(i)
               ncol(i)=1 !mark collisions
            endif
 65      continue
      goto 63
 64   continue
      close(nunit)

      open(unit=nunit,file=filenameEB,status='old',err=903)
      read(nunit,*) dumc
      i=1
 39   read(nunit,*,end=40) ebkid(i),ebper(i),dumr,ebpdepth(i),ebmorph(i)
         i=i+1
         goto 39
 40   continue
      nEB=i-1

C     No cuts..

      do 4 i=1,nmax
         nkoi1(i)=0
 4    continue

      do 20 i=1,n
           nkoi1(int(koi(i)))=nkoi1(int(koi(i)))+1
 20   continue

      do 21 i=1,6
         np(i)=0
 21   continue

      do 22 i=1,nmax
         if(nkoi1(i).gt.0)then
            np(nkoi1(i))=np(nkoi1(i))+1
         endif
 22   continue

      write(6,500) np(1),"No Cuts"
 500  format((I4,1X),A25,1X,I4)

CCCCCCCCCC FAs

      do 56 i=1,nmax
         nkoi(i)=0
 56    continue

      nblarge=0
      do 57 i=1,n
         if((sn(i).gt.7.1).and.(nkoi1(int(koi(i))).eq.1)) then
           nkoi(int(koi(i)))=nkoi(int(koi(i)))+1
           if(rdr(i)+b(i).gt.0.98) nblarge=nblarge+1
         endif
 57   continue

      do 58 i=1,6
         np(i)=0
 58   continue

      do 59 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 59   continue

      write(6,500) np(1),"FAs",nblarge

CCCCCCCCCC T0/Period collisions

      do 66 i=1,nmax
         nkoi(i)=0
 66    continue

      do 67 i=1,n
         if((sn(i).gt.7.1).and.(nkoi1(int(koi(i))).eq.1).and.
     .    (ncol(i).eq.0)) then
           nkoi(int(koi(i)))=nkoi(int(koi(i)))+1
         endif
 67   continue

      do 68 i=1,6
         np(i)=0
 68   continue

      do 69 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 69   continue

      write(6,500) np(1),"FAs,Col"

CCCCCCCCCC T0/Period collisions, P

      do 70 i=1,nmax
         nkoi(i)=0
 70    continue

      do 71 i=1,n
         if((sn(i).gt.7.1).and.(nkoi1(int(koi(i))).eq.1).and.
     .    (ncol(i).eq.0).and.(per(i).gt.1.6)) then
           nkoi(int(koi(i)))=nkoi(int(koi(i)))+1
         endif
 71   continue

      do 72 i=1,6
         np(i)=0
 72   continue

      do 73 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 73   continue

      write(6,500) np(1),"FAs,Col,P"


CCCCCCCCCC FPs

      do 23 i=1,nmax
         nkoi(i)=0
 23    continue

      nblarge=0
      do 24 i=1,n
         if((flag(i).lt.1).and.(nkoi1(int(koi(i))).eq.1).and.
     .    (sn(i).gt.7.1).and.(ncol(i).eq.0)) then
           nkoi(int(koi(i)))=nkoi(int(koi(i)))+1
         if((rdr(i)+b(i).gt.0.98).and.(rp(i).gt.15.0)) nblarge=nblarge+1
         endif
 24   continue

      do 25 i=1,6
         np(i)=0
 25   continue

      do 26 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 26   continue

      write(6,500) np(1),"FA,Col,FPs",nblarge

CCCCCCCCCC FPs,P

      do 74 i=1,nmax
         nkoi(i)=0
 74    continue

      do 75 i=1,n
         if((flag(i).lt.1).and.(nkoi1(int(koi(i))).eq.1).and.
     .    (sn(i).gt.7.1).and.(ncol(i).eq.0).and.(per(i).gt.1.6)) then
           nkoi(int(koi(i)))=nkoi(int(koi(i)))+1
         endif
 75   continue

      do 76 i=1,6
         np(i)=0
 76   continue

      do 77 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 77   continue

      write(6,500) np(1),"FA,Col,FPs,P"

CCCCCCCCCC FPs,SN

      do 27 i=1,nmax
         nkoi(i)=0
 27    continue

      do 28 i=1,n
         if((flag(i).lt.1).and.(nkoi1(int(koi(i))).eq.1).and.
     .    (sn(i).ge.10.0).and.(ncol(i).eq.0)) then
           nkoi(int(koi(i)))=nkoi(int(koi(i)))+1
         endif
 28   continue

      do 29 i=1,6
         np(i)=0
 29   continue

      do 30 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 30   continue

      write(6,500) np(1),"FA,Col,FPs,SN"

CCCCCCCCCC FPs,SN,b

      do 78 i=1,nmax
         nkoi(i)=0
 78    continue

      do 79 i=1,n
         if((flag(i).lt.1).and.(nkoi1(int(koi(i))).eq.1).and.
     .    (sn(i).ge.10.0).and.(ncol(i).eq.0).and.
     .    (rdr(i)+b(i).lt.0.98)) then
           nkoi(int(koi(i)))=nkoi(int(koi(i)))+1
         endif
 79   continue

      do 80 i=1,6
         np(i)=0
 80   continue

      do 81 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 81   continue

      write(6,500) np(1),"FA,Col,FPs,SN,b"

CCCCCCCCCC FPs,SN, Per

      do 31 i=1,nmax
         nkoi(i)=0
 31    continue

      do 32 i=1,n
         if((flag(i).lt.1).and.(nkoi1(int(koi(i))).eq.1).and.
     .    (sn(i).ge.10.0).and.(per(i).gt.1.6).and.(ncol(i).eq.0)) then
           nkoi(int(koi(i)))=nkoi(int(koi(i)))+1
         endif
 32   continue

      do 33 i=1,6
         np(i)=0
 33   continue

      do 34 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 34   continue

      write(6,500) np(1),"FA,Col,FPs,SN,P"

CCCCCCCCCC FPs,SN, Per, b

      do 35 i=1,nmax
         nkoi(i)=0
 35    continue

      do 36 i=1,n
         if((flag(i).lt.1).and.(nkoi1(int(koi(i))).eq.1).and.
     .    (sn(i).ge.10.0).and.(per(i).gt.1.6).and.
     .    (rdr(i)+b(i).lt.0.98).and.(ncol(i).eq.0)) then
           nkoi(int(koi(i)))=nkoi(int(koi(i)))+1
         endif
 36   continue

      do 37 i=1,6
         np(i)=0
 37   continue

      do 38 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 38   continue

      write(6,500) np(1),"FA,Col,FPs,SN,P,b"

CCCCCC EB counts

      do 43 i=1,nEB
         ebflag(i)=0
 43   continue

      do 41 i=1,nEB
         do 42 j=1,n
            if(kid(j).eq.ebkid(i))then
               ebflag(i)=1
            endif
 42      continue
 41   continue

      nEBc=0
      do 44 i=1,nEB
         if(ebflag(i).eq.0)then
            nEBc=nEBc+1
         endif
 44   continue

      write(6,500) nEBc, "EBs"

CCCCCC EBs, morph

      nEBc=0
      do 45 i=1,nEB
         if((ebflag(i).eq.0).and.(ebmorph(i).lt.0.5))then
            nEBc=nEBc+1
         endif
 45   continue

      write(6,500) nEBc, "EBs, morph"

CCCCCC EBs, morph, pdepth

      nEBc=0
      do 46 i=1,nEB
         if((ebflag(i).eq.0).and.(ebmorph(i).lt.0.5).and.
     .    (ebpdepth(i).lt.0.02))then
            nEBc=nEBc+1
         endif
 46   continue

      write(6,500) nEBc, "EBs, morph, pdep (2%)"

CCCCCC EBs, morph, pdepth, P

      nEBc=0
      do 47 i=1,nEB
         if((ebflag(i).eq.0).and.(ebmorph(i).lt.0.5).and.
     .    (ebpdepth(i).lt.0.02).and.(ebper(i).gt.1.6).and.
     .    (ebper(i).lt.1484.0))then
            nEBc=nEBc+1
         endif
 47   continue

      write(6,500) nEBc, "EBs, morph, pdep (2%), P"

CCCCCC EBs, pdepth

      nEBc=0
      do 48 i=1,nEB
         if((ebflag(i).eq.0).and.(ebpdepth(i).lt.0.02))then
            nEBc=nEBc+1
         endif
 48   continue

      write(6,500) nEBc, "EBs, pdep (2%)"

CCCCCC EBs, pdepth,P

      nEBc=0
      do 49 i=1,nEB
         if((ebflag(i).eq.0).and.
     .    (ebper(i).gt.1.6).and.(ebper(i).lt.1484.0))then
            nEBc=nEBc+1
         endif
 49   continue

      write(6,500) nEBc, "EBs, P"

CCCCCC EBs, pdepth,P

      nEBc=0
      do 50 i=1,nEB
         if((ebflag(i).eq.0).and.(ebpdepth(i).lt.0.02).and.
     .    (ebper(i).gt.1.6).and.(ebper(i).lt.1484.0))then
            nEBc=nEBc+1
         endif
 50   continue

      write(6,500) nEBc, "EBs, pdep (2%), P"

      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 902  write(0,*) "Error reading ",filename
      goto 999
 903  write(0,*) "Cannot open ",filenameEB
      goto 999
 904  write(0,*) "Error reading ",filename2
      goto 999
 905  write(0,*) "Cannot open ",filename2
      goto 999
 906  write(0,*) "Cannot open ",filenameFP
      goto 999
 907  write(0,*) "Cannot open ",filenameCC
      goto 999
 999  end
