      program multistatsV2
      implicit none
      integer nmax,n,i,j,nex,k,ii,jj,npmax,i3,nkex,ndyn,nunit,nsys,
     .   nblarge
      parameter(nmax=5000,nex=17,npmax=10,nkex=3)
      integer nkoi(nmax),np(6),newkoi,sump,nkoifp(nmax),nkoi1(nmax),
     .   nkoikeep(nmax),nK,planetN(npmax),nsort(npmax),id1,id2,
     .   kepid1(nmax),kepid2(nmax),istart,nKep,excl(nmax)
      real data(nmax,26),pexcl(nex),planetPer(npmax),koi,kexcl(nkex),
     .   dynkoi(nmax),dumr
      character nKc,pc(npmax),dynchar
      character*80 filename,dumc,kepnames,dynfile
      data pexcl/375.01,422.01,435.02,490.02,771.01,1032.01,1096.01,
     .   1192.01,1206.01,1274.01,1356.01,1421.01,1463.01,1174.01,
     .   1208.01,1008.01,99.01/
      data kexcl/1750.01,1750.02,1955.03/

      nunit=10
      dynfile="dynamictest.20130905.dat"
      open(unit=nunit,file=dynfile,status='old',err=904)
      read(nunit,*) dumc
      i=1
 126  read(nunit,*,end=127) dynkoi(i),dumr,dumr,dumr,dynchar
         if(dynchar.eq.'C')then
c            write(0,*) dynkoi(i)
            i=i+1
         endif
      goto 126
 127  continue
      ndyn=i-1
      close(nunit)

      pc(1)="b"
      pc(2)="c"
      pc(3)="d"
      pc(4)="e"
      pc(5)="f"
      pc(6)="g"
      pc(7)="h"
      pc(8)="i"
      pc(9)="j"
      pc(10)="k"

      do 105 i=1,nmax
         kepid1(i)=0
         kepid2(i)=0
         excl(i)=0
 105  continue

      filename="multiprops.dat"
      open(unit=10,file=filename,status='old',err=901)

      read(10,*) dumc !first line is a comment
      n=1
 10   read(10,*,end=11,err=902) (data(n,i),i=1,26) !read in data
c         write(0,*) n, data(n,1)

         do 70 j=1,nex
            if(pexcl(j).eq.data(n,1))then
c               write(0,*) pexcl(j),data(n,1)
               data(n,4)=-1.0
             endif
 70      continue

         do 125 j=1,nkex
            if(kexcl(j).eq.data(n,1))then
c               write(0,*) kexcl(j),data(n,1)
               excl(n)=1
            endif
 125     continue

         do 128 j=1,ndyn
            if(dynkoi(j).eq.data(n,1))then
c               write(0,*) dynkoi(j),data(n,1)
               excl(n)=1
            endif
 128     continue


         n=n+1
      goto 10
 11   continue
      n=n-1
      close(10)

      kepnames="keplernames.txt"
      open(unit=10,file=kepnames,status='old',err=902)
 102  read(10,*,end=103) koi,id1,id2
         do 104 j=1,n
            if(data(j,1).eq.koi)then
               kepid1(j)=id1
               kepid2(j)=id2
c               write(6,*) data(j,1),kepid1(j),kepid2(j)
            endif
 104     continue
      goto 102
 103  continue
      close(10)


CCCCCCCCCCCCCCCCC
C     no cuts..

      do 4 i=1,nmax
         nkoi(i)=0
 4    continue

      do 20 i=1,n
           nkoi(int(data(i,1)))=nkoi(int(data(i,1)))+1
 20   continue

      do 21 i=1,6
         np(i)=0
 21   continue

      do 22 i=1,nmax
         nkoi1(i)=nkoi(i)
c         if(nkoi1(i).eq.6) write(0,*) i,nkoi1(i)
         if(nkoi(i).eq.1) write(0,*) i,nkoi(i)
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 22   continue

      newkoi=0
      do 23 i=1,n
         if(data(i,26).eq.0)then
            newkoi=newkoi+1
         endif
 23   continue

      sump=0
      do 49 i=2,6
         sump=sump+i*np(i)
 49   continue

      write(6,500) (np(i),i=1,6),sump,newkoi,"No Cuts"

C     only keep 6 planet systems
c      j=0
c      do 82 i=1,nmax
c         if(nkoi1(i).eq.2)then
c            j=j+1
c            nkoikeep(j)=i
c            write(0,*) nkoikeep(j)
c         endif
c 82   continue
c      ii=0
c      do 83 i=1,n
c         do 84 k=1,j
c            if(int(data(i,1)).eq.nkoikeep(k))then
c               ii=ii+1
c               do 85 jj=1,26
c                  data(ii,jj)=data(i,jj)
c 85            continue
c            endif
c 84      continue
c 83   continue
c      n=ii


CCCCCCCCCCCCCCC
C     junk cuts

      do 69 i=1,nmax
         nkoi(i)=0
 69   continue

      do 75 i=1,n
         if(data(i,21).ge.7.1)then
           nkoi(int(data(i,1)))=nkoi(int(data(i,1)))+1
         endif
 75   continue

      do 71 i=1,6
         np(i)=0
 71   continue

      do 72 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 72   continue

      nblarge=0
      newkoi=0
      do 73 i=1,n
c         if((nkoi(int(data(i,1))).gt.0).and.(data(i,24).le.1))then
c            write(0,*) data(i,1),nkoi(int(data(i,1))),int(data(i,24)),
c     .         data(i,21)
c         endif

         if((nkoi(int(data(i,1))).gt.1).and.
     .    (data(i,26).eq.0).and.(data(i,21).ge.7.1))then
            newkoi=newkoi+1
            if(data(i,6)+data(i,7)+data(i,8).gt.1.00) nblarge=nblarge+1
         endif
 73   continue
c      write(6,*) "nblarge:",nblarge

      sump=0
      do 74 i=2,6
         sump=sump+i*np(i)
 74   continue

      write(6,500) (np(i),i=1,6),sump,newkoi,"FA"

CCCCCCCCCCCCCCC
C     cut P/T0 collisions

      do 89 i=1,nmax
         nkoi(i)=0
 89   continue

      do 90 i=1,n
         if((data(i,21).ge.7.1).and.((data(i,24).gt.1).or.
     .    (data(i,24).eq.0)))then
           nkoi(int(data(i,1)))=nkoi(int(data(i,1)))+1
         endif
 90   continue

      do 91 i=1,6
         np(i)=0
 91   continue

      do 92 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 92   continue

      newkoi=0
      do 93 i=1,n
         if((nkoi(int(data(i,1))).gt.1).and.
     .    (data(i,26).eq.0).and.(data(i,21).ge.7.1).and.
     .    ((data(i,24).gt.1).or.(data(i,24).eq.0)))then
            newkoi=newkoi+1
         endif
 93   continue

      sump=0
      do 94 i=2,6
         sump=sump+i*np(i)
 94   continue

      write(6,500) (np(i),i=1,6),sump,newkoi,"FA,Col"

CCCCCCCCCCCCCCC
C     cut P/T0 collisions and Period

      do 107 i=1,nmax
         nkoi(i)=0
 107   continue

      do 108 i=1,n
         if((data(i,21).ge.7.1).and.((data(i,24).gt.1).or.
     .    (data(i,24).eq.0)).and.(data(i,4).gt.1.6))then
           nkoi(int(data(i,1)))=nkoi(int(data(i,1)))+1
         endif
 108   continue

      do 109 i=1,6
         np(i)=0
 109  continue

      do 110 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 110  continue

      newkoi=0
      do 111 i=1,n
         if((nkoi(int(data(i,1))).gt.1).and.
     .    (data(i,26).eq.0).and.(data(i,21).ge.7.1).and.
     .    ((data(i,24).gt.1).or.(data(i,24).eq.0)).and.
     .    (data(i,4).gt.1.6))then
            newkoi=newkoi+1
         endif
 111  continue

      sump=0
      do 112 i=2,6
         sump=sump+i*np(i)
 112  continue

      write(6,500) (np(i),i=1,6),sump,newkoi,"FA,Col,P"

CCCCCCCCCCCCCCCCCCCCCCCCC
C     FP cuts

      do 24 i=1,nmax
         nkoi(i)=0
 24   continue

      do 25 i=1,n
         if((data(i,24).gt.1).and.(data(i,21).ge.7.1))then
           nkoi(int(data(i,1)))=nkoi(int(data(i,1)))+1
c         else
c            write(0,*) data(i,1), data(i,21), data(i,24)
         endif
 25   continue

      do 26 i=1,6
         np(i)=0
 26   continue

      do 27 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 27   continue

      nblarge=0
      newkoi=0
      do 28 i=1,n
         if((nkoi(int(data(i,1))).gt.1).and.(data(i,24).gt.1).and.
     .    (data(i,26).eq.0).and.(data(i,21).ge.7.1))then
            newkoi=newkoi+1
            if(data(i,6)+data(i,7)+data(i,8).gt.1.00) nblarge=nblarge+1
         endif
 28   continue
c      write(6,*) "nblarge: ",nblarge

      sump=0
      do 50 i=2,6
         sump=sump+i*np(i)
 50   continue

      write(6,500) (np(i),i=1,6),sump,newkoi,"FA,Col,FPs"

CCCCCCCCCCCCCCCCCCCCCCCCC
C     FP cuts, P

      do 119 i=1,nmax
         nkoi(i)=0
 119  continue

      do 120 i=1,n
         if((data(i,24).gt.1).and.(data(i,21).ge.7.1).and.
     .    (data(i,4).gt.1.6))then
           nkoi(int(data(i,1)))=nkoi(int(data(i,1)))+1
c           else
c           write(0,*) data(i,1)
         endif
 120  continue

      do 121 i=1,6
         np(i)=0
 121   continue

      do 122 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 122  continue

      newkoi=0
      do 123 i=1,n
         if((nkoi(int(data(i,1))).gt.1).and.(data(i,24).gt.1).and.
     .    (data(i,26).eq.0).and.(data(i,21).ge.7.1).and.
     .    (data(i,4).gt.1.6))then
            newkoi=newkoi+1
         endif
 123   continue

      sump=0
      do 124 i=2,6
         sump=sump+i*np(i)
 124   continue

      write(6,500) (np(i),i=1,6),sump,newkoi,"FA,Col,FPs,P"

CCCCCCCCCCCCCCCCCCCCCCCCCCC
C     FA,Col,FP,SN cuts

      do 34 i=1,nmax
         nkoi(i)=0
 34   continue

      do 35 i=1,n
         if((data(i,24).gt.1).and.(data(i,21).ge.10.0))then
           nkoi(int(data(i,1)))=nkoi(int(data(i,1)))+1
         endif
 35   continue

      do 36 i=1,6
         np(i)=0
 36   continue

      do 37 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 37   continue

      newkoi=0
      do 38 i=1,n
         if((nkoi(int(data(i,1))).gt.1).and.(data(i,24).gt.1).and.
     .    (data(i,26).eq.0).and.(data(i,21).ge.10.0))then
            newkoi=newkoi+1
         endif
 38   continue

      sump=0
      do 51 i=2,6
         sump=sump+i*np(i)
 51   continue

      write(6,500) (np(i),i=1,6),sump,newkoi,"FA,Col,FPs,SN"

CCCCCCCCCCCCCCCCCCCCCCCCCCC
C     FA,Col,FP,SN,b cuts

      do 113 i=1,nmax
         nkoi(i)=0
 113   continue

      do 114 i=1,n
         if((data(i,24).gt.1).and.(data(i,21).ge.10.0).and.
     .    (data(i,6)+data(i,7)+data(i,8).lt.1.00))then
           nkoi(int(data(i,1)))=nkoi(int(data(i,1)))+1
         endif
 114  continue

      do 115 i=1,6
         np(i)=0
 115  continue

      do 116 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 116  continue

      newkoi=0
      do 117 i=1,n
         if((nkoi(int(data(i,1))).gt.1).and.(data(i,24).gt.1).and.
     .    (data(i,26).eq.0).and.(data(i,21).ge.10.0).and.
     .    (data(i,6)+data(i,7)+data(i,8).lt.1.00))then
            newkoi=newkoi+1
         endif
 117  continue

      sump=0
      do 118 i=2,6
         sump=sump+i*np(i)
 118  continue

      write(6,500) (np(i),i=1,6),sump,newkoi,"FA,Col,FPs,SN,b"

C     FP,SN,P cuts

      do 39 i=1,nmax
         nkoi(i)=0
 39   continue

      do 40 i=1,n
         if((data(i,24).gt.1).and.(data(i,21).ge.10.0).and.
     .    (data(i,4).gt.1.6))then
           nkoi(int(data(i,1)))=nkoi(int(data(i,1)))+1
         endif
 40   continue

      do 41 i=1,6
         np(i)=0
 41   continue

      do 42 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 42   continue

      newkoi=0
      do 43 i=1,n
         if((nkoi(int(data(i,1))).gt.1).and.(data(i,24).gt.1).and.
     .    (data(i,26).eq.0).and.(data(i,21).ge.10.0).and.
     .    (data(i,4).gt.1.6))then
            newkoi=newkoi+1
         endif
 43   continue

      sump=0
      do 52 i=2,6
         sump=sump+i*np(i)
 52   continue

      write(6,500) (np(i),i=1,6),sump,newkoi,"FA,Col,FPs,SN,P"

C     FP,SN,P,b cuts

      do 44 i=1,nmax
         nkoi(i)=0
 44   continue

      do 45 i=1,n
         if((data(i,24).gt.1).and.(data(i,21).ge.10.0).and.
     .    (data(i,6)+data(i,7)+data(i,8).lt.1.00).and.
     .    (data(i,4).gt.1.6))then
           nkoi(int(data(i,1)))=nkoi(int(data(i,1)))+1
         endif
 45   continue

      do 46 i=1,6
         np(i)=0
 46   continue

      do 47 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 47   continue

      newkoi=0
      do 48 i=1,n
         if((nkoi(int(data(i,1))).gt.1).and.(data(i,24).gt.1).and.
     .    (data(i,26).eq.0).and.(data(i,21).ge.10.0).and.
     .    (data(i,6)+data(i,7)+data(i,8).lt.1.00).and.
     .    (data(i,4).gt.1.6))then
            newkoi=newkoi+1
         endif
 48   continue

      sump=0
      do 53 i=2,6
         sump=sump+i*np(i)
 53   continue

      write(6,500) (np(i),i=1,6),sump,newkoi,"FA,Col,FPs,SN,P,b"

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Cutting FAs and centroid non-passes
      do 76 i=1,nmax
         nkoi(i)=0
         nkoifp(i)=0
 76   continue

      do 77 i=1,n
         if((data(i,23).gt.3).and.(data(i,21).gt.7.1))then
            nkoi(int(data(i,1)))=nkoi(int(data(i,1)))+1
         else
            nkoifp(int(data(i,1)))=nkoifp(int(data(i,1)))+1
         endif
 77   continue

      do 78 i=1,6
         np(i)=0
 78   continue

      do 79 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 79   continue

      newkoi=0
      do 80 i=1,n
         if((nkoi(int(data(i,1))).gt.1).and.(data(i,23).gt.3).and.
     .    (data(i,21).gt.7.1).and.(data(i,26).eq.0))then
            newkoi=newkoi+1
         endif
 80   continue

      sump=0
      do 81 i=2,6
         sump=sump+i*np(i)
 81   continue

      write(6,500) (np(i),i=1,6),sump,newkoi,"FA,cent"

CCCCCCCCCCCCCCCCC
C     Cutting FAs, Collisions and centroid non-passes
      do 95 i=1,nmax
         nkoi(i)=0
         nkoifp(i)=0
 95   continue

      do 96 i=1,n
         if((data(i,23).gt.3).and.(data(i,21).gt.7.1).and.
     .    ((data(i,24).gt.1).or.(data(i,24).eq.0)))then
            nkoi(int(data(i,1)))=nkoi(int(data(i,1)))+1
         else
            nkoifp(int(data(i,1)))=nkoifp(int(data(i,1)))+1
         endif
 96   continue

      do 97 i=1,6
         np(i)=0
 97   continue

      do 98 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 98   continue

      newkoi=0
      do 99 i=1,n
         if((nkoi(int(data(i,1))).gt.1).and.(data(i,23).gt.3).and.
     .    (data(i,21).gt.7.1).and.(data(i,26).eq.0).and.
     .    ((data(i,24).gt.1).or.(data(i,24).eq.0)))then
            newkoi=newkoi+1
         endif
 99   continue

      sump=0
      do 101 i=2,6
         sump=sump+i*np(i)
 101  continue

      write(6,500) (np(i),i=1,6),sump,newkoi,"FA,Col,cent"

CCCCCCCCCCCCCCCCC
C     Cutting FA, collisons FPs and centroids non-passes

      do 5 i=1,nmax
         nkoi(i)=0
         nkoifp(i)=0
 5    continue

      do 12 i=1,n
         if((data(i,24).gt.1).and.(data(i,23).gt.3).and.
     .    (data(i,21).gt.7.1))then
            nkoi(int(data(i,1)))=nkoi(int(data(i,1)))+1
         else
            nkoifp(int(data(i,1)))=nkoifp(int(data(i,1)))+1
         endif
 12   continue

c      do 86 i=1,nmax
c         if((nkoi(i)+nkoifp(i).eq.3))then
c            write(0,*) i,nkoi(i),nkoifp(i)
c         endif
c 86   continue

      do 14 i=1,6
         np(i)=0
 14   continue

      do 13 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 13   continue

      newkoi=0
      do 15 i=1,n
         if((nkoi(int(data(i,1))).gt.1).and.(data(i,24).gt.1).and.
     .    (data(i,23).gt.3).and.(data(i,26).eq.0).and.
     .    (data(i,21).gt.7.1))then
            newkoi=newkoi+1
         endif
 15   continue

      sump=0
      do 54 i=2,6
         sump=sump+i*np(i)
 54   continue

      write(6,500) (np(i),i=1,6),sump,newkoi,"FA,Col,FP,cent"
 500  format(8(I4,1X),A25)

CCCCCCCCCCCC
C     Cutting FPs and centroids non-passes, SN

      do 57 i=1,nmax
         nkoi(i)=0
 57   continue

      do 58 i=1,n
         if((data(i,24).gt.1).and.(data(i,23).gt.3).and.
     .    (data(i,21).ge.10.0))then
            nkoi(int(data(i,1)))=nkoi(int(data(i,1)))+1
         endif
 58   continue

      do 59 i=1,6
         np(i)=0
 59   continue

      do 60 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 60   continue

      newkoi=0
      do 61 i=1,n
         if((nkoi(int(data(i,1))).gt.1).and.(data(i,24).gt.1).and.
     .    (data(i,23).gt.3).and.(data(i,26).eq.0).and.
     .    (data(i,21).ge.10.0))then
            newkoi=newkoi+1
         endif
 61   continue

      sump=0
      do 62 i=2,6
         sump=sump+i*np(i)
 62   continue

      write(6,500) (np(i),i=1,6),sump,newkoi,"FA,Col,FP,cent,SN"

CCCCCCCCCCCCCCCCC
C     Cutting FPs and centroids non-passes, SN, P

      do 63 i=1,nmax
         nkoi(i)=0
 63   continue

      do 64 i=1,n
         if((data(i,24).gt.1).and.(data(i,23).gt.3).and.
     .    (data(i,21).ge.10.0).and.(data(i,4).gt.1.6))then
            nkoi(int(data(i,1)))=nkoi(int(data(i,1)))+1
         endif
 64   continue

      do 65 i=1,6
         np(i)=0
 65   continue

      do 66 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 66   continue

      newkoi=0
      do 67 i=1,n
         if((nkoi(int(data(i,1))).gt.1).and.(data(i,24).gt.1).and.
     .    (data(i,23).gt.3).and.(data(i,26).eq.0).and.
     .    (data(i,21).ge.10.0).and.(data(i,4).gt.1.6))then
            newkoi=newkoi+1
         endif
 67   continue

      sump=0
      do 68 i=2,6
         sump=sump+i*np(i)
 68   continue

      write(6,500) (np(i),i=1,6),sump,newkoi,"FA,Col,FP,cent,SN,P"

CCCCCCCCCCCCCCCCC
C     Cutting FPs and centroids non-passes, SN, P, b

      do 6 i=1,nmax
         nkoi(i)=0
 6    continue

      do 16 i=1,n
         if((data(i,24).gt.1).and.(data(i,23).gt.3).and.
     .      (data(i,4).gt.1.6).and.
     .      (data(i,6)+data(i,7)+data(i,8).lt.1.00).and.
     .      (data(i,21).ge.10.0))then
            nkoi(int(data(i,1)))=nkoi(int(data(i,1)))+1
         endif
 16   continue

      do 17 i=1,6
         np(i)=0
 17   continue

      do 18 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 18   continue

      newkoi=0
      do 19 i=1,n
         if((nkoi(int(data(i,1))).gt.1).and.(data(i,24).gt.1).and.
     .    (data(i,23).gt.3).and.(data(i,4).gt.1.6).and.
     .    (data(i,6)+data(i,7)+data(i,8).lt.1.00).and.
     .    (data(i,26).eq.0).and.
     .    (data(i,21).ge.10.0))then
            newkoi=newkoi+1
         endif
 19   continue

      sump=0
      do 55 i=2,6
         sump=sump+i*np(i)
 55   continue

      write(6,500) (np(i),i=1,6),sump,newkoi,"FA,Col,FP,cent,SN,P,b"

CCCCCCC Cut dynamical short period fails

      do 129 i=1,nmax
         nkoi(i)=0
 129  continue

      do 134 i=1,n
         if((data(i,24).gt.1).and.(data(i,23).gt.3).and.
     .      (data(i,4).gt.1.6).and.
     .      (data(i,6)+data(i,7)+data(i,8).lt.1.00).and.
     .      (data(i,21).ge.10.0).and.(excl(i).eq.0))then
            nkoi(int(data(i,1)))=nkoi(int(data(i,1)))+1
         endif
 134   continue

      do 130 i=1,6
         np(i)=0
 130   continue

      do 131 i=1,nmax
         if(nkoi(i).gt.0)then
            np(nkoi(i))=np(nkoi(i))+1
         endif
 131   continue

      newkoi=0
      do 132 i=1,n
         if((nkoi(int(data(i,1))).gt.1).and.(data(i,24).gt.1).and.
     .    (data(i,23).gt.3).and.(data(i,4).gt.1.6).and.
     .    (data(i,6)+data(i,7)+data(i,8).lt.1.00).and.
     .    (data(i,26).eq.0).and.
     .    (data(i,21).ge.10.0).and.(excl(i).eq.0))then
c            write(6,*) data(i,1),data(i,6),data(i,7),data(i,8)
            newkoi=newkoi+1
         endif
 132   continue

      sump=0
      do 133 i=2,6
         sump=sump+i*np(i)
 133   continue

      write(6,500) (np(i),i=1,6),sump,newkoi,"Dynamical and SP Cuts"


      goto 900 !skip KOI number generation

C     Make list of KOIs that get Kepler numbers
      nsys=0 !count number of systems
      newkoi=0
      nKep=99
      do 86 i=1,nmax
         if(nkoi(i).gt.1)then
c            write(6,*) i,nkoi(i)

            istart=1
            nK=0
            do 106 j=1,n
               if((i.eq.int(data(j,1))).and.(kepid1(j).gt.0))then
                  nK=kepid1(j)
                  istart=max(istart,kepid2(j)+1)
c                  write(6,*) "nK: ",nK,istart
               endif
 106        continue

            if(nK.eq.0) then
               nKep=nKep+1
               nK=nKep
            endif

c            write(0,*) i,nkoi(i),nK
            k=0 !counter for adding in new planets
            do 87 j=1,n
               if((i.eq.int(data(j,1))).and.(data(j,24).gt.1).and.
     .          (data(j,23).gt.3).and.(data(j,4).gt.1.6).and.
     .          (data(j,6)+data(j,7)+data(j,8).lt.1.00).and.
     .          (data(j,26).eq.0).and.
     .          (data(j,21).ge.10.0).and.(kepid1(j).eq.0).and.
     .          (excl(j).eq.0))then
                  k=k+1
                  PlanetN(k)=j
                  PlanetPer(k)=data(j,4)
c                  write(0,*) data(j,1),data(j,23),data(j,24)
                endif
 87         continue
            call rqsort(k,PlanetPer,nsort)
            do 88 ii=1,k
               i3=ii+istart-1
               if(i.lt.10)then
                  jj=PlanetN(nsort(ii))
                  if(nK.lt.100)then
                     write(6,509) data(jj,1),",Kepler-",nK,pc(i3)
                  elseif(nK.lt.1000)then
                     write(6,501) data(jj,1),",Kepler-",nK,pc(i3)
                  else
                     write(6,502) data(jj,1),",Kepler-",nK,pc(i3)
                  endif
               elseif(i.lt.100)then
                  jj=PlanetN(nsort(ii))
                  if(nK.lt.100)then
                     write(6,510) data(jj,1),",Kepler-",nK,pc(i3)
                  elseif(nK.lt.1000)then
                     write(6,503) data(jj,1),",Kepler-",nK,pc(i3)
                  else
                     write(6,504) data(jj,1),",Kepler-",nK,pc(i3)
                  endif
               elseif(i.lt.1000)then
                  jj=PlanetN(nsort(ii))
                  if(nK.lt.100)then
                     write(6,511) data(jj,1),",Kepler-",nK,pc(i3)
                  elseif(nK.lt.1000)then
                     write(6,505) data(jj,1),",Kepler-",nK,pc(i3)
                  else
                     write(6,506) data(jj,1),",Kepler-",nK,pc(i3)
                  endif
               else
                  jj=PlanetN(nsort(ii))
                  if(nK.lt.100)then
                     write(6,512) data(jj,1),",Kepler-",nK,pc(i3)
                  elseif(nK.lt.1000)then
                     write(6,507) data(jj,1),",Kepler-",nK,pc(i3)
                  else
                     write(6,508) data(jj,1),",Kepler-",nK,pc(i3)
                  endif
               endif
 88         continue
c         read(5,*)
         endif
 86   continue

 509  format(F4.2,A8,I2,1X,A1)
 501  format(F4.2,A8,I3,1X,A1)
 502  format(F4.2,A8,I4,1X,A1)

 510  format(F5.2,A8,I2,1X,A1)
 503  format(F5.2,A8,I3,1X,A1)
 504  format(F5.2,A8,I4,1X,A1)

 511  format(F6.2,A8,I2,1X,A1)
 505  format(F6.2,A8,I3,1X,A1)
 506  format(F6.2,A8,I4,1X,A1)

 512  format(F7.2,A8,I2,1X,A1)
 507  format(F7.2,A8,I3,1X,A1)
 508  format(F7.2,A8,I4,1X,A1)

c         if((nkoi(int(data(i,1))).gt.1).and.(data(i,22).eq.0).and.
c     .    (data(i,23).gt.1).and.(data(i,4).gt.1.6).and.
c     .    (data(i,6)+data(i,8).lt.0.98).and.(data(i,26).eq.0).and.
c     .    (data(i,21).ge.10.0))then
c            write(6,501) data(i,1),"Kepler-",nK,nKc
c         endif
c 86   continue
c 501  format(F7.2,1X,A7,)

 900  continue

      goto 999
 901  write(0,*) "Cannot open: ",filename
      goto 999
 902  write(0,*) "Error on line ",n
      goto 999
 903  write(0,*) "Cannot open: ",kepnames
      goto 999
 904  write(0,*) "Cannot open: ",dynfile
      goto 999
 999  end


