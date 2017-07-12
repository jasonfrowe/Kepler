      program quarterobs
      implicit none
      integer nqmax,i,nkoi,nmax,j
      parameter(nqmax=8,nmax=5000)
      integer nq(nqmax),nunit,kid,qobs,koikid(nmax),kidmatch,nm(nqmax),
     .   nflag,npkoi(nmax)
      real koi(nmax),koimatch,data(nmax,26)
      character*80 multiprops,kidq1q8,koichar,dumc

      do 10 i=1,nqmax
         nq(i)=0  !initization of counters
         nm(i)=0
 10   continue

      multiprops="multiprops.dat" !file that contains multi-analysis
      kidq1q8="kidq1q8.txt" ! file that contains the # of Qs of observ
      koichar="koi_characteristics.20130509.dat" !KOI spreadsheet

      nunit=10  !this is the number number that gets used for I/O

C     open up multiprops file that has the multi-planet analysis
      open(unit=nunit,file=multiprops,status='old',err=902)
      read(nunit,*) dumc !first line is a header
      i=1
 14   read(nunit,*,end=15) (data(i,j),j=1,26)
         koi(i)=data(i,1)
c         write(0,*) i,koi(i)  !this line was for debugging
         i=i+1
      goto 14
 15   continue
      close(nunit)
      nkoi=i-1

C     Now we junk the FAs and make sure that there are still at least
C     2 good candidates
      do 20 i=1,nmax
         npkoi(i)=0  !counts number of planets for a given KOI, init=0
 20   continue

      do 21 i=1,nkoi
         if(data(i,21).ge.7.1)then
           npkoi(int(data(i,1)))=npkoi(int(data(i,1)))+1
         endif
 21   continue


C     Now we need to get the KOI-KID matches, so we use the KOI
C     spreadsheet
      open(unit=nunit,file=koichar,status='old',err=903)
      read(nunit,*) dumc !first line is a header
      read(nunit,*) dumc !so is the second line
 16   read(nunit,*,end=17) koimatch,kidmatch
         do 18 i=1,nkoi
c            write(6,*) i,npkoi(int(koi(i)))
            if((int(koi(i)).eq.int(koimatch)).and.
     .       (npkoi(int(koi(i))).gt.1))then
               koikid(i)=kidmatch
            endif
 18      continue
      goto 16
 17   continue
      close(nunit)

C     open up the file that has counts of number of targets with number
C     of quarters that have observations
      open(unit=nunit,file=kidq1q8,status='old',err=901)
C     each line has a KID and the number of quarters that have observs
 12   read(nunit,*,end=13) kid,qobs
         if((qobs.gt.0).and.(qobs.le.nqmax))then
            nq(qobs)=nq(qobs)+1  !count up number of KIDs
            nflag=0 !we only want to count a KOI once
            do 19 i=1,nmax
               if((kid.eq.koikid(i)).and.(nflag.eq.0))then
                  nm(qobs)=nm(qobs)+1
                  nflag=1 !flip bit to avoid recounting multi's
               endif
 19         continue
         endif
      goto 12
 13   continue
      close(nunit)

C     now we write out the number of targets that 1,2,3,..,8 Qs of obser
      do 11 i=1,nqmax
         write(6,500) i,nq(i),nm(i),100.0*real(nm(i))/real(nq(i))
 11   continue
 500  format(I1,1X,I6,1X,I6,1X,F4.2)

      goto 999
 901  write(0,*) "Cannot open ",kidq1q8
      goto 999
 902  write(0,*) "Cannot open ",multiprops
      goto 999
 903  write(0,*) "Cannot open ",koichar
      goto 999
 999  end
