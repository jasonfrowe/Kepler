      program targetlist
      implicit none
      integer nmax,tq,i,nunit,kid,j,nkoi
      parameter(nmax=99999999,tq=8)
      integer nq(nmax),qc(8),koi(nmax),qkoi(tq),q1(nmax)
      real dumr
      character*80 filename

      nunit=10

      do 10 i=1,nmax
         nq(i)=0
         koi(i)=0
         q1(i)=0
 10   continue

      do 15 i=1,tq
         qc(i)=0
         qkoi(i)=0
 15   continue

      filename="multiplanetlist.txt"
      open(unit=nunit,file=filename,status='old',err=901)
 17   read(nunit,*,end=18) dumr,nkoi
         koi(nkoi)=1
      goto 17
 18   continue
      close(nunit)

      do 11 i=1,tq
         write(filename,500) 'q0',i,'list.txt'
         write(0,*) filename
 500     format(A2,I1,A8)

         open(unit=nunit,file=filename,status='old',err=901)
         j=1
 12      read(nunit,*,end=13,err=902) kid
            if(kid.gt.nmax)then
               write(0,*) "KID Error: ",kid
               goto 12
            endif
            nq(kid)=nq(kid)+1
            if(i.eq.1) q1(kid)=1
            j=j+1
         goto 12
 13      continue

         close(nunit)
 11   continue

      do 19 i=1,nmax
         if((nq(i).eq.1).and.(q1(i).eq.1))then
            write(6,*) i
         endif
 19   continue

c      do 14 i=1,nmax
c         if(nq(i).gt.0) then
c            qc(nq(i))=qc(nq(i))+1
c            write(6,*) i,nq(i) !write out KID and number of quarters obs
c         endif
c         if(koi(i).gt.0) qkoi(nq(i))=qkoi(nq(i))+1
c 14   continue
c
c      do 16 i=1,tq
c         write(6,501) i,qc(i),qkoi(i),100.0*real(qkoi(i))/real(qc(i))
c 501     format(I2,1X,I6,1X,I5,1X,F6.3)
c 16   continue

      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 902  write(0,*) "Error on line ",j
      write(0,*) "file: ",filename
      goto 999
 999  end


