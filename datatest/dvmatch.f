      program dvmatch
      implicit none
      integer nmax,i,nkoi,j
      parameter(nmax=10000)
      integer kid(nmax)
      real koi(nmax),dv(9)
      character*80 filename,dumc
      
      filename='koi2kic.txt' !KOI 2 KIC translations
      open(unit=10,file=filename,status='old',err=901)
      read(10,*) dumc !first two lines are comments
      read(10,*) dumc
      i=1
 10   read(10,*,end=11) koi(i),kid(i)
        i=i+1
      goto 10
 11   continue
      nkoi=i-1
      close(10)
      
      filename='dv_q2.dat'
      open(unit=10,file=filename,status='old',err=901)
 12   read(10,*,end=13) (dv(i),i=1,9)
        i=1
        do while(i.le.nkoi)
            if(int(dv(1)).eq.kid(i))then
                write(6,500) kid(i),(dv(j),j=2,7),dv(9),koi(i)
 500            format(I8,'| ',4(F13.7,'| '),2(F10.3,'| '),F6.1,'| ',
     .              F7.2)
                i=nkoi+1
            endif
            i=i+1
        enddo
      goto 12
 13   continue 
      
      goto 999
 901  write(6,*) "Cannot open ",filename
      goto 999
 999  end     
   
      
      