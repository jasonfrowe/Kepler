      program readkeplerlc
C     Reads in Kepler long cadence data
      implicit none
      integer nT,nstar,nunit,nbzi,i,j
      parameter(nT=476,nstar=52496)
      integer kID(nstar)
      real time(nT),flux(nT,nstar)
      character*80 filename
      
      nunit=10
      filename="/home/rowe/Kepler/lcdata/nflux.dat"
      
C     Open the binary data file
      nbzi=nT*(nstar+1)
      open(unit=nunit,file=filename,status='old',form='unformatted',
     .  access='direct',recl=nbzi,err=901)
      read(nunit,rec=1) time,flux
      close(nunit)    
      
      filename="/home/rowe/Kepler/lcdata/kepId.txt"
      open(unit=nunit,file=filename,status='old',err=901)
      do 10 i=1,nstar
        read(nunit,*) kID(i)
 10   continue
      close(nunit)
      
      do 11 i=1,nstar
        if(kID(i).lt.10)then
            write(filename,507) "klc0000000",kID(i),".dat"
 507        format(A10,I1,A4)      
        elseif(kID(i).lt.100)then
            write(filename,506) "klc000000",kID(i),".dat"
 506        format(A9,I2,A4)      
        elseif(kID(i).lt.1000)then
            write(filename,505) "klc00000",kID(i),".dat"
 505        format(A8,I3,A4)      
        elseif(kID(i).lt.10000)then
            write(filename,504) "klc0000",kID(i),".dat"
 504        format(A7,I4,A4)      
        elseif(kID(i).lt.100000)then
            write(filename,503) "klc000",kID(i),".dat"
 503        format(A6,I5,A4)      
        elseif(kID(i).lt.1000000)then
            write(filename,502) "klc00",kID(i),".dat"
 502        format(A5,I6,A4)
        elseif(kID(i).lt.10000000)then
            write(filename,501) "klc0",kID(i),".dat"
 501        format(A4,I7,A4)
        else
            write(filename,500) "klc",kID(i),".dat"
 500        format(A3,I8,A4)
        endif
        
c        write(6,508) filename
c 508    format(A15)
        
        open(unit=nunit,file=filename)
        do 12 j=1,nT
            write(nunit,*) time(j),flux(j,i)
 12     continue       
        close(nunit)
 
 11   continue
      
      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 999  end