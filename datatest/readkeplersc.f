      program readkeplersc
C     Reads in Kepler long cadence data
      implicit none
      integer nT,nstar,nunit,nbzi,i,j
      parameter(nT=14280,nstar=310)
      integer kID(nstar)
      real flux(nT,nstar)
      real*8 time(nT)
      character*80 filename
      
      nunit=10
      filename="/home/rowe/Kepler/scdata/scFluxTimeSeries.dat"
      
C     Open the binary data file
      nbzi=nT*(nstar)
      open(unit=nunit,file=filename,status='old',form='unformatted',
     .  access='direct',recl=nbzi,err=901)
      read(nunit,rec=1) flux
      close(nunit)  
      
      filename="/home/rowe/Kepler/scdata/scCadenceTimesr8.dat"
      nbzi=2*nT
      open(unit=nunit,file=filename,status='old',form='unformatted',
     .  access='direct',recl=nbzi,err=901)
      read(nunit,rec=1) time
      close(nunit)  
      
c      do 10 i=1,nT
c        write(6,*) time(i),flux(i,1)
c 10   continue
      
      filename="/home/rowe/Kepler/scdata/cdpp.idse"
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