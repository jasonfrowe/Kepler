      program klccomb
      implicit none
      integer iargc,nunit1,nunit2
      real*8 time,mag,cdppzero
      character*80 file1,file2
      
      cdppzero=54900.0d0+53.03815236d0
      
      if(iargc().lt.2) goto 901
      call getarg(1,file1)
      call getarg(2,file2)
      
      nunit1=10
      open(unit=nunit1,file=file1,status='old',err=902)
      nunit2=11
      open(unit=nunit2,file=file2,status='old',err=903)
      
 10   read(nunit1,*,end=11) time,mag
        write(6,*) time+cdppzero,mag
        goto 10
 11   continue
 
 12   read(nunit2,*,end=13) time,mag
        write(6,*) time,mag
        goto 12
 13   continue
        
      close(nunit1)
      close(nunit2)
      
      goto 999
 901  write(0,*) "Usage: klccomb <cdppdata> <q1data>"
      goto 999
 902  write(0,*) "Cannot open ",file1
      goto 999
 903  write(0,*) "Cannot open ",file2
      goto 999
 999  end