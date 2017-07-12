CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readKeplerbandpass(nKepler,Keplerlam,Keplerpass)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nKepler,i,tempi
      double precision Keplerlam(nKepler),Keplerpass(nKepler)
      character*80 filename
      
      filename='kepler_trans.dat'
      
      open(unit=11,file=filename,status='old',err=901)
      
      do 10 i=1,nKepler
        read(11,*) tempi,Keplerpass(i)
        Keplerlam(i)=dble(tempi)/10.0 !convert from A to nm
 10   continue
 
      close(11)
      goto 999
 901  write(0,*) "Cannot open: ",filename
      stop
      goto 999
 999  return
      end
