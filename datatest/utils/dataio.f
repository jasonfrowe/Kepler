CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readkeplc(nunit,nmax,npt,dtime,flux,ferr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nmax,npt,i
      double precision dtime(nmax),flux(nmax),ferr(nmax)

      i=1
      
  9   continue
 10   read(nunit,*,err=9,end=20) dtime(i),flux(i),ferr(i)
c        mag(i)=-2.5*log10(mag(i)+1.0d0)
c        ferr(i)=0.00005
        i=i+1
      goto 10
 20   continue
        
      npt=i-1
      write(0,*)   "-------------------------"
      write(0,500) "Observations read: ",npt
      write(0,*)   "-------------------------"
 500  format(1X,A19,I8)
 
      return
      end