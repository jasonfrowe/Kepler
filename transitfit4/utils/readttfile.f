      subroutine readttfile(nunit,nmax,ntt,tobs,omc)
      implicit none
      integer nunit,ntt,i,nmax
      double precision tobs(nmax),omc(nmax)
      
      i=1
 10   read(nunit,*,end=11) tobs(i),omc(i)
        i=i+1
      goto 10
 11   continue
      ntt=i-1
      
      return
      end