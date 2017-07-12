      program limbtest
      implicit none
      integer nmax,i
      parameter(nmax=50)
      double precision p,cn(4),b(nmax),flux(nmax),flag,bmin,bmax,db
      
      p=1.2116222821E-02
      cn(1)=1.0860000000E+00
      cn(2)=-1.3660000000E+00
      cn(3)=1.8230000000E+00
      cn(4)=-6.7200000000E-01
      
      bmin=-1.2
      bmax= 1.2
      db=(bmax-bmin)/dble(nmax-1)
      do 10 i=1,nmax
        b(i)=abs(bmin+dble(i-1)*db)
 10   continue
      
      call occultnlv2(p,cn(1),cn(2),cn(3),cn(4),nmax,b,flux)
      
      end