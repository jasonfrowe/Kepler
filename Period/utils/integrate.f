CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function integrate(npt,time,mag,period,k,sc)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This function evaluates the fourier co-efficents though integral
C     approximations.  These parameters are then used a first guess for 
C     least square minimization
C
C     npt - number of data points
C     time - times of observations
C     mag - the observations
C     period - fixed period 
C     k - order of fourier series
C     sc - sine (=1) or cosine (=2) term
      implicit none
      integer nmax
      parameter (nmax=2000000)
      integer npt,iord(nmax),i,k,sc
      double precision time(npt),mag(npt),period,fint,phase(nmax),pi

      pi=3.141592654

      call phasept(npt,time,mag,phase,period)
      call rqsort(npt,phase,iord)
      
      do 5 i=1,npt
         phase(i)=phase(i)-0.5
 5    continue
      
      if(k.lt.1) then
         fint=(phase(iord(1))+0.5)*mag(iord(1))
         do 10 i=2,npt
            fint=fint+(phase(iord(i))-phase(iord(i-1)))*mag(iord(i))
 10      continue
      else
         if(sc.eq.2) then
            fint=(phase(iord(1))+0.5)*mag(iord(1))*cos(2.0*dble(k)*pi*
     .           phase(iord(1)))
            do 20 i=2,npt
               fint=fint+(phase(iord(i))-phase(iord(i-1)))*mag(iord(i))*
     .              cos(2.0*dble(k)*pi*phase(iord(i)))
 20         continue
            fint=2.0*fint
         elseif(sc.eq.1) then
            fint=(phase(iord(1))+0.5)*phase(iord(1))*sin(2.0*dble(k)*pi*
     .           mag(iord(1)))
            do 30 i=2,npt
               fint=fint+(phase(iord(i))-phase(iord(i-1)))*mag(iord(i))*
     .              sin(2.0*dble(k)*pi*phase(iord(i)))
 30         continue
            fint=2.0*fint
         else
            fint=-1.0
         endif
      endif

      integrate=fint

      return
      end
