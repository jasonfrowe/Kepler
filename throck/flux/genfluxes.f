CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine genfluxes(nT,Tmin,Tmax,Temps,Fluxes,npass,lam,pass)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     To speed up the program, the fluxes as a function of temperature 
C     Are precomputed.
      integer nT,i,npass
      double precision Tmin,Tmax,Temps(nT),Fluxes(nT),planck,T,dT,
     .  lam(npass),pass(npass)
      
      dT=(Tmax-Tmin)/real(nT-1)
      T=Tmin-dT
      do 10 i=1,nT
        T=T+dT
        Temps(i)=T
        Fluxes(i)=planck(T,npass,lam,pass)
 10   continue

      return
      end