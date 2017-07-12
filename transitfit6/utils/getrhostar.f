      subroutine getrhostar(nunit,np,nmax,mstar,age,z,rstar,rhostar,
     .  temp,lum)
      implicit none
      integer np,nmax,nunit,i
      double precision mstar(nmax),age(nmax),z(nmax),rstar(nmax),
     .  rhostar(nmax),temp(nmax),lum(nmax)
     
      i=1
 10   read(nunit,*,end=11) mstar(i),age(i),z(i),rstar(i),rhostar(i),
     .  temp(i),lum(i)
        lum(i)=10.0**lum(i)
        rhostar(i)=rhostar(i)/1.0d3
        i=i+1
      goto 10
 11   continue
     
      np=i-1
      write(0,*) "Stellar parameters read:",np
     
      return
      end