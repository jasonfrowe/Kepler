CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function obsflux(nlat,nlon,nR,nT,vangle,theta,
     .  lambda,T,Tgridmin,Tgridmax,pi,pid2,tpi,Temps,Fluxes,yT,sarea,Rp,
     .  aveT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nlat,nlon,nR,i,j,k,nT
      double precision vangle(2),theta(nlat),lambda(nlon),a1,a2,
     .  T(nlat,nlon,nR),Tmin,Tmax,Tgridmin,Tgridmax,dflux,pi,Pid2,alpha,
     .  calpha,Temps(nT),Fluxes(nT),yT(nT),flux,sarea(nlat),Rp,aveT,tpi
      
      aveT=0.0d0
      flux=0.0d0
      k=1 !surface element
      Tmax=T(1,1,1)
      Tmin=T(1,1,1)
      do 10 i=1,nlat
        a1=theta(i)
        do 11 j=1,nlon
            Tmax=max(T(i,j,k),Tmax) !maximum temperature
            Tmin=min(T(i,j,k),Tmin)
            if(Tmin.lt.Tgridmin) goto 902 !interpolation error
            if(Tmax.gt.Tgridmax) goto 903
            a2=lambda(j)
C     Calculate angle of surface relative to observer norm
            alpha=acos(sin(a1-Pi/2.0e0)*sin(vangle(1)-Pi/2.0e0)
     .          +cos(a1-Pi/2.0e0)*cos(vangle(1)-Pi/2.0e0)
     .          *cos(a2-vangle(2)))
            calpha=cos(abs(alpha))
                
C     This line is for direct calculation of fluxes (aka.. ssslllooowww)
c                dflux=planck(T(i,j),nMOST,MOSTlam,MOSTpass)*calpha
C     Next two lines use interpolation to get fluxes
c            write(6,*) Temps
            call splint(Temps,Fluxes,yT,nT,T(i,j,k),dflux)
            dflux=dflux*calpha
C     making sure surface element is visible            
            if((alpha.ge.pid2).or.(calpha.lt.0.0)) dflux=0.0
            flux=flux+dflux*sarea(i)
C     Calculate average surface temperature
            if((alpha.lt.pid2).and.(calpha.ge.0.0)) 
     .          aveT=aveT+T(i,j,k)*sarea(i) ! [K]
 11     continue
 10   continue
C     -4pi stuff (actually 2 Pi) is taken care of my the surface element 
C       areas
C     -convert from ergs s^-1 cm^-2 to W m^-2
      obsflux=flux*Rp*Rp*1.0d-3  ![W m^-2]
C     Average surface Temperature
      aveT=aveT/tpi ![K]


      goto 999
 902  write(0,*) "Decrease Tgridmin to at least ",Tmin
      goto 999
 903  write(0,*) "Increase Tgridmax to at least ",Tmax
      goto 999        
 999  return
      end