CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine transitmodel(nrad,ntheta,nfit,nplanet,sol,npt,time,
     .  exptime,dtype,tmodel,rbb,nplot)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nrad,ntheta,np,col,nplanet,nfit,npt,nplot,dtype(npt)
      real rbb(3,4)
      double precision F0,TFlux,f,feff,obl,logg,omega,Teff,beta,
     .  sol(nfit),Fp,time(npt),tmodel(npt),exptime(npt),rhostar


      np=1 !planet number
      col=11*(np-1)+15
      
      rhostar=sol(1) !mean stellar density
      f    =sol(12)  !fpole/feq
      obl  =sol(13)  !obquity 
      logg =sol(10)  !log(g) (star)
      Teff =sol(9)   !Teff (K) (star)
      omega=sol(14)  !stellar rotation
      beta =sol(15)  !gravity darkening 
      
      feff=1.0d0-sqrt((1.0d0-f)*(1.0d0-f)*cos(obl)*cos(obl)+
     .  sin(obl)*sin(obl))
      
      F0=TFlux(nrad,ntheta,feff,f,obl,logg,omega,Teff,beta,rhostar,rbb,
     .  nplot)
c      write(0,*) "F0:",F0
      
      call PFlux(nrad,ntheta,nfit,nplanet,sol,npt,time,exptime,dtype,F0,
     .  tmodel,rbb,nplot)
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine Pflux(nrad,ntheta,nfit,np,sol,npt,time,exptime,dtype,
     .  F0,tmodel,rbb,nplot)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nrad,ntheta,nfit,np,col,i,j,ii,ntheta2,npt,k,nintg,jj,
     .  nplot,kk,dtype(npt)
      parameter(nintg=11)
      real rbb(3,4)
      double precision sol(nfit),b,rdr,xc,yc,dR,R,twon,ringarea,Pi,
     .  twoPi,dtheta,dntheta,x,y,xp,yp,f,obl,logg,Teff,omega,beta,feff,
     .  In,Intensity,u,limb,nl(4),zonearea,theta,Rp,Per,adrs,ecw,
     .  esw,eccn,w,G,dist,time(npt),tmodel(npt),dnintg,dnintgm1,Eanom,
     .  epoch,Manom,phi0,trueanomaly,zc,zed,phase,incl,yreal,rot(2),
     .  rangle,t(nintg),exptime(npt),phi(nintg),Tanom(nintg),
     .  arad(nintg),x2(nintg),y2(nintg),distance,F0,phi1,ell,pid2,ag,
     .  alb,albedomod,ted,parea,vt,Kr,dop,Cs,fDB,dn(2),rhostar,Vr

      Pi=acos(-1.d0)!define Pi
      twoPi=2.0d0*Pi
      pid2=Pi/2.0d0
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      Cs=2.99792458e8 !Speed of light
      fDB=-1.896 !Doppler Boosting factor

C     Stellar parameters
      f    =sol(12)
      obl  =sol(13)
      logg =sol(10)
      Teff =sol(9)
      omega=sol(14)
      beta =sol(15)
      feff=1.0d0-sqrt((1.0d0-f)*(1.0d0-f)*cos(obl)*cos(obl)+
     .  sin(obl)*sin(obl))
      rhostar=sol(1)

      nl(1)=sol(2)  !limb darkening
      nl(2)=sol(3)
      nl(3)=sol(4)
      nl(4)=sol(5)
      
      dnintg=dble(nintg) !convert integer to double
      dnintgm1=2.0*dnintg-2.0

CCCCCCCCC start multiplanet CCCCCCCC
      do 12 k=1,npt !walk through each time-step
       tmodel(k)=0.0d0 !initialize flux for each time step
       do 16 kk=1,np  !loop through each planet
      
       col=11*(kk-1)+15 !offset to get planet parameters

       per=sol(col+2) !period in days
       adrs=1000.0*sol(1)*G*(Per*86400.0d0)**2/(3.0d0*Pi)
       adrs=adrs**(1.0d0/3.0d0) !a/R*
c      write(0,*) 'adrs:',adrs
       b=abs(sol(col+3))  !impact parameter
       rdr=abs(sol(col+4)) !scaled planetary radius
       ag=sol(col+10) !amplitude of thermal/reflection lightcurve.
       ted=sol(col+8)*1.0d-6*F0 !eclipse depth
       parea=pi*rdr*rdr
       Kr=sol(col+7)!abs(sol(col+7)) !radial velocity (m/s)
       rangle=sol(col+11)*pi/180.0 !rotation of orbital plane 
      
       ecw=sol(col+5)
       esw=sol(col+6)
       eccn=sqrt(ecw*ecw+esw*esw) !eccentricity
       if(eccn.ge.1.0) eccn=0.99
       if(eccn.eq.0.0d0)then
          w=0.0d0
       else
          if(ecw.eq.0.0d0)then
            w=Pi/2.0d0
          else
            w=atan(esw/ecw)
          endif
          if((ecw.gt.0.0d0).and.(esw.lt.0.0d0))then
              w=twoPi+w
          elseif((ecw.lt.0.0d0).and.(esw.ge.0.0d0))then 
              w=Pi+w
          elseif((ecw.le.0.0d0).and.(esw.lt.0.0d0))then
              w=Pi+w
          endif
       endif      
c      w=Pi/2.0d0
      
C      Find phase at centre of transit
       epoch=sol(col+1)   !center of transit time (days)
       phi1=(epoch/per-int(epoch/per))*twopi !mean anonomaly of transit
       Eanom=tan(w/2.0d0)/sqrt((1.0d0+eccn)/(1.0d0-eccn))
       Eanom=2.0d0*atan(Eanom) !ecc anomaly
       phi0=Eanom-eccn*sin(Eanom) !phase offset to center transit
      
C      Find inclination
       incl=acos(b/adrs)      

C       set up arrays for sub-time steps to integrate exposure time
        do 13 jj=1,nintg  !get array of times spanning exposure time
C           These times are centered on time(i)
            t(jj)=time(k)+exptime(k)*(2.0*dble(jj)-dnintg-1.0)/dnintgm1-
     .          epoch
            phi(jj)=t(jj)/per-floor(t(jj)/per)
            phi(jj)=phi(jj)*twoPi+phi0  
            Manom=phi(jj) !mean anomaly
            if(Manom.gt.twoPi) Manom=Manom-twoPi
            if(Manom.lt.0.0d0) Manom=Manom+twoPi
            call kepler(Manom,Eanom,eccn) !eccentric anomaly 
            Tanom(jj)=trueanomaly(eccn,Eanom) !get true anonaly
            if(phi(jj).gt.Pi) phi(jj)=phi(jj)-twoPi            
            arad(jj)=distance(adrs,eccn,Tanom(jj)) 
c            arad(jj)=(x2(jj)*x2(jj)+y2(jj)+y2(jj)) !star-planet sep.
            x2(jj)=arad(jj)*Sin(Tanom(jj)-w)  !x position of planet
            y2(jj)=arad(jj)*Cos(Tanom(jj)-w)  !y position of planet

 13     continue !end sub-time array assignment

C       Plotting
        if(nplot.eq.1)then !nplot flags whether to plot
            xc=x2(nintg/2+1)    !get mid-time of x,y position
            yreal=y2(nintg/2+1)      
            call pgpanl(2,1)  !switch to correct panel
            call pgwindow(rbb(2,1),rbb(2,2),rbb(2,3),rbb(2,4))
            call pgsci(kk+1) !choose colour
c            write(0,*) xc,yc
            call pgpt1(real(xc),real(yreal),-1)  !plot position.
c            call pgcirc(real(xc),real(yreal),real(rdr)) !plot position
            call pgsci(1) !revert back to default colour
            yc=yreal*cos(incl)
            rot(1)= xc*cos(rangle)+yc*sin(rangle)
            rot(2)=-xc*sin(rangle)+yc*cos(rangle)
            xc=rot(1)
            yc=rot(2)
        endif !end of plotting

c        tmodel(k)=0.0d0
        ell=0.0d0  !initialize sums for ellipsoidal
        alb=0.0d0  !                    reflection/emission
        vt=0.0d0   !                    radial velocity
        dop=0.0d0  !                    Doppler boosting
        do 15 jj=1,nintg !loop over sub-time-steps
            xc=x2(jj)    !x-position of orbit 
            yreal=y2(jj) !y-position of orbit
            yc=y2(jj)*cos(incl) !projected y-position of orbit
            rot(1)= xc*cos(rangle)+yc*sin(rangle) !rotated on sky
            rot(2)=-xc*sin(rangle)+yc*cos(rangle) !rotated on sky
            xc=rot(1) !projected x,y positions
            yc=rot(2)

C           plotting transit position
            if((abs(xc).le.1.1+rdr).and.(abs(yc).le.1.1+rdr)
     .        .and.(yreal.ge.0.0))then  !if position is plotable.. 
                if(nplot.eq.1)then
                    call pgpanl(1,1) !choose panel
                    call pgwindow(rbb(1,1),rbb(1,2),rbb(1,3),rbb(1,4))
                    call pgsci(kk+1) !choose colour
c                    call pgcirc(real(xc),real(yc),real(rdr))
                    call pgpt1(real(xc),real(yc+rdr),1) !plot position
                    call pgpt1(real(xc),real(yc-rdr),1)
                    call pgsci(1) !default colour
                endif
            endif !end of plotting 

C           Radial velocities
            Vr=Kr*(cos(Tanom(jj)-w+pid2)+eccn*cos(-w+pid2))
            vt=vt+Vr

CCCCCCCCCCC Determine whether we have photometry or RVs

C           Ellipsoidal variations
            ell=ell+sol(col+9)*(arad(jj)/adrs)**(1.0d0/3.0d0)*
     .          cos(2.0d0*(Pid2+phi(jj)))
     
C           Phasecurve - assumes Lambertian model
            alb=alb+albedomod(Pi,ag,phi(jj))*adrs/arad(jj)
            
C           Doppler-boosting, controlled by fDB
            dop=dop+fDB*Vr/Cs

C           check to see if there could be a transit
            dist=sqrt(xc*xc+yc*yc) !transformed space

C           transformed space distance
            dist=sqrt(xc*xc+(yc/(1.0d0-feff))**2.0d0)  
            if(dist.le.1.0+rdr)then  !check for transit/occultation

                dR=rdr/dble(nrad)/2.0d0 !number of divisions in radius  
                twon=dble(2*nrad) !preconvert to dble
      
                do 10 i=1,nrad !loop over radius divisions
                    R=rdr*dble(2*i-1)/twon !R (ring for integration)
                    !piR2^2 - piR1^2
                    ringarea=pi*((R+dR)*(R+dR)-(R-dR)*(R-dR)) 
c                  write(0,*) R,dR,ringarea
                    ntheta2=max(int(dble(ntheta)*R),10) !ring divisions
                    dtheta=2.0d0*Pi/dble(ntheta2) !angular size
                    dntheta=dble(ntheta2) !preconvert to double
                    do 11 j=1,ntheta2
                        theta=twopi*dble(j)/dntheta !theta
                        zonearea=ringarea*dtheta/twopi !dtheta/2pi
                        x=R*cos(theta)+xc !x-pos of surface element
                        y=R*sin(theta)+yc !y-pos of surface element
                        xp=x              !transformed
                        yp=y/(1.0d0-feff) !transformed 
c            write(6,*) xp,yp,theta

                        Rp=sqrt(xp*xp+yp*yp) !transformed distance
                        if(Rp.le.1.0d0)then  !check for transit/occult

                            if(yreal.ge.0.0)then !if y>=0 then transit
C                               call gravity darkening routine
                                In=Intensity(xp,yp,feff,f,obl,logg,
     .                              omega,Teff,beta,rhostar)
c                          write(6,*) xp,yp,In
            
                                u=sqrt(1-(xp*xp+yp*yp)) !limb-darkening
                                limb=0.0d0
                                do 14 ii=1,4
                                    limb=limb+nl(ii)*(1.0d0-u**
     .                                  dble(ii/2.0d0))
 14                             continue
                                limb=1.0d0-limb !amount of limb-darken

C                               Total flux obscured by element            
                                tmodel(k)=tmodel(k)+zonearea*In*limb

                            else  !if y<0 then occult

C                               The occultation is modeled here
                                tmodel(k)=tmodel(k)+zonearea*ted/parea
                            
                            endif !end of transit-or-occult flux est.
                             
                            
                        endif !end of in-transit/-in-occult logic
            
 11                 continue !end of loop over theta
 10             continue !end of loop over radius
                
            else 
                tmodel(k)=tmodel(k)+0.0d0 !no occultation or transit
            endif !end flux-integration logic
 15     continue  !end of loop for integration time (sub-time steps)
 16    continue !end of loop for each planet
 
       vt=vt/dnintg          !integration time averaged RV
       ell=ell/dnintg*1.0d-6 !convert from ppm to flux and average
       alb=alb/dnintg*1.0d-6 !convert from ppm to flux and average

       if(dtype(k).eq.0)then        
        dn(1)=F0+parea*ted
        dn(2)=dn(1)*sol(6)/(1.0d0-sol(6))  !add dilution
        tmodel(k)=(dn(1)+dn(2)-tmodel(k)/dnintg)/(dn(1)+dn(2))
        dn(1)=ell+alb+sol(8)+dop/dnintg 
        dn(2)=dn(1)*sol(6)/(1.0d0-sol(6)) !add dilution
        tmodel(k)=tmodel(k)+dn(1)+dn(2) !return model-flux
       elseif(dtype(k).eq.1)then
        tmodel(k)=vt
       endif
 12   continue !end of loop for each time-step

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function albedomod(Pi,ag,phi)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision Pi,phi,alpha,phase,ag

      phi=phi+Pi
      if(phi.gt.2.0*Pi) phi=phi-2.0*Pi


      alpha=abs(phi)      
c      alpha=2.0*Pi*t/Per+phi
      alpha=alpha-2.0*Pi*int(alpha/(2.0*Pi))
      if(alpha.gt.Pi) alpha=abs(alpha-2.0*pi)
c      write(6,*) t,alpha
c      phase=(1.0d0+cos(alpha))/2.0d0
      phase=(sin(alpha)+(Pi-alpha)*cos(alpha))/Pi  !Lambertian Sphere
      
      albedomod=ag*phase
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function TFlux(nrad,ntheta,feff,f,obl,logg,omega,
     .  Teff,beta,rhostar,rbb,nplot)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nrad,ntheta,i,j,CI1,CI2,icolour,k,ii,ntheta2,nplot
      real rx(nrad*ntheta),ry(nrad*ntheta),rz(nrad*ntheta),minz,maxz,
     .  rbb(3,4)
      double precision ringarea,dR,dtheta,Pi,R,twon,zonearea,Intensity,
     .  twopi,theta,xp,yp,feff,f,obl,logg,omega,dntheta,Teff,beta,In,
     .  nl(4),limb,u,T,Rp,rhostar
     
      nl(1)=0.5199
      nl(2)=0.1761
      nl(3)=-0.0091
      nl(4)=-0.0411

      CI1 = 16 
      CI2 = 86
      
      Pi=acos(-1.d0)!define Pi
      twoPi=2.0d0*Pi
      dR=1.0/dble(nrad)/2.0d0 !number of divisions in radius
c      dtheta=2.0d0*Pi/dble(ntheta) !number divisions in theta
      twon=dble(2*nrad) !preconvert to dble
c      dntheta=dble(ntheta)
      
      
      TFlux=0.0d0 !Initialize sum 
      k=0
      do 10 i=1,nrad
        R=1.0d0*dble(2*i-1)/twon !R
        ringarea=pi*((R+dR)*(R+dR)-(R-dR)*(R-dR)) !piR2^2 - piR1^2
c        write(0,*) R,dR,ringarea
        ntheta2=max(int(dble(ntheta)*R),10)
        dtheta=2.0d0*Pi/dble(ntheta2)
        dntheta=dble(ntheta2)
        do 11 j=1,ntheta2
            theta=twopi*dble(j)/dntheta !theta
            zonearea=ringarea*dtheta/twopi !dtheta/2pi
            xp=R*cos(theta)
            yp=R*sin(theta)
c            write(6,*) xp,yp,theta
            In=Intensity(xp,yp,feff,f,obl,logg,omega,Teff,beta,rhostar)
            
            u=sqrt(1-(xp*xp+yp*yp))
            limb=0.0d0
            do 14 ii=1,4
                limb=limb+nl(ii)*(1.0d0-u**dble(ii/2.0d0))
 14         continue
            limb=1.0d0-limb
            
            TFlux=TFlux+zonearea*In*limb
            
            k=k+1
            rx(k)=real(xp)
            ry(k)=real((1.0d0-feff)*yp)
            rz(k)=real(In*limb)
            
 11     continue
 10   continue

C     Now we go ahead and plot what we have..   
      if(nplot.eq.1)then
        call pgpanl(1,1) 
        call pgwindow(rbb(1,1),rbb(1,2),rbb(1,3),rbb(1,4))  
        minz=rz(1)
        maxz=rz(1)
        do 12 i=2,k
            minz=min(rz(i),minz) !scan for max-min range of intensity
            maxz=max(rz(i),maxz)
 12     continue
c      write(0,*) "maxmin",minz,maxz
        call pgbbuf() !we turn on graphics buffer for speed!
        do 13 i=1,k
            icolour=int((CI2-CI1)*(rz(i)-minz)/(maxz-minz))+CI1 
            call pgsci(icolour) !change colour
            if(ntheta.lt.200)call pgsch(3.0)
            call pgpt1(rx(i),ry(i),17) !plot point
            call pgsch(1.5)
 13     continue
        call pgebuf()
        call pgsci(1)
      endif
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function Intensity(xp,yp,feff,f,obl,logg,omega,
     .  Teff,beta,rhostar)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision xp,yp,x,y,z,feff,zed,f,obl,x0,y0,z0,R,gh(3),logg,
     .  Rt,omega,om2,g,Teff,T,beta,Fratio,grav,rhostar,Gr,fourthirdspi,
     .  pi,gpole
      
      Pi=acos(-1.d0)!define Pi
      fourthirdspi=4.0d0/3.0d0*pi
      om2=omega*omega
      Gr=6.674d-11 !N m^2 kg^-2  Gravitation constant
      
      x=xp
      y=(1.0d0-feff)*yp
      z=zed(x,y,f,obl)
c      Rp=sqrt(xp*xp+yp*yp)

c      write(6,500) x,y,z,Rp
 500  format(4(F6.3,1X))

      x0= x
      y0= y*cos(obl)+z*sin(obl)
      z0=-y*sin(obl)+z*cos(obl)
c      write(0,*) sqrt(x0*x0+y0*y0+z0*z0)
      
      R=sqrt(x0*x0+y0*y0+z0*z0)
      Rt=sqrt(x0*x0+z0*z0)
      grav=10.d0**logg/(R*R)
      gh(1)=-Gr*fourthirdspi*rhostar*1.0d3*x0/R + om2*x0/Rt
      gh(2)=-Gr*fourthirdspi*rhostar*1.0d3*y0/R
      gh(3)=-Gr*fourthirdspi*rhostar*1.0d3*z0/R + om2*z0/Rt
c      gh(1)=-grav*x0/R+om2*Rt*x0/Rt
c      gh(2)=-grav*y0/R
c      gh(3)=-grav*z0/R+om2*Rt*z0/Rt

      g=sqrt(gh(1)*gh(1)+gh(2)*gh(2)+gh(3)*gh(3))
      gpole=Gr*fourthirdspi*rhostar*1.0d3
      T=Teff*g**beta/(gpole**beta)
c      write(6,*) x0,y0,z0
c      write(6,*) gh(1),gh(2),gh(3)
c      write(6,*) g,gpole,T
c      write(6,*) '---'
c      read(5,*)
c      T=Teff*g**beta/((10.0d0**logg)**beta)
      Fratio=(T/Teff)**4.0d0 !bolometric.. for now.
c      write(6,*) log10(g),z,R
c      write(6,*) x,y,T

      Intensity=Fratio

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function zed(x,y,f,obl)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision x,y,f,obl,temp(3),sin2obl,cos2obl,Req,det,R2
      
      Req=1.0d0
      R2=Req*Req
      sin2obl=sin(obl)*sin(obl)
      cos2obl=cos(obl)*cos(obl)
      
      temp(1)=-f*f*R2*cos2obl+f*f*x*x*cos2obl+2.0d0*f*R2*cos2obl-
     .  2.0d0*f*x*x*cos2obl-R2*sin2obl-R2*sin2obl-R2*cos2obl+
     .  x*x*sin2obl+x*x*cos2obl+y*y*sin2obl*sin2obl+y*y*cos2obl*cos2obl+
     .  2.0d0*y*y*sin2obl*cos2obl
      temp(2)=f*f*y*sin(obl)*cos(obl)-2.0d0*f*y*sin(obl)*cos(obl)
      temp(3)=(f*f*cos2obl-2.0d0*f*cos2obl+sin2obl+cos2obl)
      zed=(sqrt(-(f-1)*(f-1)*temp(1))+temp(2))/temp(3)
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function zedold(x,y,f,obl)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision x,y,f,obl,temp(3),sin2obl,cos2obl,Req,det
      
      Req=1.0d0
      
      sin2obl=sin(obl)*sin(obl)
      cos2obl=cos(obl)*cos(obl)
      
      temp(1)=4.0d0*y*y*(1.0d0-(1.0d0-f*f))**2.0d0*sin2obl*cos2obl
    
      temp(2)=cos2obl*(1.0d0-f)*(1.0d0-f)+sin2obl
     
      temp(3)=(y*y*sin2obl-Req*Req+x*x)*(1.0d0-f*f)+y*y*cos2obl
     
      det=temp(1)-4.0d0*temp(2)*temp(3)
      
      temp(1)=-2.0d0*y*(1.0d0-(1.0d0-f)*(1.0d0-f))*sin(obl)*cos(obl)+
     .  sqrt(det)
      temp(2)=2.0d0*((1.0d0-f)*(1.0d0-f)*cos2obl+sin2obl)
      
      zedold=temp(1)/temp(2)
      
      return
      end
      