      program onedrock
      implicit none
      integer nr,i,nplot,niter,np
      parameter(nr=520)
      real px(nr),py(nr),px2(nr),py2(nr)
      double precision T(nr),T0,asemi,time,dtime,Fstar,L,Lsun,Pi,
     .  Lstar,Ab,dist,AU,sbolt,sarea(nr),dEn(nr),tcond(nr),drad,
     .  Cp(nr),M(nr),rho,fpi,rad(nr),rad0,REarth
      logical loop
      
C     Constants
      Lsun=3.839d26 !Solar luminosity [W]
      Pi=acos(-1.d0)   !Pi
      fPi=4.0d0*Pi ! 4Pi
      AU=1.49598d11 !Astronomical Unit [m]
      sbolt=5.670d-8 !Stefan-Boltzmann [W m-2 K-4]
      REarth=6.371*10d6 !Earth radius [m]
      
C     Model parameters      
      T0=1600.0d0 !temperature of inner boundary (core) [K]
      asemi=1.00d0 !distance to star [AU]
      dist=asemi*AU !distance to star [m]
      dtime=3.0d7 !change in time [seconds]
      Lstar=1.0d0 !star luminosity [Lsun]
      L=Lstar*Lsun !star luminosity [W]
      Ab=0.367d0 !albedo of surface
c      sarea=1.0d0 !surface area of cell [m^2]
      drad=250.0d0 !length of cell [m]
c      tcond=1.1d0 !thermal conductivity [W m-1 K-1]
c      Cp=0.84d0*1000.0d0 !specific heat capacity [J kg-1 K-1]
      rho=2700.0d0 !density of material [kg m-3]
      rad0=1.0d0*REarth
      
      nplot=10000      
      
C     Initialize data arrays
      call init(nr,T,T0,dEn,M,rad,rad0,rho,drad,sarea,fpi)
      time=0.0d0 !init time counter
      loop=.true. !init loop conditional
      niter=0 !counter iterations
      np=0 !number of plots made
     
      call opengraphics()
      
      do while(loop)
C       Step 1 - update time
        time=time+dtime
        niter=niter+1
      
        do 10 i=1,nr
            dEn(i)=0.0
 10     continue
            
C       Step 2a irradiate surface element
        Fstar=L/(4.0*Pi*dist*dist)*(1.0-Ab) ![J/s/m^2]
C       absorb only on half of sphere
        dEn(nr)=dEn(nr)+Fstar*dtime*sarea(nr)/2.0d0     ![J]
C       Step 2b radiate Energy away from surface
        dEn(nr)=dEn(nr)-dtime*sarea(nr)*sbolt*T(nr)**4.0d0 ![J]

C       update Cp and tcond
        call kcp(nR,tcond,Cp,T,rho) 

C       radiogenic heat production        
        call Arad(nR,dEn,drad,sarea,dtime)
        
C       thrusting
        if(time.lt.1.25d15) call As(nR,dEn,drad,sarea,dtime)

C       Step 3 apply heat conduction
        call cond(nr,dEn,T,sarea,drad,tcond,dtime)
 
C       Step 4 - update temperatures from heat capacity
        call hcap(nr,T,dEn,Cp,M)
c        write(6,500) time,(En(i),i=1,nr)
        if(mod(niter,nplot).eq.0)then
c            write(6,500) time/(24.0*60.0*60.0),(T(i),i=1,nr)
            call plot(nr,px,py,px2,py2,drad,T,np,Time)
        endif
c        read(5,*)
 500    format(F7.2,10(1X,1PE9.2))
      enddo
      
      call closegraphics()
      
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine kcp(nR,tcond,Cp,T,rho)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nR,i
      double precision tcond(nR),Cp(nR),T(nR),rho,kappa
     
      do 11 i=1,nR
C       calculate diffusivity and specific heat
C       [mm^2 s^-1]     [J mole^-1 K^-1]
        if(T(i).lt.846.0)then
            kappa=567.3/T(i)-0.062
            Cp(i)=199.50+0.0857*T(i)-5.0d-6/T(i)/T(i)
        else
            kappa=0.732-0.000135*T(i)
            Cp(i)=229.32+0.0323*T(i)-47.9d-6/T(i)/T(i)
        endif
C       average molar mass : 221.78 g mol^-1
        Cp(i)=Cp(i)*1000.0d0/221.78d0
        tcond(i)=kappa*1.0d-6*rho*Cp(i) ![W m-1 K-1]
 11   continue
 
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine As(nR,dEn,drad,sarea,dtime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nR,i
      double precision dEn(nR),drad,sarea(nR),depth,sEn,dtime,zc,tau,v,
     .  area,dz
      
      zc=35000.0 !center of sheer zone
      dz= 3000.0 !width of sheer zone
      tau=30.0d6 !sheer stress [Pa]
      v=9.506d-10 !thrusting velocity [m/s]
      
      do 10 i=1,nR
        depth=drad*dble(nR-i)
        area=sarea(i)*drad
        sEn=tau*v*exp(-(depth-zc)*(depth-zc)/(dz*dz))/dz
        dEn(i)=dEn(i)+sEn*area*dtime
 10   continue
      
      return
      end
                
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine Arad(nR,dEn,drad,sarea,dtime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nR,i
      double precision dEn(nR),drad,sarea(nR),depth,aEn,dtime,d1,d2,d3,
     .  ds,a1,a2,a3,area
      
      a1=2.0d-6 !/mu W m^-3
      a2=0.2d-6 !/mu W m^-3
      a3=0.02d-6 !/mu W m^-3
      d1=15000.0
      d2=35000.0
      d3=70000.0
      ds= 8500.0
      
      do 10 i=1,nR
        depth=drad*dble(nR-i)
        area=sarea(i)*drad
        if(depth.lt.d1)then
            aEn=a1
        elseif((depth.ge.d1).and.(depth.lt.d2))then
            aEn=a1*exp(-(depth-d1)/ds)
        elseif((depth.ge.d2).and.(depth.lt.d3))then
            aEn=a2
        else
            aEn=a3
        endif
        dEn(i)=dEn(i)+aEn*dtime*area
 10   continue
 
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plot(nr,px,py,px2,py2,drad,T,np,Time)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nr,i,np,nred
      real px(nr),py(nr),rdrad,px2(nr),py2(nr),R,G,B,dC
      double precision drad,T(nr),time,tyear
      character*80 ttext,ttextold
      
      np=np+1
      nred=10
      
      tyear=time/31556700.00
      rdrad=real(drad)
      do 10 i=1,nr
        py(i)=rdrad*real(nr-i)/1000.0
        px(i)=real(T(i))
 10   continue
c      write(6,*) "np:",np
      if((np.eq.1).or.(mod(np,nred).eq.0))then
        call panel(3)
        call pgwindow(0.0,2000.0,py(1),py(nR))
        call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
        call pglabel("Temperature (K)","Depth (km)","")
        write(ttext,500) tyear
        if(np.gt.1)then
            call pgscr(0,0.0,0.0,0.0)
            call pgsci(0)
            call pgptxt(1500.0,25.0,0.0,0.0,ttextold)
            call pgsci(1)
        endif
        call pgptxt(1500.0,25.0,0.0,0.0,ttext)
        ttextold=ttext
      endif
 500  format(F11.1)
      if(np.gt.1) then
        dc=0.7*mod(real(tyear/1.0d7),1.0)
        R=0.0+dc
        G=0.3+dc
        B=0.2+dc
        write(0,*) R,G,B
        call pgscr(0,R,G,B)
        call pgsci(0)
        call pgline(nr,px2,py2)
      endif
      call pgsci(2)
      call pgline(nr,px,py)
      call pgsci(1)
      do 11 i=1,nr
        px2(i)=px(i)
        py2(i)=py(i)
 11   continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine panel(n)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer n
      
      if (n.eq.1) then
        call pgsvp(0.1,0.95,0.83,0.98)
      elseif(n.eq.2) then
        call pgsvp(0.1,0.95,0.60,0.75)
      elseif(n.eq.3) then
        call pgsvp(0.1,0.50,0.10,0.50)
      elseif(n.eq.4) then
        call pgsvp(0.55,0.95,0.10,0.50)
      endif
      
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine opengraphics()
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Calls PGPLOT routines to open graphic window for plotting
      implicit none
      
      call pgopen('?')
c      call pgpap(6.0,1.0)
c      call pgask(.false.)
      call pgsch(0.8)
c      call pgsubp(1,4)
c      call pgvport(0.1,0.95,0.2,0.8) !gives enough room for labels
      
c      call pgpage()
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine closegraphics()
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     close up graphics device
      implicit none
      
      call pgclos()
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine hcap(nr,T,dEn,Cp,M)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nr,i
      double precision T(nr),dEn(nr),Cp(nR),M(nr)
      
      do 10 i=1,nr
        T(i)=T(i)+dEn(i)/Cp(i)/M(i)
 10   continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine cond(nr,dEn,T,sarea,drad,tcond,dtime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nr,i
      double precision dEn(nr),T(nr),sarea(nR),drad,tcond(nR),dtime
      
C     bottom element
      dEn(1)=dEn(1)+tcond(i)*sarea(1)*(T(2)-T(1))/drad*dtime ! [J]
C     middle layers
      do 10 i=2,nr-1
        dEn(i)=dEn(i)+tcond(i)*sarea(i)*(T(i+1)-T(i))/drad*dtime 
        dEn(i)=dEn(i)+tcond(i)*sarea(i)*(T(i-1)-T(i))/drad*dtime
 10   continue
C     top element
      dEn(nr)=dEn(nr)+tcond(i)*sarea(i)*(T(nr-1)-T(nr))/drad*dtime
      
      return
      end

      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine init(nr,T,T0,En,M,rad,rad0,rho,drad,sarea,fpi)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nr,i
      double precision T(nr),T0,En(nr),M(nr),rho,drad,sarea(nr),rad(nr),
     .  rad0,fpi,T1,T2,d1,d2,d3,depth
      
      T1=200.0d0
      T2=800.0d0
      d1=0.0d0
      d2=35000.0d0
      d3=drad*dble(nR-1)
      
      T(1)=T0  !inside temperature
      En(1)=0.0d0
      rad(1)=rad0-drad*dble(nR-1)
      sarea(1)=fpi*rad(1)*rad(1)
      M(1)=sarea(1)*drad*rho
      
      do 10 i=2,nR
        T(i)=T1  !temperature [K]
        depth=drad*dble(nR-i)
        if((depth.ge.d1).and.(depth.lt.d2))then
            T(i)=T1+(T2-T1)/(d2-d1)*depth
        else
            T(i)=T1+(T0-T1)/(d3-d2)*(depth-d2)
        endif
        En(i)=0.0d0 !Energy [J]
        rad(i)=rad0-dble(nR-i)*drad
        sarea(i)=fpi*rad(1)*rad(1)
        M(i)=sarea(i)*drad*rho !mass [kg]
 10   continue
 
      return
      end


