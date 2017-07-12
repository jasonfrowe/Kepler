CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program wdchandra
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Calculate polytropic WD models
      implicit none      
      integer n,iargc
      double precision Pi,me,c,h,mu,ue,C1,C2,Msun,Rsun,onedzc,y(2),
     .  dydx(2),htry,eps,hdid,hnext,yscal(2),z,x,zc,zeta1,dphidzeta1,
     .  alpha,G,radius,mass,onedzc2
      character*80 cline
      common /fitting/ onedzc
      
      Pi=acos(-1.d0)    !Pi
      me=9.10938188d-31 !electron mass (kg)
      mu=1.66053886d-27 !atomic mass unit (kg)
      c =299792458d0   !speed of light (m/s)
      h =6.626068d-34   !Planck constant (m^2 kg / s)    
      Msun=1.9891d30 !kg  mass of Sun
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      
      ue=2.0d0  !electron density
      
      C1=Pi*me**4*c**5/(3.0d0*h**3)
      C2=8.0d0*Pi*ue*mu*me**3*c**3/(3.0d0*h**3)
      
c      write(0,*) "C1:",C1
c      write(0,*) "C2:",C2
      
      if(iargc().lt.1) goto 901
      call getarg(1,cline)
      read(cline,*) onedzc2
      if((onedzc2.le.0).or.(onedzc2.ge.1.0)) goto 902
      
C     define model onedzc=(0,inf)
c      onedzc=sqrt(0.1)
      onedzc=sqrt(onedzc2)
      
      n=2
      y(1)=0.0d0! dphi/dzeta=a (initial conditions 35.9)
      y(2)=1.0d0! phi
      
      zc=1.0d0/onedzc !central z
      x=1.0d-30 !zeta start
      
      dydx(2)=y(1)
      dydx(1)=2.0*y(1)/x-(y(2)*y(2)-onedzc*onedzc)**(3.0d0/2.0d0)

      htry=1.0d-6 !why not try this.
      eps=1.0d-14 !accuracy
      yscal(1)=1.0d0
      yscal(2)=1.0d0
      z=zc
c      write(0,*) "z:",z
      
      do 10 while(z.gt.1.0d0)
        zeta1=x
        dphidzeta1=y(1)
        call derivs(x,y,dydx)
        call rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext)
        htry=hnext
        z=y(2)*zc
c        write(6,*) x,y(1),z
c        write(6,*) hdid
c        write(6,*) hdid,hnext
c        read(5,*)
 10   continue
 
      alpha=Sqrt(2.0d0*C1/(Pi*G))/(C2*zc)
      Radius=alpha*zeta1
      Mass=(4.0d0*Pi/(C2*C2))*(2.0d0*C1/(Pi*G))**(3.0d0/2.0d0)*
     .  (-zeta1*zeta1*dphidzeta1)
c      write(6,*) "alpha:    ",alpha
c      write(6,*) "zeta:     ",zeta1
c      write(6,*) "-z^2dPdz: ",-zeta1*zeta1*dphidzeta1
c      write(6,*) Mass/Msun*ue**2,Radius/1.0d3*ue
      write(6,*) Mass/Msun,Radius/Rsun
      
      goto 999
 901  write(6,*) "Usage: wdchandra [1/zc^2]"
      write(6,*) "   where 0 < 1/zc^2 < 1  "
 902  write(6,*) " 1/zc^2 must greater than 0 and less than 1"
      goto 999
 999  end