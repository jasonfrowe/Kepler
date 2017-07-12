CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function mandelagol(nintg,R1,R2,x1,x2,y1,y2,c,
     .  b0,mu,mulimb0,mulimbf,dist)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Adapted from Mandel and Agol, 2002, ApJ 580, L171
      implicit none
      integer i,nintg,sflag
      double precision R1,R2,x1,x2(nintg),y1,y2(nintg),c(4),
     .  c1,c2,c3,c4,mulimb0(nintg),mulimbf(nintg,5),rl,b0(nintg),
     .  mu(nintg),dist(nintg)
      
      mu=0

      c1=c(1)
      c2=c(2)
      c3=c(3)
      c4=c(4)
      rl=R2/R1
      sflag=0
      do 10 i=1,nintg
       dist(i)=Sqrt((x2(i)-x1)*(x2(i)-x1)+(y2(i)-y1)*(y2(i)-y1))/(R1+R2)
c        if(dist(i).ge.1.0d0)then
c            sflag=sflag+1
c            b0(i)=2.0
c        else
            b0(i)=(R1+R2)*dist(i)/R1
c        endif
 10   continue
      
c      write(6,500) "hello",(b0(i),i=1,nintg)
 500  format(A5,11(1X,F5.3))
      if(sflag.lt.nintg) then
c          call occultnl(rl,c1,c2,c3,c4,b0,mulimb0,mulimbf,nintg)
c          call occultsmall(rl,c1,c2,c3,c4,nintg,b0,mulimb0)
      endif
c      if(b0(1).le.1.0) write(6,*)b0(1),mu

C      mandelagol=mulimb0(1)
      mandelagol=0.0
      do 11 i=1,nintg
c        mandelagol=mandelagol+mu(i)
c        if(dist(i).ge.1.0d0) mulimb0(i)=1.0d0
        mandelagol=mandelagol+mulimb0(i)
 11   continue
      mandelagol=mandelagol/dble(nintg)
c      write(6,*) "hello2",mandelagol
   
      return
      end

      subroutine occultsmall(p,c1,c2,c3,c4,nz,z,mu)
      implicit none
      integer i,nz
c      parameter (nz=201)
      real*8 p,c1,c2,c3,c4,z(nz),mu(nz),i1,norm,
     &       x,tmp,iofr,pi
C This routine approximates the lightcurve for a small 
C planet. (See section 5 of Mandel & Agol (2002) for
C details):
C Input:
C  p      ratio of planet radius to stellar radius
C  c1-c4  non-linear limb-darkening coefficients
C  z      impact parameters (positive number normalized to stellar 
C        radius)- this is an array which MUST be input to the routine
C  NOTE:  nz must match the size of z & mu in calling routine
C Output:
C  mu     flux relative to unobscured source for each z
C
      pi=acos(-1.d0)
      norm=pi*(1.d0-c1/5.d0-c2/3.d0-3.d0*c3/7.d0-c4/2.d0)
      i1=1.d0-c1-c2-c3-c4
      do i=1,nz
        mu(i)=1.d0
        if(z(i).gt.1.d0-p.and.z(i).lt.1.d0+p) then
          x=1.d0-(z(i)-p)**2
          tmp=(1.d0-c1*(1.d0-0.8d0*x**0.25)
     &             -c2*(1.d0-2.d0/3.d0*x**0.5)
     &             -c3*(1.d0-4.d0/7.d0*x**0.75)
     &             -c4*(1.d0-0.5d0*x))
          mu(i)=1.d0-tmp*(p**2*acos((z(i)-1.d0)/p)
     &        -(z(i)-1.d0)*sqrt(p**2-(z(i)-1.d0)**2))/norm
        endif
        if(z(i).le.1.d0-p.and.z(i).ne.0.d0) then
          mu(i)=1.d0-pi*p**2*iofr(c1,c2,c3,c4,z(i),p)/norm
        endif
        if(z(i).eq.0.d0) then
          mu(i)=1.d0-pi*p**2/norm
        endif
      enddo
      return
      end

      function iofr(c1,c2,c3,c4,r,p)
      implicit none
      real*8 r,p,c1,c2,c3,c4,sig1,sig2,iofr
      sig1=sqrt(sqrt(1.d0-(r-p)**2))
      sig2=sqrt(sqrt(1.d0-(r+p)**2))
      iofr=1.d0-c1*(1.d0+(sig2**5-sig1**5)/5.d0/p/r)
     &         -c2*(1.d0+(sig2**6-sig1**6)/6.d0/p/r)
     &         -c3*(1.d0+(sig2**7-sig1**7)/7.d0/p/r)
     &         -c4*(p**2+r**2)
      return
      end

      subroutine occultnl(rl,c1,c2,c3,c4,b0,mulimb0,mulimbf,nb)
c; Please cite Mandel & Agol (2002) if making use of this routine.
      implicit none
      integer i,j,nb,nr,i1,i2,nmax
      parameter (nmax=2**16)
      real*8 mulimbf(nb,5),pi,c1,c2,c3,c4,rl,bt0(nb),b0(nb),mulimb0(nb),
     &       mulimb(nb),mulimbp(nb),dt,t(nmax+1),th(nmax+1),r(nmax+1),
     &       sig,mulimb1(nb),mulimbhalf(nb),mulimb3half(nb),mulimb2(nb),
     &       sig1,sig2,omega,dmumax,fac,mu(nb),f1,f2
      pi=acos(-1.d0)
C  This routine uses the results for a uniform source to
C  compute the lightcurve for a limb-darkened source
C  (5-1-02 notes)
C Input:
C   rl        radius of the lens   in units of the source radius
C   c1-c4     limb-darkening coefficients
C   b0        impact parameter normalized to source radius
C Output:
C  mulimb0 limb-darkened magnification
C  mulimbf lightcurves for each component
C  
C  First, make grid in radius:
C  Call magnification of uniform source:
      call occultuniform(b0,rl,mulimb0,nb)
      i1=nb
      i2=1
      fac=0.d0
      do i=1,nb
        bt0(i)=b0(i)
        mulimbf(i,1)=1.d0
        mulimbf(i,2)=0.8d0
        mulimbf(i,3)=2.d0/3.d0
        mulimbf(i,4)=4.d0/7.d0
        mulimbf(i,5)=0.5d0
        mulimb(i)=mulimb0(i)
        if(mulimb0(i).ne.1.d0) then
          i1=min(i1,i)
          i2=max(i2,i)
        endif
        fac=max(fac,abs(mulimb0(i)-1.d0))
      enddo
C print,rl
      omega=4.*((1.d0-c1-c2-c3-c4)/4.+c1/5.+c2/6.+c3/7.+c4/8.)
      nr=2
      dmumax=1.d0
c      write(6,*) 'i1,i2 ',i1,i2
      do while (dmumax.gt.fac*1.d-3)
        do i=i1,i2
          mulimbp(i)=mulimb(i)
        enddo
        nr=nr*2
c        write(6,*) 'nr ',nr
        dt=0.5d0*pi/dble(nr)
        if(nr+1.gt.nmax) then
         write(0,*) "M&A Seg: ",nr+1,b0(1)
         return
        endif
        do j=1,nr+1
          t(j) =dt*dble(j-1)
          th(j)=t(j)+0.5d0*dt
          r(j)=sin(t(j))
        enddo
        sig=sqrt(cos(th(nr)))
        do i=i1,i2
          mulimbhalf(i) =sig**3*mulimb0(i)/(1.d0-r(nr))
          mulimb1(i)    =sig**4*mulimb0(i)/(1.d0-r(nr))
          mulimb3half(i)=sig**5*mulimb0(i)/(1.d0-r(nr))
          mulimb2(i)    =sig**6*mulimb0(i)/(1.d0-r(nr))
        enddo
        do j=2,nr
          do i=1,nb
            b0(i)=bt0(i)/r(j)
          enddo
C  Calculate uniform magnification at intermediate radii:
          call occultuniform(b0,rl/r(j),mu,nb)
C  Equation (29):
          sig1=sqrt(cos(th(j-1)))
          sig2=sqrt(cos(th(j)))
          dmumax=0.d0
          do i=i1,i2
            f1=r(j)*r(j)*mu(i)/(r(j)-r(j-1))
            f2=r(j)*r(j)*mu(i)/(r(j+1)-r(j))
            mulimbhalf(i) =mulimbhalf(i) +f1*sig1**3-f2*sig2**3
            mulimb1(i)    =mulimb1(i)    +f1*sig1**4-f2*sig2**4
            mulimb3half(i)=mulimb3half(i)+f1*sig1**5-f2*sig2**5
            mulimb2(i)    =mulimb2(i)    +f1*sig1**6-f2*sig2**6
            mulimb(i)=((1.d0-c1-c2-c3-c4)*mulimb0(i)+c1*mulimbhalf(i)*dt
     &        +c2*mulimb1(i)*dt+c3*mulimb3half(i)*dt+c4*mulimb2(i)*dt)
     &        /omega
            if(mulimb(i)+mulimbp(i).ne.0.d0) then 
              dmumax=max(dmumax,abs(mulimb(i)-mulimbp(i))/(mulimb(i)+
     &               mulimbp(i)))
            endif
          enddo
        enddo
      enddo
      do i=i1,i2
        mulimbf(i,1)=mulimb0(i)
        mulimbf(i,2)=mulimbhalf(i)*dt
        mulimbf(i,3)=mulimb1(i)*dt
        mulimbf(i,4)=mulimb3half(i)*dt
        mulimbf(i,5)=mulimb2(i)*dt
        mulimb0(i)=mulimb(i)
      enddo
      do i=1,nb
        b0(i)=bt0(i)
      enddo
      return
      end
      
      subroutine occultuniform(b0,w,muo1,nb)
      implicit none
      integer i,nb
      real*8 muo1(nb),w,b0(nb),z,pi,lambdae,kap0,kap1
      if(abs(w-0.5d0).lt.1.d-3) w=0.5d0
      pi=acos(-1.d0)
C  This routine computes the lightcurve for occultation
C  of a uniform source without microlensing  (Mandel & Agol 2002).
C Input:
C 
C  rs   radius of the source (set to unity)
C  b0   impact parameter in units of rs
C  w    occulting star size in units of rs
C 
C Output:
C  muo1 fraction of flux at each b0 for a uniform source
C 
C  Now, compute pure occultation curve:
      do i=1,nb
C  substitute z=b0(i) to shorten expressions
        z=b0(i)
C  the source is unocculted:
C  Table 3, I.
        if(z.ge.1.d0+w) then
          muo1(i)=1.d0
          goto 1
        endif
C  the  source is completely occulted:
C  Table 3, II.
        if(w.ge.1.d0.and.z.le.w-1.d0) then
          muo1(i)=0.d0
          goto 1
        endif
C  the source is partly occulted and the occulting object crosses the limb:
C  Equation (26):
        if(z.ge.abs(1.d0-w).and.z.le.1.d0+w) then
          kap1=acos(min((1.d0-w*w+z*z)/2.d0/z,1.d0))
          kap0=acos(min((w*w+z*z-1.d0)/2.d0/w/z,1.d0))
          lambdae=w*w*kap0+kap1
          lambdae=(lambdae-0.5d0*sqrt(max(4.d0*z*z-(1.d0+z*z-w*w)**2,
     &            0.d0)))/pi
          muo1(i)=1.d0-lambdae
        endif
C  the occulting object transits the source star (but doesn't
C  completely cover it):
        if(z.le.1.d0-w) muo1(i)=1.d0-w*w
 1      continue
      enddo
C muo1=1.d0-lambdae
      return
      end
      
