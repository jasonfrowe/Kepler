      subroutine occultnlv2(p,c1,c2,c3,c4,npt,z,flux)
      implicit none
      integer i,j,npt
      double precision p,c1,c2,c3,c4,z(npt),flux(npt),F,a,b,fN,omega,
     .  cn(4),pi,tpi
     
      pi=acos(-1.d0)
      tpi=2.0d0*pi
     
     
      cn(1)=c1
      cn(2)=c2
      cn(3)=c3
      cn(4)=c4
      omega=(1.0d0-cn(1)-cn(2)-cn(3)-cn(4))/4.0d0+
     .  cn(1)/5.0d0+cn(2)/6.0d0+cn(3)/7.0d0+cn(4)/8.0d0
      
      do 10 i=1,npt
      
        a=(z(i)-p)*(z(i)-p)
        b=(z(i)+p)*(z(i)+p)
      
C     There are 11 cases
        F=0.0d0
C     Case 1
        if((p.eq.0.0d0).or.(z(i).ge.1.0d0+p))then
            F=1.0d0
C     Case 2
        elseif((p.gt.0.0d0).and.(z(i).gt.0.5d0+abs(p-0.5d0)).and.
     .      (z(i).lt.1.0d0+p))then
            F=1.0d0
            do 11 j=1,4
                F=F-1.0d0/(tpi*omega)*fN(a,b,j,z(i),p)*cn(j)/
     .              dble(j+4)
 11         continue
        endif
 
        flux(i)=F
        write(0,*) z(i),flux(i)
 10   continue
           
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function fN(a,b,n,z,p)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer n
      double precision a,b,z,p,temp(5),pe,beta
      
      pe=dble(n+6)/4.0d0
      temp(1)=(1.0d0-a)**pe/sqrt(b-a)
      
      pe=dble(n+8)/4.0d0
      temp(2)=beta(pe,0.5d0)
      
      temp(3)=(z*z-p*p)/a
      
      fN=temp(1)*temp(2)*(temp(3))
      
      write(6,*) "fn: ",fN
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      double precision FUNCTION beta(z,w)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision w,z
C     USES gammln
C     Returns the value of the beta function B(z,w).
      double precision gammln
      beta=exp(gammln(z)+gammln(w)-gammln(z+w))
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   
      double precision FUNCTION gammln(xx)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision xx
C     Returns the value ln[Γ(xx)] for xx > 0.
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
C     Internal arithmetic will be done in double precision, a nicety 
C     that you can omit if five-figure accuracy is good enough.
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     * 24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     * -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
 11   continue
      gammln=tmp+log(stp*ser/x)
      return
      END
      
    !-------------------------------------------------------------------
    ! function  : f1 (complex*16)
    !
    ! package   : F1
    !
    ! Language  : Fortran 90
    !
    ! author    : F. Colavecchia (flavioc@lanl.gov)
    !                            (flavioc@cab.cnea.gov.ar)
    !
    ! date      : 3/26/97      version: 0.1
    ! revision  : 6/25/02      version: 1.0
    ! revision  : 7/25/02      version: 1.0.1  Patch eq. (23) and (24)
    !                                          check parameters included.
    !
    ! purpose   :  Computes the F1 function following the ideas in 
    !              the paper CPC 138 (2001) 29.
    !
    ! input     :    a  -> complex parameter of Appell's F1
    !                b1 -> complex parameter of Appell's F1
    !                b2 -> complex parameter of Appell's F1
    !                c  -> complex parameter of Appell's F1
    !                x  -> real variable
    !                y  -> real variable
    !
    !
    ! output    :    f1  -> F1 Appell's hypergeometric function
    !
    ! modules   :
    !
    !
    ! common    :
    !
    !
    ! notes     :  most of the code is f77, few things come from f90.
    !
    !-------------------------------------------------------------------

      complex*16 function f1(a,b1,b2,c,x,y)
      implicit none 
      logical test,isthere,ispossible
      logical debug
      integer*4 flag,isneg,pflag
      !integer get_transformation
      integer get_transformation2  
      integer check_parms
      real*8 x,y
      real*8 zero
      parameter(zero = 1d-18)
      real*8 tmax1,tmax2,d1
      complex*16 a,b1,b2,c,cgamma,f21,s,g2,cx,cy,cgammar
      complex*16 f1conv,sa,sb,sc,f1aux,g2aux,hypgeof1,f1bnl  
      complex*16 gc,g1,g4,g3,g5,g6,g7
      complex*16 g2gam
      character*30 outformat

      outformat = '(a20,f17.8,a4,f17.8)'
  
      debug = .false.

      test = .true.
      cx=x*(1.0,0.0)
      cy=y*(1.0,0.0)
 
      flag=-1
      ispossible=.false.
      !
      !   Some particular cases
      !
      if (dabs(x-1d0).lt.zero .and. dabs(y-1d0).lt.zero) then 
        f1 = cgamma(c)*cgamma(c-a-b1-b2)*(cgammar(c-a)*cgammar(c-b1-b2))
        ispossible=.true.
      else if(dabs(x-y).lt.zero) then
        f1 = f21(a,b1+b2,c,cx)
        ispossible=.true.
      else if(cdabs(a-c).lt.zero) then
        f1 = (1-x)**(-b1)*(1-y)**(-b2)
        ispossible=.true.
      else if(dabs(y-1d0).lt.zero) then
        f1 = cgamma(c)*cgamma(c-a-b2)*f21(a,b1,c-b2,cx)*(cgammar(c-a)
     .      *cgammar(c-b2))
        ispossible=.true.
      else if(dabs(x-1d0).lt.zero) then
        f1 = cgamma(c)*cgamma(c-a-b1)*f21(a,b2,c-b1,cy)*(cgammar(c-a)
     .      *cgammar(c-b1))
        ispossible=.true.
      else if(x.eq.0d0.or.dabs(x).lt.zero) then
        f1 = f21(a,b2,c,cy)
        ispossible=.true.
      else if(y.eq.0d0.or.dabs(y).lt.zero) then
        f1 = f21(a,b1,c,cx)
        ispossible=.true.
      else if (cdabs(b1).lt.zero) then
        f1 = f21(a,b2,c,cy)
        ispossible=.true.
      else if (cdabs(b2).lt.zero) then
        f1 = f21(a,b1,c,cx)
        ispossible=.true.
      elseif (isneg(b1).eq.-1.or.isneg(b2).eq.-1.or.isneg(a).eq.-1) then
        ispossible=.true.
        f1 = f1conv(a,b1,b2,c,x,y)
      endif  
      if(ispossible) flag = 5
    
      if(debug) then
        inquire(file='flag.dat',exist=isthere)
        if(isthere) then
            open(unit=23,file='flag.dat',position='append')
        else
            open(unit=23,file='flag.dat',status='new')
        end if
        if(ispossible) then
            write(23,'(2f7.3,i4)') x,y,flag
        end if
        close(23)
      endif
      if(ispossible) then
        if(debug) write(*,*) "Simple transformations: "
        if(debug) call writef1(6,a,b1,b2,c,cx,cy,f1)
        return
      end if
      !
      !   Polynomial cases
      !
      if(isneg(a).le.0.or.isneg(b1).le.0.or.isneg(b2).le.0.or.isneg(c-a)
     .   .le.0) then
        ispossible = .true.
        f1 = f1bnl(a,b1,b2,c,x,y)
      end if
      !if(isneg(c-b1-b2).le.0) then
      !    ispossible = .true.    
      !    f1 = f1bnl(a,c-b1-b2,b2,c,x/(-1 + x),(x-y)/(-1 + x))/(1-x)**a            
      !end if
      if(ispossible) then
        if(debug) write(*,*) "Polynomial transformations: "
        if(debug) call writef1(6,a,b1,b2,c,cx,cy,f1)
        return
      end if    
    
      !
      !
      ! Check regions
      !  
      flag = get_transformation2(x,y)
      !
      ! Check parameters
      !
      pflag = check_parms(a,b1,b2,c,flag)

      if(flag.eq.0) then
        if (debug) write(*,outformat) 'B&Ch Series:     x= ',x,' y= ',y
        s= f1conv(a,b1,b2,c,x,y)
      else if(flag.eq.1) then  
        if (debug) write(*,outformat) 'ODE Integr.:     x= ',x,' y= ',y
        s= hypgeof1(a,b1,b2,c,cx,cy)
      else if(flag.eq.15) then
      !
      !   Eq. (15)      
      !
      if (debug) write(*,outformat) '                 x= ',x,' y= ',y
      if (debug) write(*,outformat) 'Transf. (15):    u= ',x/(x-1),' w= 
     .',y/(y-1)
        s= f1conv(-a+c,b1,b2,c,x/(-1+x),y/(-1+y))/((1-x)**b1*(1-y)**b2)
      else if(flag.eq.16) then
      !
      !   Eq. (16)      
      !
      if (debug) write(*,outformat) '                 x= ',x,' y= ',y
      if (debug) write(*,outformat) 'Transf. (16):    u= ',x/(x-1),' w= 
     .',(x-y)/(x-1)
        s= f1conv(a,-b1-b2 + c,b2,c,x/(-1 + x),(x-y)/(-1 + x))/(1-x)**a
      else if(flag.eq.17) then
      !
      !   Eq. (17)      
      !
      if (debug) write(*,outformat) '                 x= ',x,' y= ',y    
      if (debug) write(*,outformat) 'Transf. (17):    u= ',y/(y-1),' w= 
     .',(y-x)/(y-1)
      s= f1conv(a,b1,-b1-b2 + c,c,(-x + y)/(-1 + y),y/(-1 + y))/(1-y)**a
      else if(flag.eq.21) then
      !
      !   Eq. (21)       Transformation in (1,1)x 
      !
      if (debug) write(*,outformat) '                 x= ',x,' y= ',y    
      if (debug) write(*,outformat) 'Transf. (21):    u= ',1-x,' w= ',1-
     .  y    
      s=cgamma(c)*cgamma(c-a-b1-b2)*cgammar(c-a)*cgammar(c-b1-b2)                           
     .   *f1conv(a,b1,b2,1+a+b1+b2-c,1-x,1-y)+                                               
     .   cgamma(c)*cgamma(a+b2-c)*(cgammar(a)*cgammar(b2))                                   
     .   *(1-x)**(-b1)*(1-y)**(c-a-b2)*f1conv(c-a,b1,c-b1-b2,c-a-b2+1,
     .   (1-y)/(1-x),1-y)+       
     .   cgamma(c)*cgamma(c-a-b2)*cgamma(a+b1+b2-c)*(cgammar(a)*
     .   cgammar(b1)*cgammar(c-a))   
     .   *(1-x)**(c-a-b1-b2)*g2(c-b1-b2,b2,a+b1+b2-c,c-a-b2,x-1,(1-y)/
     .   (x-1))
      else if(flag.eq.22)  then
      !
      !   Eq. (22)      Transformation in (1,1)y
      !
      if (debug) write(*,outformat) '                 x= ',x,' y= ',y    
      if (debug) write(*,outformat) 'Transf. (22):    u= ',1-x,' w= ',1-
     .  y    
      sa = cgamma(c)*cgamma(c-a-b1-b2)*cgammar(c-a)*cgammar(c-b1-b2)*
     .  f1conv(a,b1,b2,1+a+b1+b2-c, 1-x,1-y)
      sb = cgamma(c)*cgamma(a+b1-c)*(1-y)**(-b2)*(1-x)**(c-a-b1)                    
     .      *f1conv(c-a,c-b2-b1,b2,c-a-b1+1, 1-x,(1-x)/(1-y))*
     .      cgammar(a)*cgammar(b1)
      sc = cgamma(c)*cgamma(c-a-b1)*cgamma(a+b1+b2-c)*(1-y)**(c-a-b1-b2)             
     .      *g2(b1,c-b1-b2,c-a-b1,a+b1+b2-c,(1-x)/(y-1),y-1)*cgammar(a)
     .      *cgammar(b2)*cgammar(c-a)
      s = sa+sb+sc
      else if(flag.eq.23) then 
      !
      !   Eq. (23)      Transformation in (0,Inf.)
      !
      if (debug) write(*,outformat) '                 x= ',x,' y= ',y    
      if (debug) write(*,outformat) 'Transf. (23):    u= ',x/y,' w= ',
     .  1/y
      if(pflag.eq.0) then 
          ! Normal case
          s= cgamma(c)*cgamma(b2-a)*(cgammar(b2)*cgammar(c-a))        
     .        *(-y)**(-a)*f1conv(a,b1,1+a-c,a-b2+1,x/y,1/y)+           
     .        (cgamma(c)*cgamma(a-b2)*(cgammar(a)*cgammar(c-b2)))      
     .        *(-y)**(-b2)*g2(b1,b2,1+b2-c,a-b2,-x,-1/y)
      else
          ! c-b2=neg. int.
          s= cgamma(c)*cgamma(b2-a)*(cgammar(b2)*cgammar(c-a))        
     .        *(-y)**(-a)*f1conv(a,b1,1+a-c,a-b2+1,x/y,1/y)+           
     .        (cgamma(c)*cgamma(a-b2)*cgammar(a))      
     .        *(-y)**(-b2)*g2gam(b1,b2,1+b2-c,a-b2,-x,-1/y)
      end if      
      else if(flag.eq.24) then
      !
      !   Eq. (24)      Transformation in (Inf.0)
      !
      if (debug) write(*,outformat) '                 x= ',x,' y= ',y
      if (debug) write(*,outformat) 'Transf. (24):    u= ',y/x,' w= ',
     .  1/x
      if(pflag.eq.0) then
          ! Normal case
          s=cgamma(c)*cgamma(b1-a)*cgammar(b1)*cgammar(c-a)     
     .     *(-x)**(-a)*f1conv(a,1+a-c,b2,a-b1+1,1/x,y/x)+          
     .     cgamma(c)*cgamma(a-b1)*(cgammar(a)* cgammar(c-b1))     
     .     *(-x)**(-b1)*g2(b1,b2,a-b1,1+b1-c,-1/x,-y)      
      else   
          ! c-b1=neg. int.
          s=cgamma(c)*cgamma(b1-a)*cgammar(b1)*cgammar(c-a)     
     .     *(-x)**(-a)*f1conv(a,1+a-c,b2,a-b1+1,1/x,y/x)+          
     .     cgamma(c)*cgamma(a-b1)*cgammar(a)                    
     .     *(-x)**(-b1)*g2gam(b1,b2,a-b1,1+b1-c,-1/x,-y)      
      end if
      else if(flag.eq.25) then
      !
      !   Eq. (25)      Transformation in (1,Inf.)
      !
      if (debug) write(*,outformat) '                 x= ',x,' y= ',y
      if (debug) write(*,outformat) 'Transf. (25):    u= ',(1-x)/(1-y),
     .' w= ',1/(1-y)    
      s=cgamma(c)*cgamma(b2-a)*(cgammar(c-a)*cgammar(b2))                                  
     .   *(1-y)**(-a)*f1conv(a,b1,c-b1-b2,1+a-b2,(1-x)/(1-y),1/(1-y))+                      
     .   cgamma(c)*cgamma(a+b1-c)*(cgammar(a)*cgammar(b1))                                  
     .   *(1-x)**(c-a-b1)*(1-y)**(-b2)*f1conv(c-a,b2,c-b2-b1,c-a-b1+ 1,
     .   (1-x)/(1-y),1-x)+    
     .   cgamma(c)*cgamma(a-b2)*cgamma(c-a-b1)*(cgammar(a)*
     .   cgammar(c-b1-b2)*cgammar(c-a))  
     .   *(1-y)**(-b2)*g2(b1,b2,c-a-b1,a-b2,x-1,1/(y-1))
      else if(flag.eq.26) then
      !
      !   Eq. (26)      Transformation in (Inf.,1)
      !
      if (debug) write(*,outformat) '                 x= ',x,' y= ',y
      if (debug) write(*,outformat) 'Transf. (26):    u= ',(1-y)/(1-x),
     .  ' w= ',1/(1-x)        
      s=cgamma(c)*cgamma(b1-a)*(cgammar(c-a)*cgammar(b1))                                    
     .   *(1-x)**(-a)*f1conv(a,c-b1-b2,b2,1+a-b1,1/(1-x),(1-y)/(1-x))+                        
     .   cgamma(c)*cgamma(a+b2-c)*(cgammar(a)*cgammar(b2))                                    
     .   *(1-y)**(c-a-b2)*(1-x)**(-b1)*f1conv(c-a,b1,c-b1-b2,c-a-b2+1,
     .   (1-y)/(1-x),1-y)+      
     .   cgamma(c)*cgamma(a-b1)*cgamma(c-a-b2)*(cgammar(a)*
     .   cgammar(c-b1-b2)*cgammar(c-a))    
     .   *(1-x)**(-b1)*g2(b1,b2,a-b1,c-a-b2,1/(x-1),y-1)
      else if(flag.eq.27) then
      !
      !   Eq. (27)      Transformation in (Inf.Inf.)x
      !
      if (debug) write(*,outformat) '                 x= ',x,' y= ',y    
      if (debug) write(*,outformat) 'Transf. (27):    u= ',x/y,' w= ',
     .  1/x 
      s=cgamma(c)*cgamma(b2-a)*(cgammar(c-a)*cgammar(b2))                                
     .   *(-y)**(-a)*f1conv(a,b1,1+a-c,1+a-b2,x/y,1/y)+                                  
     .   cgamma(c)*cgamma(a-b1-b2)*(cgammar(a)*cgammar(c-b1-b2))                          
     .   *(-x)**(-b1)*(-y)**(-b2)*
     .   f1conv(1+b1+b2-c,b1,b2,1+b1+b2-a,1/x,1/y)+              
     .   cgamma(c)*cgamma(a-b2)*cgamma(b1+b2-a)*(cgammar(a)*cgammar(b1)*
     .   cgammar(c-a))    
     .   *(-x)**(b2-a)*(-y)**(-b2)*g2(1+a-c,b2,b1+b2-a,a-b2,-1/x,-x/y)
      else if(flag.eq.28) then
      !
      !   Eq. (28)      Transformation in (Inf.Inf.)y
      !
      if (debug) write(*,outformat) '                 x= ',x,' y= ',y    
      if (debug) write(*,outformat) 'Transf. (28):    u= ',1/x,' w= ',
     .  1/y 
      s= cgamma(c)*cgamma(b1-a)*(cgammar(c-a)*cgammar(b1))                            
     .   *(-x)**(-a)*f1conv(a,1+a-c,b2,1+a-b1,1/x,y/x)+                                
     .   cgamma(c)*cgamma(a-b1-b2)*(cgammar(a)*cgammar(c-b1-b2))                        
     .   *(-x)**(-b1)*(-y)**(-b2)*
     .   f1conv(1+b1+b2-c,b1,b2,1+b1+b2-a,1/x,1/y)+            
     .   cgamma(c)*cgamma(a-b1)*cgamma(b1+b2-a)*(cgammar(a)*cgammar(b2)*
     .   cgammar(c-a))  
     .   *(-x)**(-b1)*(-y)**(b1-a)*g2(b1,1+a-c,a-b1,b1+b2-a,-y/x,-1/y)

      else if(flag.eq.29) then
      !
      !   Eq. (29)      Transformation in (Inf.Inf.)xx
      !
      if (debug) write(*,outformat) '                 x= ',x,' y= ',y    
      if (debug) write(*,outformat) 'Transf. (29):    u= ',
     .  (x-y)/(y*(x-1)),' w= ',1/y 
      f1aux = f1conv(c-a,b1,1-a,1+b1+b2-a,(x-y)/(y*(x-1)),1/y)
      g2aux = g2(b1,c-b1-b2,1-b1-b2,b1+b2-a,(x-y)/(1-x),-1/y)
      s= cgamma(c)*cgamma(a-b1-b2)*(cgammar(a)*cgammar(c-b1-b2))       
     .   *(-y)**(a-c)*(1-x)**(-b1)*(1-y)**(c-a-b2)*f1conv(c-a,b1,1-a,
     .   1+b1+b2-a,(x-y)/(y*(x-1)),1/y)+    
     .   cgamma(c)*cgamma(b1+b2-a)*(cgammar(c-a)*cgammar(b1+b2))          
     .   *(-y)**(b1+b2-c)*(1-x)**(-b1)*(1-y)**(c-a-b2)*g2(b1,c-b1-b2,
     .   1-b1-b2,b1+b2-a,(x-y)/(1-x),-1/y)  
      else if(flag.eq.30) then
      !
      !   Eq. (30)      Transformation in (Inf.Inf.)xx
      !
      if (debug) write(*,outformat) '                 x= ',x,' y= ',y
      if (debug) write(*,outformat) 'Transf. (29):    u= ',(y-x)/
     .  (x*(y-1)),' w= ',1/x
      s=cgamma(c)*cgamma(a-b1-b2)*(cgammar(a)*cgammar(c-b1-b2))                                        
     .   *(-x)**(a-c)*(1-y)**(-b2)*(1-x)**(c-a-b1)*f1conv(c-a,1-a,b2,
     .  1+b1+b2-a,1/x,(y-x)/(x*(y-1)))+    
     .   cgamma(c)*cgamma(b1+b2-a)*(cgammar(c-a)*cgammar(b1+b2))                                        
     .   *(-x)**(b1+b2-c)*(1-y)**(-b2)*(1-x)**(c-a-b1)*g2(c-b1-b2,b2,
     .  b1+b2-a,1-b1-b2,-1/x,(y-x)/(1-y))
      end if 
      f1 = s

      if(flag.ge.0) then
        if(debug) call writef1(6,a,b1,b2,c,cx,cy,f1)
      else
        write(*,*) 'Not Possible'
        call writef1(6,a,b1,b2,c,cx,cy,f1)
      end if




      return
      end 

      