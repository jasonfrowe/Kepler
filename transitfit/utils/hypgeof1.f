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
      !
      ! purpose   :  Computes the F1 function following the ideas in
      !              the paper CPC 138 (2001) 29, eq. (33).
      !                  Calculates the Appell's hypergeometric function 
      !                          f1(a,b,b',c,x,y) in the convergence
      !                          region
      !                          through the numerical integration of
      !                          ODE 
      !                          equation representing the PDE system of 
      !                          F1.
      !
      ! input     :  a  -> complex parameter of Appell's F1  
      !              b  -> complex parameter of Appell's F1  
      !              bp -> complex parameter of Appell's F1  
      !              c  -> complex parameter of Appell's F1  
      !              u  -> complex variable                  
      !              v  -> complex variable                  
      !
      !
      ! output    :    f1  -> F1 Appell's hypergeometric function
      !
      ! modules   :    uses Felhberg fifth-four Runge-Kutta method for
      !                integration
      !
      !
      !-----------------------------------------------------------------
      function hypgeof1(a,b,bp,c,u,v)
          implicit none
      complex*16 hypgeof1,a,b,bp,c,u,v,f21
      double precision eps,fact,zero,abserr
      parameter (eps=1.d-8)
          parameter (fact=.1)
          parameter (zero=1e-12)                  
          parameter (abserr=1d-40)                    
      integer kmax,nbad,nok,neqn
      integer iwork(5)
      integer iflag
      external bsstep,hypdrvf1
          real*8 au,av,t0,one,su,sv,w(6)
      real*8 start,finish
      real*8 work(50)
      complex*16 z0,dz,aa,bb,bbp,cc,y(3),xx,yy,z,c1
      common /hypgf1/ aa,bb,bbp,cc,xx,yy,z0,dz
      common /path/ kmax
      kmax=0
          one = 1d0
          c1 = (1d0,0d0)
          au = cdabs(u)
          av = cdabs(v)
          if(au.lt.zero.and.av.lt.zero) then
                    hypgeof1 = c1
                    return
          end if
    !   
    !     There still exists a divergence behavior when u=v 
    !     even for |u,v| < 1, so we skip this case
    !
          if(cdabs(u-v).lt.zero) then
                    hypgeof1 = f21(a,b+bp,c,u)
                    return
          endif

          su = dreal(u)/au
          sv = dreal(v)/av
    !
    !     Select z0, starting integration point
    !
          if(su.lt.0d0.and.dreal(v).lt.one.or.sv.lt.0d0.and.
     .      dreal(u).lt.one) then
                t0 = 1d0/(16d0*max(au,av))
      elseif(au.lt.one.and.av.lt.one) then
                    t0 = 1d0/(5d0*max(au,av))
          elseif(au.gt.one.and.av.lt.one) then
                    t0 = (1d0+fact)/au            
          elseif(av.gt.one.and.au.lt.one) then
                    t0 = (1d0+fact)/av            
          endif
          z0 = dcmplx(t0,0d0)
          z  = dcmplx(1.d0,0.d0)
      aa=a
      bb=b
          bbp=bp
      cc=c
          xx=u
          yy=v
      dz=z-z0
    !
    !     Initial values of the function f1
    !
          call hypstartf1(a,b,bp,c,u,v,t0*c1,y)

      neqn = 6
      iflag = 1
      start = 0d0
      finish = 1d0
    !
    !   Proceed with the integration
    !
      call rkf45(hypdrvf1, neqn, y, start,finish, EPS, abserr, iflag,
     .  work, iwork )
      hypgeof1=y(1)

      if(iflag.ne.2) then
        write(*,*) 'iflag =',iflag 
        write(*,*) 'hypgeof1: Problems with the integration'
      end if
      return
      end 
