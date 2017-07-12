      !-----------------------------------------------------------------
      ! function  : g2gam (complex*16)
      !
      ! package   : F1
      !
      ! Language  : Fortran 90
      !
      ! author    : F. Colavecchia (flavioc@lanl.gov)
      !                            (flavioc@cab.cnea.gov.ar)
      !
      ! date      : 7/25/02      version: 1.0
      !
      ! purpose   :  Computes the G2 function times 1/Gamma[b1].
      !
      !
      ! input     :    a1  -> complex parameter of Appell's F2
      !                a2 -> complex parameter of Appell's F2
      !                b1 -> complex parameter of Appell's F2
      !                b2 -> complex parameter of Appell's F2
      !                x  -> real variable
      !                y  -> real variable
      !
      ! output    :    g2  ->G2 function
      !
      !-----------------------------------------------------------------
      complex*16 function g2gam(a1,a2,b1,b2,x,y)

      real*8 x,y  
      complex*16 a1,a2,b1,b2,f2gam
      complex*16 cx,cy
      logical debug

      debug = .false.
      cx = x*(1d0,0d0)
      cy = y*(1d0,0d0)
      if (debug) write(*,*) " Computing g2gam"
      if (debug) write(*,*) " x =",x
      if (debug) write(*,*) " y =",y
      if (debug) write(*,*) " a1=",a1
      if (debug) write(*,*) " a2=",a2
      if (debug) write(*,*) " b1=",b1
      if (debug) write(*,*) " b2=",b2

      g2gam=(1+cx)**(-a1)*(1+cy)**(-a2)*
     .  f2gam(1-b1-b2,a1,a2,1-b1,1-b2,cx/(cx+1),cy/(cy+1))
      if (debug) write(*,*) g2gam
      if (debug) write(*,*) " End computing g2gam"
      return
      end

      !-----------------------------------------------------------------
      ! function  : f2gam (complex*16)
      !
      ! package   : F1
      !
      ! Language  : Fortran 90
      !
      ! author    : F. Colavecchia (flavioc@lanl.gov)
      !                            (flavioc@cab.cnea.gov.ar)
      !
      ! date      : 7/25/02      version: 0.1
      !
      ! purpose   :  Computes the  F2 series times 1/Gamma[c1]
      !              in the convergence region.
      !
      !
      ! input     :    a  -> complex parameter of Appell's F2
      !                b1 -> complex parameter of Appell's F2
      !                b2 -> complex parameter of Appell's F2
      !                c1 -> complex parameter of Appell's F2
      !                c2 -> complex parameter of Appell's F2
      !                x  -> real variable
      !                y  -> real variable
      !
      ! output    :    f2  -> nth sum of F2
      !
      !-----------------------------------------------------------------
      complex*16 function f2gam(a,b1,b2,c1,c2,cx,cy)
      implicit none 
      logical debug, ispossible
      complex*16 a,b1,b2,c1,c2,f21,s,f2s_gamma,c0,cx,cy
      real*8 tmax1,tmax2
      integer flag
      
      ispossible =.false.
      debug = .false.

      c0 = dcmplx(0.,0.)
      if (debug) write(*,*) " Computing f2gam"
      if (debug) write(*,*) " cx =",cx
      if (debug) write(*,*) " cy =",cy
      if (debug) write(*,*) " a  =",a
      if (debug) write(*,*) " b1 =",b1
      if (debug) write(*,*) " b2 =",b2
      if (debug) write(*,*) " c1 =",c1
      if (debug) write(*,*) " c2 =",c2

      !
      !   Simple cases for 1/Gamma[c2] F2 need
      !   the computation of 1/Gamma[c] 2F1[a,b,c,z]
      !   which  is not yet implemented.
      !
      !
      !if(b1.eq.c1) then
      !    if (debug) write(*,*) "b1=c1 " 
      !    f2 = (1-cx)**(-a)*f21(a,b2,c2,cy/(1-cx))
      !    return
      !else if(b2.eq.c2) then  
      !    if (debug) write(*,*) "b2=c2 " 
      !    f2 =  (1-cy)**(-a)*f21(a,b1,c1,cx/(1-cy))
      !    return
      !else if(b1.eq.0) then
      !    if (debug) write(*,*) "b1=0  " 
      !    f2 = f21(a,b2,c2,cy)
      !    return
      !else if(b2.eq.0) then 
      !    if (debug) write(*,*) "b2=0  " 
      !    f2 =  f21(a,b1,c1,cx)
      !    return
      !else if(cx.eq.0) then 
      !    f2 =  f21(a,b2,c2,cy)
      !    return
      !else if(cy.eq.0) then
      !    f2 = f21(a,b1,c1,cx)
      !    return
      !end if 

      flag=  0 
      tmax1=  1.0 
      tmax2=  cdsqrt(cx**2+cy**2) 
      if (debug) write(*,*) flag," ", tmax1,"   ",tmax2
      if(tmax2.lt.tmax1) then
          flag=  1 
          ispossible = .true. 
          tmax1=tmax2
      end if
      if (debug) write(*,*) flag," ", tmax1,"   ",tmax2

      tmax2=  cdsqrt((cx/(1-cy))**2+(cy/(cy-1))**2) 
      if(tmax2.lt.tmax1) then
          flag=  22
          ispossible = .true. 
          tmax1=tmax2
      end if
      if (debug) write(*,*) flag," ", tmax1,"   ",tmax2 

      tmax2=  cdsqrt((cx/(cx-1))**2+(cy/(1-cx))**2) 
      if(tmax2.lt.tmax1) then
          flag=  21 
          ispossible = .true. 
          tmax1=tmax2
      end if

      if (debug) write(*,*) flag," ", tmax1,"   ",tmax2 
      tmax2=  cdsqrt((cx/(cx+cy-1))**2+(cy/(cx+cy-1))**2) 
      if(tmax2.lt.tmax1) then
          flag=  3
          ispossible = .true. 
          tmax1=tmax2
      end if

      if (debug) write(*,*) flag," ", tmax1,"   ",tmax2
      if(flag.eq.1)  then
           if (debug) write(*,*) "Series in f21(cx)     :",cx," ",cy
          s=f2s_gamma(a,b1,b2,c1,c2,cx,cy)
      else if(flag.eq.21)  then
          if (debug) write(*,*) "Transformation 2x:",cx/(cx-1)," ",
     .      cy/(1-cx)
          s=(1-cx)**(-a)*f2s_gamma(a,c1-b1,b2,c1,c2,cx/(cx-1),cy/(1-cx))
      else if(flag.eq.22)  then  
          if (debug) write(*,*) "Transformation 2y:",cx/(1-cy)," ",
     .      cx/(1-cy) 
          s=(1-cy)**(-a)*f2s_gamma(a,b1,c2-b2,c1,c2,cx/(1-cy),cy/(cy-1))
      else if(flag.eq.3)  then
          if (debug) write(*,*) "Transformation 3: ",cx/(cx+cy-1)," ",
     .      cy/(cx+cy-1) 
          s=(1-cx-cy)**(-a)*f2s_gamma(a,c1-b1,c2-b2,c1,c2,cx/(cx+cy-1),
     .      cy/(cx+cy-1))   
      else  
          if (debug) write(*,*) "Not Possible"
          f2gam = c0
          return
      end if
      f2gam = s
      return
      end
      !-----------------------------------------------------------------
      ! function  : f2s_gamma (complex*16)
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
      ! purpose   :  Computes the  F2 series times 1/Gamma[c2]
      !              in the convergence region.
      !
      !
      ! input     :    a  -> complex parameter of Appell's F2
      !                b1 -> complex parameter of Appell's F2
      !                b2 -> complex parameter of Appell's F2
      !                c1 -> complex parameter of Appell's F2
      !                c2 -> complex parameter of Appell's F2
      !                cx  -> complex variable
      !                cy  -> complex variable
      !
      ! output    :    f2s_gamma  -> sum of F2 series
      !
      !-----------------------------------------------------------------
      complex*16 function f2s_gamma(a,b1,b2,c1,c2,cx,cy)
      implicit none 
      complex*16 a,b1,b2,c1,c2,cx,cy,f2sx_gamma
      integer isneg

      if(isneg(c1).le.0) then 
          f2s_gamma = f2sx_gamma(a,b2,b1,c2,c1,cy,cx)
      else if(isneg(c2).le.0) then
          f2s_gamma = f2sx_gamma(a,b1,b2,c1,c2,cx,cy)
      end if
      
      !if(cdabs(cy).lt.cdabs(cx)) then
      !    f2s_gamma = f2sx_gamma(a,b2,b1,c2,c1,cy,cx)
      !else
      !    f2s_gamma = f2sx_gamma(a,b1,b2,c1,c2,cx,cy)
      !end if
      return
      end
      !-----------------------------------------------------------------
      ! function  : f2sx_gamma (complex*16)
      !
      ! package   : F1
      !
      ! Language  : Fortran 90
      !
      ! author    : F. Colavecchia (flavioc@lanl.gov)
      !                            (flavioc@cab.cnea.gov.ar)
      !
      ! date      : 7/25/02      version: 1.0
      !
      ! purpose   :  Computes the  F2 series times 1/Gamma[c2]
      !              in the convergence region, single index,
      !              single Gauss function.
      !
      !
      ! input     :    a  -> complex parameter of Appell's F2
      !                b1 -> complex parameter of Appell's F2
      !                b2 -> complex parameter of Appell's F2
      !                c1 -> complex parameter of Appell's F2
      !                c2 -> complex parameter of Appell's F2
      !                x  -> real variable
      !                y  -> real variable
      !
      ! output    :    f2sx  -> converged sum F2
      !
      ! notes     :  most of the code is f77, few things come from f90.
      !
      !-----------------------------------------------------------------
      complex*16 function f2sx_gamma(a,b1,b2,c1,c2,cx,cy)
      implicit none 
      logical debug
      integer*4 i,isneg
      complex*16 a,b1,b2,c1,c2,cx,cy
      complex*16 coef,f21,suma,tmp,c,cgammar
      real*8 rtest


      debug = .false.
      if (debug) write(*,*) " Computing f2sx_gamma"
      if (debug) write(*,*) " cx =",cx
      if (debug) write(*,*) " cy =",cy
      if (debug) write(*,*) " a  =",a
      if (debug) write(*,*) " b1 =",b1
      if (debug) write(*,*) " b2 =",b2
      if (debug) write(*,*) " c1 =",c1
      if (debug) write(*,*) " c2 =",c2

      if(isneg(c1).le.0) then
          stop 'Error in f2sx_gamma, c1 is zero or negative integer'
      end if
      
      suma = f21(a,b1,c1,cx)*cgammar(c2) 
      tmp = suma
      coef = (1d0,0d0)
      c = (0d0,0d0)
      i=0 
      rtest = tmp/suma
      do while(i.lt.300 .and. (rtest.gt.1e-5))
          i=i+1 
          coef  = coef*(a+i-1)*(b2+i-1)/(i*1d0)
          tmp   = coef*cgammar(c2+i)*f21(a+i,b1,c1,cx)*cy**(i)
          suma  = suma+tmp
          rtest = cdabs(tmp/suma)
      end do
      if(i.gt.200) then 
      !  if (debug) write(*,*) 'Max. number of terms obtained'
      end if
      f2sx_gamma = suma
      !  if (debug) write(*,*) i
      return
      end
