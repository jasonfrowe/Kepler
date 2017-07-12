      !-----------------------------------------------------------------
      ! function  : g2 (complex*16)
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
      ! purpose   :  Computes the G2 function, Eq. (34) of paper CPC 138
      ! (2001) 29.
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
      complex*16 function g2(a1,a2,b1,b2,x,y)

          real*8 x,y  
          complex*16 a1,a2,b1,b2,f2
          complex*16 cx,cy
          logical debug

          debug = .false.
          cx = x*(1d0,0d0)
          cy = y*(1d0,0d0)
          if (debug) write(*,*) " Computing g2"
          if (debug) write(*,*) " x=",x
          if (debug) write(*,*) " y=",y

          if (debug) write(*,*) a1," ",a2," ",b1," ",b2
          g2=(1+cx)**(-a1)*(1+cy)**(-a2)*f2(1-b1-b2,a1,a2,1-b1,1-b2,
     .      cx/(cx+1),cy/(cy+1))
          if (debug) write(*,*) g2
          if (debug) write(*,*) " End computing g2"
          return
          end

      !-----------------------------------------------------------------
      ! function  : f2 (complex*16)
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
      ! purpose   :  Computes the  F2 series 
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
          complex*16 function f2(a,b1,b2,c1,c2,cx,cy)
          implicit none 
      logical debug, ispossible
          complex*16 a,b1,b2,c1,c2,f21,s,f2s,c0,cx,cy
          real*8 tmax1,tmax2
          integer flag
    
      ispossible =.false.
      debug = .false.

          c0 = dcmplx(0.,0.)
          if(b1.eq.c1) then
                  if (debug) write(*,*) "b1=c1 " 
                  f2 = (1-cx)**(-a)*f21(a,b2,c2,cy/(1-cx))
                  return
          else if(b2.eq.c2) then  
                    if (debug) write(*,*) "b2=c2 " 
                    f2 =  (1-cy)**(-a)*f21(a,b1,c1,cx/(1-cy))
                    return
          else if(b1.eq.0) then
                  if (debug) write(*,*) "b1=0  " 
                  f2 = f21(a,b2,c2,cy)
                    return
          else if(b2.eq.0) then 
                  if (debug) write(*,*) "b2=0  " 
                  f2 =  f21(a,b1,c1,cx)
                    return
          else if(cx.eq.0) then 
                  f2 =  f21(a,b2,c2,cy)
                  return
          else if(cy.eq.0) then
              f2 = f21(a,b1,c1,cx)
          return
          end if 

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
      if(flag.eq.1)       then
                  if (debug) write(*,*) "Series in f21(cx)     :",cx,
     .              " ",cy
                  s=f2s(a,b1,b2,c1,c2,cx,cy)
          else if(flag.eq.21)   then
                  if (debug) write(*,*) "Transformation 2x:",cx/(cx-1),
     .              " ",cy/(1-cx)
                  s=(1-cx)**(-a)*f2s(a,c1-b1,b2,c1,c2,cx/(cx-1),
     .              cy/(1-cx))
          else if(flag.eq.22)   then  
                  if (debug) write(*,*) "Transformation 2y:",cx/(1-cy),
     .              " ",cx/(1-cy) 
                  s=(1-cy)**(-a)*f2s(a,b1,c2-b2,c1,c2,cx/(1-cy),
     .              cy/(cy-1))
          else if(flag.eq.3)    then
                  if (debug) write(*,*) "Transformation 3:",
     .              cx/(cx+cy-1)," ",cy/(cx+cy-1) 
                  s=(1-cx-cy)**(-a)*f2s(a,c1-b1,c2-b2,c1,c2,cx/
     .              (cx+cy-1),cy/(cx+cy-1))   
          else  
                  if (debug) write(*,*) "Not Possible"
                  f2 = c0
                  return
          end if
          f2 = s
          return
          end
      !-----------------------------------------------------------------
      ! function  : f2s (complex*16)
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
      ! purpose   :  Computes the  F2 series 
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
      ! output    :    f2s  -> sum of F2 series
      !
      !-----------------------------------------------------------------
          complex*16 function f2s(a,b1,b2,c1,c2,cx,cy)
          implicit none 
          complex*16 a,b1,b2,c1,c2,cx,cy,f2sx

          if(cdabs(cy).lt.cdabs(cx)) then
                  f2s = f2sx(a,b2,b1,c2,c1,cy,cx)
          else
                    f2s = f2sx(a,b1,b2,c1,c2,cx,cy)
          end if
          return
          end
      !-----------------------------------------------------------------
      ! function  : f2sx (complex*16)
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
      ! purpose   :  Computes the  F2 series 
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
          complex*16 function f2sx(a,b1,b2,c1,c2,cx,cy)
          implicit none 
          integer*4 i
          complex*16 a,b1,b2,c1,c2,cx,cy
          complex*16 coef,f21,suma,tmp,c
          real*8 rtest

          ! if (debug) write(*,*) a," ",b1," ",b2," ",c1," ",c2
          !if (debug) write(*,*) "Calculando f2 con serie en f21 :",x,"
          !",y

          suma = f21(a,b1,c1,cx) 
          tmp = suma
          coef = (1d0,0d0)
          c = (0d0,0d0)
          i=0 
          rtest = tmp/suma
          do while(i.lt.300 .and. (rtest.gt.1e-5))
                    i=i+1 
                    coef = coef*(a+i-1)*(b2+i-1)/((c2+i-1)*i)
                    tmp  = coef*f21(a+i,b1,c1,cx)*cy**(i)
                  suma=suma+tmp
                  rtest = cdabs(tmp/suma)
          end do
          if(i.gt.200) then 
          !     if (debug) write(*,*) 'Max. number of terms obtained'
          end if
          f2sx = suma
  !     if (debug) write(*,*) i
          return
          end


      !-----------------------------------------------------------------        
          !
          ! Next functions are included for debugging purposes.
          !
      !-----------------------------------------------------------------
      ! function  : coef2 (complex*16)
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
      ! purpose   :  Computes the n-th coefficient of the F2 series 
      !              of Burchnall and Chaundy, Eq. (36) of paper CPC 138
      !              (2001) 29.
      !
      !
      ! input     :    a  -> complex parameter of Appell's F2
      !                b1 -> complex parameter of Appell's F2
      !                b2 -> complex parameter of Appell's F2
      !                c1 -> complex parameter of Appell's F2
      !                c2 -> complex parameter of Appell's F2
      !                n  -> order
      !
      ! output    :    coef2  -> coefficient
      !
      ! notes     :  most of the code is f77, few things come from f90.
      !
      !-----------------------------------------------------------------
          complex*16 function coef2(a,b1,b2,c1,c2,n)
          implicit none 
          complex*16 a,b1,b2,c1,c2,pochhammer
          real*8 fact
          integer*4 n
             coef2 = pochhammer(a,n)*pochhammer(b1,n)*pochhammer(b2,n)/
     .          (pochhammer(c1,n)*pochhammer(c2,n)*fact(n))
             return
          end
      !-----------------------------------------------------------------
      ! function  : f2term (complex*16)
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
      ! purpose   :  Computes the n-th term of the F2 series 
      !              of Burchnall and Chaundy, Eq. (36) of paper CPC 138
      !              (2001) 29.
      !
      !
      ! input     :    a  -> complex parameter of Appell's F2
      !                b1 -> complex parameter of Appell's F2
      !                b2 -> complex parameter of Appell's F2
      !                c1 -> complex parameter of Appell's F2
      !                c2 -> complex parameter of Appell's F2
      !                x  -> real variable
      !                y  -> real variable
      !                n  -> order
      !
      ! output    :    coef2  -> nth term of F2
      !
      ! notes     :  most of the code is f77, few things come from f90.
      !
      !-----------------------------------------------------------------
          complex*16 function f2term(a,b1,b2,c1,c2,x,y,n)
          implicit none 
          complex*16 a,b1,b2,c1,c2,coef2,f21
          real*8 x,y 
          complex*16 cx,cy
          integer*4 n

          cx = x*(1d0,0d0)
          cy = y*(1d0,0d0)
          f2term = coef2(a,b1,b2,c1,c2,n)*(x*y)**n*
     .      f21(a+n,b1+n,c1+n,cx)*f21(a+n,b2+n,c2+n,cy)
          return
          end
      !-----------------------------------------------------------------
      ! function  : f2n (complex*16)
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
      ! purpose   :  Computes the n-th sum of the F2 series 
      !              of Burchnall and Chaundy, Eq. (36) of paper CPC 138
      !              (2001) 29.
      !
      !
      ! input     :    a  -> complex parameter of Appell's F2
      !                b1 -> complex parameter of Appell's F2
      !                b2 -> complex parameter of Appell's F2
      !                c1 -> complex parameter of Appell's F2
      !                c2 -> complex parameter of Appell's F2
      !                x  -> real variable
      !                y  -> real variable
      !                n  -> order
      !
      ! output    :    f2n  -> nth sum of F2
      !
      ! notes     :  most of the code is f77, few things come from f90.
      !
      !-----------------------------------------------------------------
          complex*16 function f2n(a,b1,b2,c1,c2,x,y,N)
          implicit none 
          complex*16 a,b1,b2,c1,c2,f2term
          real*8 x,y 
          integer*4 i,n
          do i=0,n
                  f2n = f2n + f2term(a,b1,b2,c1,c2,x,y,i)
          end do
          return 
          end
