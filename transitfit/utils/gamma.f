      !#################################################################
      !gamma.for
      !2/3/98
      !_________
      !
      !   Function Gamma and their relatives.
      !
      !fortran 77 source file
      !
      !by : f. colavecchia 
      !
      !
      !data:      
      !   factorial : The factorial of the natural numbers 
      !                           up to 40 as a global variable fct.
      !
      !functions: 
      !   Pochhammer: Calculates the pochhammer symbols 
      !   fact      : Calculates the factorial.
      !   clgam     : Logarithm of the gamma function
      !   clgama    : Logarithm of the gamma function (used by clgam)
      !   coulf     : Coulomb factor
      !
      !subroutines:
      !
      !   chk: 26/02/98
      !                                                                                    
      !#################################################################
      !

      block data fact40
      real*8 fct(0:40)
      common /factorial/ fct
      data fct /1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,
     .  362880.0,3.6288d6,3.99168d7,4.790016d8,6.2270208d9,
     .  8.71782912d10,1.307674368d12,2.0922789888d13,3.55687428096d14,
     .  6.402373705728d15,1.21645100408832d17,2.43290200817664d18,
     .  5.109094217170944d19,1.12400072777760768d21,
     .  2.585201673888497664d22,
     .  6.2044840173323943936d23,1.5511210043330985984d25,
     .  4.03291461126605635584d26,
     .  1.0888869450418352160768d28,3.04888344611713860501504d29,
     .  8.84176199373970195454362d30,2.652528598121910586363085d32,
     .  8.22283865417792281772556d33,2.631308369336935301672180d35,
     .  8.68331761881188649551819d36,2.952327990396041408476186d38,
     .  1.033314796638614492966665d40,3.71993326789901217467999d41,
     .  1.376375309122634504631598d43,5.23022617466601111760007d44,
     .  2.039788208119744335864028d46,8.15915283247897734345611d47 /
      end 

      !-----------------------------------------------------------------
      !  function Pochhammer                                     20/2/97
      !   ______________
      !
      !   by:             f. colavecchia
      !           
      !                   Calculates the PochHammer symbols
      !   input:  
      !                   complex*16 x
      !                   integer*4 n
      !   output: Pochhammer (x,n)
      !
      !   limitations: 
      !
      !-----------------------------------------------------------------
      !
      complex*16 function Pochhammer(x,n)
      implicit none

      complex*16 x
      integer*4 n,i

      Pochhammer = (1d0,0d0)
      if(n.eq.0) then
            return
      end if

      do i=1,n
            Pochhammer = Pochhammer*(x+i-1)
      end do

      return
      end


      !-----------------------------------------------------------------
      !   function fact                                   20/2/97
      !   ______________
      !
      !   by:             f. colavecchia
      !           
      !                   Return n!
      !   input:  
      !                   integer*4 n
      !   output: n!
      !
      !   limitations: n<40
      !
      !-----------------------------------------------------------------
      !
      real*8 function fact(n)
      implicit none

      integer*4 n
      real*8 a(0:40) 
      data a /1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,
     .  362880.0,3.6288d6,3.99168d7,4.790016d8,6.2270208d9,
     .  8.71782912d10,1.307674368d12,2.0922789888d13,3.55687428096d14,
     .  6.402373705728d15,1.21645100408832d17,2.43290200817664d18,
     .  5.109094217170944d19,1.12400072777760768d21,
     .  2.585201673888497664d22,
     .  6.2044840173323943936d23,1.5511210043330985984d25,
     .  4.03291461126605635584d26,
     .  1.0888869450418352160768d28,3.04888344611713860501504d29,
     .  8.84176199373970195454362d30,2.652528598121910586363085d32,
     .  8.22283865417792281772556d33,2.631308369336935301672180d35,
     .  8.68331761881188649551819d36,2.952327990396041408476186d38,
     .  1.033314796638614492966665d40,3.71993326789901217467999d41,
     .  1.376375309122634504631598d43,5.23022617466601111760007d44,
     .  2.039788208119744335864028d46,8.15915283247897734345611d47 /

      if(n.gt.40) then
            write(*,*) 'Too Large Factorial'
            return
      end if

      fact = a(n)
      return
      end


      !   function cgammar
      !   ----------------
      !   Reciprocal of the gamma function
      !   to avoid poles.
      complex*16 function cgammar(z)
      implicit none
      complex*16 z,gamm,cgamma
      real*8 az,err

      err = 1.d-10
      az = cdabs(z)
      if(az.eq.0d0) then
              cgammar = (0d0,0d0)
            return
      end if
      az = dabs(dreal(z)-dnint(dreal(z)))
      if(az.lt.err.and.dabs(dimag(z)).lt.err.and.dreal(z).lt.0d0) then
      !   write(*,*) "z is considered as a pole"
              cgammar = (0d0,0d0)
              return
      end if
      gamm = cgamma(z)
      cgammar = 1d0/gamm
      !write(*,*) cgammar
      return
      end





      !
      !
      !
      !*
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !*     CGAMMA(C)= LOGARITHM OF THE GAMMA FUNCTION (C)
      !.
      !*               COMPLEX DOUBLE PRESITION
      !*
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          complex*16 FUNCTION CGAMMA(C)
          IMPLICIT REAL*8(A-B,D-H,O-Z), COMPLEX*16(C)
              DIMENSION H(6)
          DATA H/1.D0,+5.772157D-01,-6.558781D-01,-4.2002635D-02,
     .      +1.665386D-01,-4.219773D-02/
          DATA ER/1.D-07/
          XC=CDABS(C)
          IF(XC.LT.0.1D0) GOTO 1

 2        DR=XC*ER
          DI=DR
          CALL CLGAMA(C,0,DR,DI,IPOLO,CLGAM)
              CGAMMA=CDEXP(CLGAM)
          RETURN

 1        IF(XC.LT.ER) GOTO 2
          CGAM=C*(H(1)+C*(H(2)+C*(H(3)+C*(H(4)+C*(H(5)+C*H(6))))))
          CGAMMA=1.D0/CGAM

          RETURN
          END

      ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !     CLGAM(C)= LOGARITHM OF THE GAMMA FUNCTION (C)
      !     .
      !               COMPLEX DOUBLE PRESITION
      ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          FUNCTION CLGAM(C)
          IMPLICIT REAL*8(A-B,D-H,O-Z), COMPLEX*16(C)
          DIMENSION H(6)
          DATA H/1.D0,+5.772157D-01,-6.558781D-01,-4.2002635D-02,
     .      +1.665386D-01,-4.219773D-02/
          DATA ER/1.D-07/
          XC=CDABS(C)
          IF(XC.LT.0.1D0) GOTO 1

 2        DR=XC*ER
          DI=DR
          CALL CLGAMA(C,0,DR,DI,IPOLO,CLGAM)
          RETURN

 1        IF(XC.LT.ER) GOTO 2
          CGAM=C*(H(1)+C*(H(2)+C*(H(3)+C*(H(4)+C*(H(5)+C*H(6))))))
          CGAM=1.D0/CGAM
          CLGAM=CDLOG(CGAM)
          RETURN
          END
      ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          SUBROUTINE CLGAMA(ZZ,IWRITE,DR,DI,IPOLO,LOGAM)
          IMPLICIT REAL*8(A-H,O-Z)
          COMPLEX*16 Z,Z2,LOGAM,C,ZA,ZZ
          COMPLEX*16 C0,C1,R,W
          DIMENSION A(2),B(12),G(7)
          EQUIVALENCE (ZA,A(1))
          PARAMETER(C0=(0.D0,0.D0), C1=(1.D0,0.D0))
          PARAMETER(NZERO=0, NONE=1, NSIX=6)
          PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, FOUR=4.D0)
          PARAMETER(A0=9.189385D-1, A8=0.5D0)
          DATA B/+8.333333D-2,-2.777778D-3,+7.936508D-4,          
     .           -5.952381D-4,+8.417508D-4,-1.917527D-3,   
     .           +6.410256D-3,-2.955065D-2,+1.796444D-1,   
     .           -1.392432D+0,+1.340286D+1,-1.568483D+2/
          DATA G/+8.333333D-2,+3.333333D-2,+2.523809D-1,        
     .           +5.256065D-1,+1.011523D+0,+1.517474D+0, 
     .           +2.269489D+0/

          ZA=ZZ
          Z=ZZ
          LOGAM=C0
          IRZ=INT(A(1))
          D=TWO*DR
          IF(A(1)-D) 2,2,1
 2        DM=DABS(A(1)-DFLOAT(IRZ))
          IF(DM.GT.DR) GOTO 1
          IF(DABS(A(2)).GT.DI) GOTO 1
          WRITE(*,3) Z
 3        FORMAT(2X,'Z=',2(1PD14.7,2X),'IS CONSIDERED AS A POLE')
          IPOLE=NONE
          RETURN

 1        CONTINUE
          IPOLE=NZERO
          IF(A(1).GE.ONE) GOTO 20

          NG=NSIX-IRZ
          IF(NG) 4,4,5
 5        C=CMPLX(DFLOAT(NG),ZERO)
          Z=Z+C
 4        Z2=Z*Z
          R=ONE/Z2
          LOGAM=(B(1)+R*(B(2)+R*(B(3)+R*(B(4)+R*(B(5)+R*(B(6)+R*(B(7)+R*
     .        (B(8)+R*(B(9)+R*B(10))))))))))/Z

 8        LOGAM=LOGAM+A0-Z+(Z-A8)*CDLOG(Z)
          IF(NG) 100,100,9

 9        C=C1
          DO 10 I=1,NG
          Z=Z-C
          LOGAM=LOGAM-CDLOG(Z)
 10       CONTINUE
          GOTO 100

 20       CONTINUE
          LOGAM=(Z-A8)*CDLOG(Z)+A0-Z
          LOGAM=LOGAM+G(1)/(Z+G(2)/(Z+G(3)/(Z+G(4)/(Z+G(5)/(Z+G(6)/  
     .                (Z+G(7)/Z))))))

 100      CONTINUE
          IF(IWRITE.EQ.1) WRITE(5,200) ZZ,LOGAM
 200      FORMAT(2X,'Z=',2(1PD14.7,2X),5X,'LOG GAMMA(Z)=',2(1PD14.7,2X))
          RETURN
          END
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !     COULOMB FACTOR
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          FUNCTION COUFA(A)
          IMPLICIT REAL*8(A-B,D-H,O-Z), COMPLEX*16(C)
          PARAMETER(PI05=1.570796D0, ONE=1.D0)
          PARAMETER(CI=(0.D0,1.D0), C1=(1.D0,0.D0))
          C=A*PI05+CLGAM(C1*ONE-CI*A)
          COUFA=CDEXP(C)
          RETURN
          END
