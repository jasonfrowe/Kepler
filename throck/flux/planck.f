CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function planck(T,npass,lam,pass)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     I believe that lam input is in nm, then converted to cm below.
      implicit none
      integer npass,i
      double precision T,l,h,c,k,temp(2),Pi,lam(npass),pass(npass),dl
      data k /1.380658d-16/
      data c /2.99792458d10/
      data h /6.6260755d-27/
      data Pi /3.141592654/

      planck=0. !initialize for integration
      do 10 i=2,npass
       
        l=lam(i)*1.0d-7 !convert nm -> cm
        dl=(lam(i)-lam(i-1))*1.0d-7 !calculate delta lambda
      
        temp(1)=2.0d0*Pi*h*c*c/(l**5.0d0) !planck at specific wavelength 
        temp(2)=exp(h*c/(l*k*T))-1.0d0
      
        planck=planck+temp(1)/temp(2)*dl*(pass(i)+pass(i-1))/2.0d0
 10   continue
      
      return
      end