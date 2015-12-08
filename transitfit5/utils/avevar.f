CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE avevar(data,n,ave,var)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n,m
      REAL*8 ave,var,data(n)
C Given array data(1:n), returns its mean as ave and its variance as var.
      INTEGER j
      REAL*8 s,ep
      ave=0.0
      m=min(1000,n)
      do 11 j=1,m
         ave=ave+data(j)
 11   continue
      ave=ave/m
      var=0.0
      ep=0.0
      do 12 j=1,n
         s=data(j)-ave
         ep=ep+s
         var=var+s*s
 12   continue
      var=(var-ep**2/n)/(n-1) !Corrected two-pass formula (14.1.8).
      return
      END 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE avevar2(data,data2,n,ave,var)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n,m,n2
      REAL*8 ave,var,data(n),std,data2(n),i
C Given array data(1:n), returns its mean as ave and its variance as var.
      INTEGER j
      REAL*8 s,ep
c      ave=0.0
c      m=min(1000,n)
c      do 11 j=1,m
c         ave=ave+data(j)
c 11   continue
c      ave=ave/m
      var=0.0
      ep=0.0
      do 12 j=1,n
         s=data(j)-ave
         ep=ep+s
         var=var+s*s
 12   continue
      var=(var-ep**2/n)/(n-1) !Corrected two-pass formula (14.1.8).
      std=sqrt(var)
c      write(0,*) "ave: ",ave,var,std

      do 15 i=1,1!3

      n2=0
      do 13 j=1,n
c         write(0,*) "c: ",data(j),ave,std
         if(abs(data(j)-ave).lt.3.0*std)then
            n2=n2+1
            data2(n2)=data(j)
         endif
 13   continue

c      write(0,*) "n2: ",n2

      var=0.0
      ep=0.0
      do 14 j=1,n2
         s=data2(j)-ave
         ep=ep+s
         var=var+s*s
 14   continue
      var=(var-ep**2/n2)/(n2-1) !Corrected two-pass formula (14.1.8).

 15   continue

      return
      END
