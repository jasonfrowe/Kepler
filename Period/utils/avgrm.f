CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine avgrm(npt,mag,rmavg)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Removes average value from data points
C     This can help reduce numerical overflow for large data sets 
      implicit none
      integer npt
      double precision mag(npt),rmavg

      integer i
      double precision ave,adev,sdev,sigma2,skew,curt

      call moment(mag,npt,ave,adev,sdev,sigma2,skew,curt)

      do 10 i=1,npt
         mag(i)=mag(i)-ave
 10   continue

      rmavg=ave
      
      return
      end