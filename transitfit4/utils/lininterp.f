      subroutine lininterp(x,y,npt,xin,yout)
      implicit none
      integer npt,i
      double precision x(npt),y(npt),xin,yout
      
C     Default is zero
      yout=0.0d0

c      write(0,*) (rhoin(i),i=1,npt)
c      write(0,*) (rhoierr(i),i=1,npt)
      if(xin.lt.x(1))then
        yout=y(1)+(xin-x(1))/(x(2)-x(1))*(y(2)-y(1))
      elseif(xin.gt.x(npt))then
        yout=y(npt)+(xin-x(npt))/(x(npt)-x(npt-1))*(y(npt)-y(npt-1))
      else
        do 10 i=1,npt-1
            if((xin.gt.x(i)).and.(xin.le.x(i+1)))then
                yout=y(i)+(xin-x(i))/(x(i+1)-x(i))*(y(i+1)-y(i))
            endif
 10     continue
      endif
c      write(0,*) drho,dsig
c      read(5,*)
 
      return
      end