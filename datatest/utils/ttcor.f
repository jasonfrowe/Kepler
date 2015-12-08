CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      subroutine readttfile(nunit,nplanetmax,nmax,nplanet,ntt,tobs,omc)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,ntt(nplanet),i,nmax,nplanetmax,nplanet
      double precision tobs(nplanetmax,nmax),omc(nplanetmax,nmax),err
      
      i=1
 10   read(nunit,*,end=11) tobs(nplanet,i),omc(nplanet,i),err
        if(err.eq.0.0d0) goto 10
        i=i+1
      goto 10
 11   continue
      ntt(nplanet)=i-1
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine lininterp(x,y,npmax,nmax,np,npt,xin,yout)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer i,npmax,nmax,npt(npmax),np
      double precision x(npmax,nmax),y(npmax,nmax),xin,yout
      
C     Default is zero
      yout=0.0d0
      if(npt(np).eq.0) return

c      write(0,*) (rhoin(i),i=1,npt)
c      write(0,*) (rhoierr(i),i=1,npt)
      if(xin.lt.x(np,1))then
        yout=y(np,1)+(xin-x(np,1))/(x(np,2)-x(np,1))*(y(np,2)-y(np,1))
      elseif(xin.gt.x(np,npt(np)))then
        yout=y(np,npt(np))+(xin-x(np,npt(np)))/(x(np,npt(np))-
     .      x(np,npt(np)-1))*(y(np,npt(np))-y(np,npt(np)-1))
      else
        do 10 i=1,npt(np)-1
            if((xin.gt.x(np,i)).and.(xin.le.x(np,i+1)))then
                yout=y(np,i)+(xin-x(np,i))/(x(np,i+1)-x(np,i))*
     .              (y(np,i+1)-y(np,i))
            endif
 10     continue
      endif
c      write(0,*) drho,dsig
c      read(5,*)
 
      return
      end