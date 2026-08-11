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
      integer npmax,nmax,npt(npmax),np,jlow,jmid,jhigh
      double precision x(npmax,nmax),y(npmax,nmax),xin,yout
      
C     Default is zero
      yout=0.0d0
      if(npt(np).eq.0) return

      if(xin.le.x(np,1))then
        if(npt(np).gt.1) then
          yout=y(np,1)+(xin-x(np,1))/(x(np,2)-x(np,1))*(y(np,2)-y(np,1))
        else
          yout=y(np,1)
        endif
      elseif(xin.ge.x(np,npt(np)))then
        if(npt(np).gt.1) then
          yout=y(np,npt(np))+(xin-x(np,npt(np)))/(x(np,npt(np))-
     .        x(np,npt(np)-1))*(y(np,npt(np))-y(np,npt(np)-1))
        else
          yout=y(np,npt(np))
        endif
      else
        jlow=1
        jhigh=npt(np)
 10     if(jhigh-jlow.gt.1) then
          jmid=(jhigh+jlow)/2
          if(xin.gt.x(np,jmid)) then
            jlow=jmid
          else
            jhigh=jmid
          endif
          goto 10
        endif
        yout=y(np,jlow)+(xin-x(np,jlow))/(x(np,jhigh)-x(np,jlow))*
     .       (y(np,jhigh)-y(np,jlow))
      endif
 
      return
      end