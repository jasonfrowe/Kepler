      subroutine findex(xi,xbar,n,klo,L_one)
      integer n,klo,khi,k,L_one
      real*8 xi(n),xbar
      L_one=0
C find the grid point,  klo, such that xi(klo)<=xbar, and
C abs(xi(klo)-xbar)<1.
      if(xi(1).ge.xbar)then
        klo=1
        khi=2
        go to 522
      endif
      if(xi(n).le.xbar)then
        klo=n
        khi=n
        go to 522
      endif
      klo=1
      khi=n  
    2 if((khi-klo).gt.1)then
         k=(khi+klo)/2.0
           if(xi(k).gt.xbar)then
              khi=k
           else
              klo=k
           endif
         go to 2
      endif
      if((khi-klo).le.0)then
                write(0,'(21h  interpolation error)')
                stop
      endif 
  522 continue
C now, (klo,khi) is sub-range of xi which contains xbar. 
      if(xi(klo).ne.0.0e0)then
c        if(abs((xi(klo)-xbar)/xi(klo)).le.0.5e-5)L_one=1
c        if(abs((xi(klo)-xbar)/xi(klo)).le.0.5e-15)L_one=1
      else
c        if(abs(xi(klo)-xbar).le.0.5e-7)L_one=1
c        if(abs(xi(klo)-xbar).le.0.5e-17)L_one=1
      endif
      return
      end
C track1 nline=24
C track2 nline=150
      subroutine read_track(filename,work,nline,*)
      integer nline,i,k
      real*8 work(5,nline)
      character*72 filename
      character*60 dummy
      common /header/dummy
      open(4, file=filename,status='old',err=999)
      read(4,'(a)')dummy
      do i=1,nline
      read(4,*,end=99,err=99)line,(work(k,i),k=1,5)
      if(i.ne.line)goto 99
      enddo
      close(4)
C      write(0,*)' Input  -> ',filename
      return
  999 close(4)
      return 1
   99 write(0,*)' error... file read in'
      write(0,*)' at ',i
      close(4)
      stop
      end
      subroutine write_track(filename,work,nline)
      integer nline,i,k
      real*8 work(5,nline)
      character*72 filename
      character*60 dummy
      common /header/dummy
      open(4, file=filename,status='new')
      write(4,'(a)')dummy
      do i=1,nline
      write(4,9)i,(work(k,i),k=1,5)
    9 format(i3,1x,5f12.8)
      enddo
      close(4)
      write(0,'(1x,8hOutput: ,a60)')filename
      return
      end 


      subroutine interp_z(z,vmass,pwork,nline,ka,kz,km,kf)
      parameter (inda=3)
      parameter (indz=11)
      parameter (indm=35)
      real*8 avalue(inda)
      real*8 zvalue(indz)
      real*8 xmass(indm)
      common /grid/avalue,zvalue,xmass
      integer ka,kz,km,kf,kvar,nline,i,k,m,mm
       real*8 workz(5,150,4),pwork(5,nline),z,vmass
       real*8 x(4),y(4)
      integer a_one,z_one,m_one
      common /single/a_one,z_one,m_one
      
      kvar=kz
        if(z_one.eq.1)then
         call interp_m(vmass,pwork,nline,ka,kvar,km,kf)
        else
          if(kz+2.gt.indz)kvar=indz-2
          if(kz-1.lt.1)kvar=2
         call interp_m(vmass,workz(1,1,1),nline,ka,kvar-1,km,kf)
         call interp_m(vmass,workz(1,1,2),nline,ka,kvar,km,kf)
         call interp_m(vmass,workz(1,1,3),nline,ka,kvar+1,km,kf)
         call interp_m(vmass,workz(1,1,4),nline,ka,kvar+2,km,kf)
         do 10 i=1,nline
         do 20 k=1,5
            mm=kvar-2
         do 1 m=1,4
          mm=mm+1
          x(m)=dlog(zvalue(mm))
          y(m)=workz(k,i,m)
    1    continue
         call newton(x,y,4,dlog(z),pwork(k,i))
   20 continue
   10 continue
        endif
       return
      end
      subroutine interp_m(vmass,pwork,nline,ka,kz,km,kf)
      parameter (inda=3)
      parameter (indz=11)
      parameter (indm=35)
      real*8 avalue(inda)
      real*8 zvalue(indz)
      real*8 xmass(indm),y(4)
      common /grid/avalue,zvalue,xmass
      real*8 workm(5,150,4),vmass,pwork(5,nline),strx(4)
      character*72 trackname
      integer a_one,z_one,m_one
      common /single/a_one,z_one,m_one
      integer ka,kz,km,kf,kvar,nline,kkip,i,j,k,m

      kkip=km
      if(m_one.eq.1)then
         call make_name(ka,kz,kkip,kf,trackname)
         call read_track(trackname,pwork,nline,*999)
         return
      endif
   99 continue
      kvar=kkip
      if(kkip+2.gt.indm)kvar=indm-2
      if(kkip-1.lt.1)kvar=2
      kvar=kvar-2
      do 8 j=1,4
    9    kvar=kvar+1
         if(kvar.gt.indm)goto 999
         call make_name(ka,kz,kvar,kf,trackname)
         call read_track(trackname,workm(1,1,j),nline,*9)
         strx(j)=xmass(kvar)
    8 continue

         do 10 i=1,nline
         do 20 k=1,5
         do 1 m=1,4
          y(m)=workm(k,i,m)
    1    continue
         call newton(strx,y,4,vmass,pwork(k,i))
   20 continue
   10 continue
      return
  999 continue
C      m_one=0
      kkip=kkip-1
      goto 99
c      stop ' no such file exist...'
      end

      subroutine newton(x,yy,n,xbar,ybar)
      real*8 x(n), y(10),xbar,ybar,yy(n)
      integer n,k,l ,m
c ,np
      do k=1,n
       y(k)=yy(k)
      enddo
      do 10 k=1,n-1
      do 20 l=1,(n-k)
        y(l)=(y(l+1)-y(l))/(x(l+k)-x(l))
   20 continue
   10 continue
C      do 30 np=1,n
C      write(0,*)y(np)
C   30 continue
      ybar=y(1)
      do 100 m=2,n
       ybar=ybar*(xbar-x(m))+y(m)
  100 continue
C      write(0,*)xbar,ybar
      return
      end
      subroutine make_name(ka,kz,km,kf,filename)
      parameter (inda=3)
      parameter (indz=11)
      parameter (indm=35)
      character*4 aname(inda)
      data aname/'a0o2','a2o2','a4o2'/
      integer zlength(indz)
      data zlength/12,2*10,3*8,5*6/
      character*12 zname(indz)
      data zname/ 'x76997z00001', 'x7697z0001', 'x7688z0004',
     +'x767z001', 'x758z004', 'x749z007', 'x74z01',
     +'x71z02', 'x65z04', 'x59z06', 'x53z08'/
      character*3 xmname(indm)
      data xmname/'m04','m05','m06','m07','m08','m09','m10',
     +'m11','m12','m13','m14','m15','m16','m17','m18','m19','m20',
     +'m21','m22','m23','m24','m25','m26','m27','m28','m29','m30',
     +'m32','m34','m36','m38','m40','m42','m45','m50'/
      character*7 ftype(2)
      data ftype/'.track1','.track2'/
      character*80 fmt(1)

      character*72 filename,trackname
      integer ka,kz,km,kf

      if(kz.le.3)then
      write(fmt,11)zlength(kz),zlength(kz)
   11 format('(2h./,a4,1h/,a',i2,',1h/,a3,a',i2,',a7)')
      else
      write(fmt,12)zlength(kz),zlength(kz)
   12 format('(2h./,a4,1h/,a',i1,',1h/,a3,a',i1,',a7)')
      endif
C
      write(trackname,fmt)aname(ka),zname(kz),
     +                    xmname(km),zname(kz),ftype(kf)
      filename=trackname
      return
      end
