CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plottrans(n,x,y,per,phase,qtran)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer i,n
      real bb(4),px,py
      double precision x(n),y(n),per,phase,qtran

      call pgvport(0.5,0.95,0.2,0.8) !gives enough room for labels

      bb(1)=real(max(0.5-2.0*qtran,0.0))
      bb(2)=real(min(0.5+2.0*qtran,1.0))
      bb(3)=y(1)
      bb(4)=y(1)
      do 10 i=1,n
        bb(3)=min(bb(3),y(i))
        bb(4)=max(bb(4),y(i))
 10   continue

      call pgwindow(bb(1),bb(2),bb(3),bb(4))
c      call pgbox('BNTS1',0.0,0,'BCNTSV1',0.0,0)
      call pgptxt(bb(1)-0.15*(bb(2)-bb(1)),(bb(3)+bb(4))/2.0,90.0,0.5,
     .   "flux")
      call pgptxt((bb(2)+bb(1))/2.0,bb(3)-0.25*(bb(4)-bb(3)), 0.0,0.5,
     .   "phase (hours)")

      call pgbbuf()
      do 30 i=1,n
        px=real(x(i)/per-int(x(i)/per)-phase+0.5)
        if(px.lt.0.0)px=px+1
        py=real(y(i))
        call pgpt1(px,py,1)
 30   continue
      call pgebuf()

      call pgwindow(real(-48.0*qtran*per),real(48.0*qtran*per),
     .   bb(3),bb(4))
      call pgbox('BNTS1',0.0,0,'BCNTSV1',0.0,0)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotph(n,x,y,per,phase)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer n,i,stepmax
      parameter(stepmax=1200)
      real px,py,bb(4),lx(2),ly(2)
      double precision x(n),y(n),per,phase
      
      bb(1)=0.0
      bb(2)=1.0
      bb(3)=y(1)
      bb(4)=y(1)
      do 10 i=1,n
        bb(3)=min(bb(3),y(i))
        bb(4)=max(bb(4),y(i))
 10   continue
   
      call pgwindow(bb(1),bb(2),bb(3),bb(4))
      call pgbox('BNTS1',0.0,0,'BCNTSV1',0.0,0)
      call pgptxt(bb(1)-0.08*(bb(2)-bb(1)),(bb(3)+bb(4))/2.0,90.0,0.5,
     .   "flux")
      call pgptxt((bb(2)+bb(1))/2.0,bb(3)-0.25*(bb(4)-bb(3)), 0.0,0.5,
     .   "phase")
      
      lx(1)=real(phase)
      lx(2)=real(phase)
      ly(1)=bb(3)
      ly(2)=bb(4)
      call pgsci(2)
      call pgline(2,lx,ly)
      call pgsci(1)
      
      call pgbbuf()
      do 30 i=1,n
        px=real(x(i)/per-int(x(i)/per))
        if(px.lt.0.0)px=px+1
        py=real(y(i))
        call pgpt1(px,py,-1)
 30   continue
      call pgebuf()
 
      return
      end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plot(n,x,y,title,xlabel,ylabel,nmark,epo,per)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer n,i,stepmax,nmark
      parameter(stepmax=1200)
      real px(3),py(3),bb(4),xp,yp
      double precision x(n),y(n),epo,per
      character*(*) title,xlabel,ylabel
      
      bb(1)=x(1)
      bb(2)=x(1)
      bb(3)=y(1)
      bb(4)=y(1)
      do 10 i=1,n
        bb(1)=min(bb(1),x(i))
        bb(2)=max(bb(2),x(i))
        bb(3)=min(bb(3),y(i))
        bb(4)=max(bb(4),y(i))
 10   continue
   
      call pgwindow(bb(1),bb(2),bb(3),bb(4))
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
c      call pglabel(xlabel,ylabel,title)
      call pgptxt(bb(1)-0.08*(bb(2)-bb(1)),(bb(3)+bb(4))/2.0,90.0,0.5,
     .   ylabel)
      call pgptxt((bb(2)+bb(1))/2.0,bb(3)-0.23*(bb(4)-bb(3)), 0.0,0.5,
     .   xlabel)
      call pgptxt((bb(2)+bb(1))/2.0,bb(4)+0.10*(bb(4)-bb(3)), 0.0,0.5,
     .   title)
      
      call pgbbuf()
      do 30 i=2,n
        px(1)=real(x(i-1))
        py(1)=real(y(i-1))
        px(2)=real(x(i))
        px(3)=real(x(i))
        py(2)=real(y(i))
        py(3)=real(y(i))
        call pgline(3,px,py)
 30   continue
      call pgebuf()
 
      if(nmark.eq.1)then
         xp=bb(1)
         i=0
         yp=bb(3)+0.05*(bb(4)-bb(3))
         do while((xp.lt.bb(2)).or.(i.lt.10000))
            xp=real(epo+per*dble(i))
            call pgsci(2)
            call pgpt1(xp,yp,13)
            call pgsci(1)
            i=i+1
         enddo
      endif

      return
      end
