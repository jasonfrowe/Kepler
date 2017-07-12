CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function eclmod2(R1,R2,x1,x2,y1,y2,zarea,norm,Pi,
     .  ted)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer sflag
      double precision R1,R2,x1,x2,y1,y2,sinth,x,y,dist,d1,d2,ratio,
     .  x01,x02,ax1,ax2,ay1,ay2,zarea,norm,Pi,ted,y2p
      
C     Put everything in upper quadrant to make life simple
      ax1=abs(x1)
      ax2=abs(x2)
      ay1=abs(y1)
      ay2=abs(y2)
      
C     distance between projected center of planet and star
      dist=sqrt( (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1) )
      
      call radsolve(R1,R2,x1,x2,y1,y2,x01,x02,sflag)
c      if(dist.lt.R1+R2)then !we have an eclipse 
      if(sflag.eq.0)then
C     get x,y which is the intercept of the stellar radius with a 
C     diameter of the stellar radius pointed towards the stellar radius
C     origin
        sinth=abs(ax2-ax1)/(dist)
        x=ax1+abs(R1*sinth)
        y=ay1+sqrt(R1*R1-(x-ax1)*(x-ax1))
        d1=sqrt( (x-ax1)*(x-ax1)+(y-ay1)*(y-ay1) )
        d2=sqrt( (x-ax2)*(x-ax2)+(y-ay2)*(y-ay2) )
c        write(6,*) "s:",sinth,x/R1,y/R1
c        write(6,*) "r:",d2/R2
c        write(6,*) "d:",d1/R1,d2/R1
        if(d1.le.dist)then
            ratio= (R2-d2)/(2.0*R2)
        else
            ratio= (R2+d2)/(2.0*R2)
        endif
c        read(5,*)
c        eclmod2=(Pi+zarea*ratio)/norm-ted*ratio
        eclmod2=(Pi+zarea)/norm-ted*ratio
      else
         y2p=Sqrt(x2*x2+y2*y2)
         if(y2p/R1+R2/R1.lt.1.0d0)then
            eclmod2=(Pi+zarea)/norm-ted!inside transit
c            write(6,*) 0,y2p/R1,eclmod
         else
c            eclmod2=(Pi+zarea)/norm!outside transit
            eclmod2=(Pi+zarea)/norm!outside transit
         endif     

      endif
        
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function eclmod(R1,R2,x1,x2,y1,y2,zarea,norm,
     .   Pi,xintold,y2pold,ted)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer sflag
      double precision cirfunc,th,getangle,x01,x02,y01,y02,xint,
     .   x2p,y2p,xintold,y2pold,xint1,xint2,isign,area1,area2,areaint,
     .   darea,zarea,norm,Pi,R1,R2,x1,x2,y1,y2,ted,tedoff
      
      call radsolve(R1,R2,x1,x2,y1,y2,x01,x02,sflag)
      if(sflag.eq.0)then !if there is a transit or eclipse.. then
         y01=cirfunc(R1,x1,y1,x01)!get y co-ordinates of intersection
         y02=cirfunc(R1,x1,y1,x02)! points
         th=getangle(R1,y01,y02)! angle between intesection points
         xint=Cos(Pi/2.0-th) !rotation angle
         x2p=0. !rotated x-co-ordinate
         y2p=Sqrt(x2*x2+y2*y2) !rotated y co-ordinate
C        case for minimum area of planet transiting/eclipsing
         if(((xint.gt.xintold).and.(y2p.lt.y2pold)).or.
     .      ((xint.lt.xintold).and.(y2p.gt.y2pold)))then
c         if(y2p/R1.lt.1.0d0)then
            xint1=-xint  !set integration range
            xint2= xint
            isign=1.0d0  !integrate over top of circle
C           integrate to find star area
            area1=areaint(1.0d0,x1,y1,xint1,xint2,isign)
            isign=-1.0d0 !integrate over bottom of circle
C           integrate to find planet area
            area2=areaint(R2/R1,x2p/R1,y2p/R1,xint1,xint2,isign)
C           gives effective eclipse area
            darea=Pi+zarea*(1.0-(area1-area2)/(Pi*R2*R2/(R1*R1)))
            tedoff=ted*(1.0-(area1-area2)/(Pi*R2*R2/(R1*R1)))
C        case for maximum area of planet transiting/eclisping
         else
            xint1=-xint !set integration range
            xint2= xint
            isign=1.0d0 !integrate over top of curve
            area1=areaint(1.0d0,x1,y1,xint1,xint2,isign)
            area2=areaint(R2/R1,x2p/R1,y2p/R1,xint1,xint2,isign)
C           calculate effective transit area
            darea=Pi+zarea*(area2-area1)/(Pi*R2*R2/(R1*R1))
            tedoff=ted*(area2-area1)/(Pi*R2*R2/(R1*R1))
         endif          
         eclmod=darea/norm-tedoff
c         write(6,500) x2/R1,y2/R1,x01/R1,y01/R1,x02/R1,y02/R1,th,eclmod
 500     format(28(F7.4,1X))
c         write(6,*) tedoff,xintold,y2pold
         xintold=xint
         y2pold=y2p
      else
         y2p=Sqrt(x2*x2+y2*y2)
         if(y2p/R1+R2/R1.lt.1.0d0)then
            eclmod=(Pi)/norm-ted!inside transit
c            write(6,*) 0,y2p/R1,eclmod
         else
            eclmod=(Pi+zarea)/norm!outside transit
         endif     
         xintold=0.   
      endif
         
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine radsolve(R1,R2,x1,x2,y1,y2,x01,x02,sflag)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer sflag
      double precision R1,R2,x1,x2,y1,y2,x01,x02,a,b,temp(3)
      
      a=x1-x2
      b=y1-y2
      
      temp(1)=b*b*((R1-R2)*(R1-R2)-a*a-b*b)*(-(R1+R2)*(R1+R2)+a*a+b*b)
     
      if(temp(1).gt.0)then
         sflag=0
         temp(2)=2.0d0*(a*a+b*b)
         temp(3)=-R1*R1*a+R2*R2*a+(x1+x2)*(a*a+b*b)
         x01=(temp(3)+Sqrt(temp(1)))/temp(2)
         x02=(temp(3)-Sqrt(temp(1)))/temp(2)
      else
         sflag=1
      endif
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function cirfunc(R1,x1,y1,x)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision R1,x1,y1,x
      
      cirfunc=Sqrt(R1*R1-(x-x1)*(x-x1))+y1
      
      return
      end   
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function getangle(R1,y01,y02)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision R1,y01,y02,th1,th2
    
      th1=Asin(y01/R1)
      th2=Asin(y02/R1)
      
      getangle=Abs(th2-th1)/2.0d0
    
      return
      end 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function areaint(R,x1,y1,xint1,xint2,isign)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision R,x1,y1,xint1,xint2,int1,int2,x0,isign,temp,Pi,
     .   isign2
      Pi=acos(-1.d0)
      
c      x0=xint1
c      ch1=(Sqrt(R**2 - (x0 - x1)**2)*(x0 - x1))/
c     -       (-R**2 + (x0 - x1)**2)
c      x0=xint2
c      ch2=(Sqrt(R**2 - (x0 - x1)**2)*(x0 - x1))/
c     -       (-R**2 + (x0 - x1)**2)
      
      x0=xint1
      temp=-R**2 + (x0 - x1)**2
      if(temp.eq.0.0d0)then
         isign2=(x0-x1)/abs(x0-x1)
         int1=isign*(Sqrt(R**2-(x0 - x1)**2)*(x0 - x1))/2.0d0 + x0*y1 - 
     -  isign*(isign2*R**2*Pi/4.0d0)/2.0d0
      else
         int1=isign*(Sqrt(R**2-(x0 - x1)**2)*(x0 - x1))/2.0d0 + x0*y1 - 
     -  isign*(R**2*Atan((Sqrt(R**2 - (x0 - x1)**2)*(x0 - x1))/
     -       (-R**2 + (x0 - x1)**2)))/2.0d0
      endif
      x0=xint2
      temp=-R**2 + (x0 - x1)**2
      if(temp.eq.0.0d0)then
         isign2=(x0-x1)/abs(x0-x1)
         int2=isign*(Sqrt(R**2-(x0 - x1)**2)*(x0 - x1))/2.0d0 + x0*y1 - 
     -  isign*(isign2*R**2*Pi/4.0d0)/2.0d0
         areaint=int2-int1
      else
         int2=isign*(Sqrt(R**2-(x0 - x1)**2)*(x0 - x1))/2.0d0 + x0*y1 - 
     -  isign*(R**2*Atan((Sqrt(R**2 - (x0 - x1)**2)*(x0 - x1))/
     -       (-R**2 + (x0 - x1)**2)))/2.0d0
      endif
      
      
      areaint=int2-int1
      
      return
      end