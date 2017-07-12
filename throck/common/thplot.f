CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine opengraphics(nlon,nlat,naxes)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nlon,nlat,naxes(2)

      naxes(1)=nlon  !defines the plotting surface.
      naxes(2)=nlat
            
C     opens plotting device
      call pgopen('/xserve')
c      call pgopen('test.ps/vcps')
      call PGPAP ( 6.0 ,1.0)  
      call PGENV ( -0.5 , 0.5 , -0.5 , 0.5 , 1 , -2 )  
C     turns on prompts
      call pgask(.false.)
C     sets the line width
      call pgslw(3)
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine closegraphics()
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      
      call pgclos()
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine displayfits(nlat,nlon,ia,pts,fitsdata,naxes,datamin,
     .   datamax,time)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer naxes(2),xmax,ymax,i,j,npt,nx,ny,ncol,dumi,nlat,nlon
      parameter(ncol=64)
      integer ia(nlat,nlon)
      real fitsdata(nlat,nlon),datamin,datamax,pts(nlat*nlon),z1,z2,R,G,
     .     B,x1,x2,y1,y2,dumr,time
      character*80 tstring
      xmax=nlat
      ymax=nlon
      
c      call pgpage()

      npt=0
      do 8 i=1,naxes(1),4
         do 9 j=1,naxes(2),4
            npt=npt+1
            pts(npt)=fitsdata(i,j)
 9       continue
 8    continue
c      call autoscale(npt,pts,z1,z2,0.95)
      z1=datamin
      z2=datamax
c      write(0,*) "Pl:",z1,z2
 500  format(A16,2(1X,1PE12.5))
      
      nx=naxes(1)
      ny=naxes(2)

      do 10 i=1,nx
         do 11 j=1,ny
            IA(i,j)=int((fitsdata(i,j)-z1)/(z2-z1)*(NCOL-1))+16
            if(fitsdata(i,j).lt.z1) ia(i,j)=16
            if(fitsdata(i,j).gt.z2) ia(i,j)=ncol+15
 11      continue
 10   continue

C     set up pgplot window
      call pgscr(0,0.0,0.3,0.2)
      call pgsvp(0.07,0.95,0.05,0.95)
      call pgwnad(0.0,1.0,0.0,1.0)
     
      open(unit=12,file="/iraf/iraf/unix/sun/heat.lut",
     .     status='old')

      read(12,*) dumi

C     Setup colour-map
      do 300 i=1,ncol
         read(12,*) r,g,b
         read(12,*) dumr,dumr,dumr
         read(12,*) dumr,dumr,dumr
         read(12,*) dumr,dumr,dumr
c         R = REAL(I-1)/REAL(NCOL-1)*0.8 + 0.2
c         G = MAX(0.0, 2.0*REAL(I-1-NCOL/2)/REAL(NCOL-1))
c         B = 0.2 + 0.4*REAL(NCOL-I)/REAL(NCOL)
c         write(6,*) "cc:",ncol-i+16,r,g,b
c         CALL PGSCR(ncol-I+16, R, G, B)
          CALL PGSCR(I+15, R, G, B)
 300  CONTINUE

      close(12)

 
CC     Setup colour-map
c      do 300 i=1,ncol
c         R = REAL(I-1)/REAL(NCOL-1)*0.8 + 0.2
c         G = MAX(0.0, 2.0*REAL(I-1-NCOL/2)/REAL(NCOL-1))
c         B = 0.2 + 0.4*REAL(NCOL-I)/REAL(NCOL)
c         CALL PGSCR(I+15, R, G, B)
c 300  CONTINUE

      call pgpixl(ia,xmax,ymax,1,nx,1,ny,
     .     0.0,real(nx)/real(max(nx,ny)),
     .     0.0,real(ny)/real(max(nx,ny)))

      call pgqvp(0,x1,x2,y1,y2)
      call pgvport(x1,x1+(x2-x1)*real(nx)/real(max(nx,ny)),y1,
     .   y1+(y2-y1)*real(ny)/real(max(nx,ny)))
c      call pgwindow(1.0,real(nx)+1.0,1.0,real(ny)+1.0)
c      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pgwindow(-180.0*3600.0,180.0*3600.0,-90.0*3600.0,90.0*3600.0)
      call PGTBOX('BCSTNZD',0.0,0,'BCSTNZVD',0.0,0)
      write(tstring,501) time/86400.0
 501  format(F7.2)
      call pgsci(4)
      call PGPTXT (130.0*3600.0,75.0*3600.0, 0.0, 0.0, tstring)
      call pgsci(1)
      
      return
      end
