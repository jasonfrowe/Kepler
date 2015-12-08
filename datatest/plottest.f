      program plottest
      implicit none
      integer nmax,i,seed
      parameter(nmax=64)
      integer ia(1,64)
      real x,y


      call pgopen('?') !open plotting device

      call PGPAP(8.0,1.0) !set paper size and aspect ratio

      call pgsubp(2,2) !make four panels

      do 10 i=1,nmax
         call pgscr(i+15,real(i)/64.0,0.0,0.0) !define colours
 10   continue

      call pgpage() !start with first panel

      call pgvport(0.3,0.9,0.2,0.8) !how much of the panel to use

      call pgwindow(0.0,1.0,0.0,1.0) !sets x,y scale
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0) !draw ticks
      call pglabel("X data","Y data","")

      seed=10
      x=ran(seed) !initialize random numbers
C     Here I just plot some points..
      do 11 i=1,nmax
         call pgsci(i+15)
         x=ran(seed)
         y=ran(seed)
         call pgpt1(x,y,17)
 11   continue
      call pgsci(1) !reset colour


C     Okay, lets move to the next panel..
      call pgpage()

      call pgvport(0.05,0.15,0.2,0.8) !how much of the panel to use

C     Generate ia array
      do 12 i=1,nmax
         ia(1,i)=i+15
 12   continue

      call pgpixl(ia,1,nmax,1,1,1,nmax,0.0,1.0,0.0,1.0)

      call pgwindow(0.0,1.0,0.0,10.0) !sets x,y scale
      call pgbox('BC',0.0,0,'BCNTS1',0.0,0) !draw ticks
      call pglabel("","Values","") !put some labels on it

      call pgclos() !close plotting device

      end

