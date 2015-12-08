CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function modecal(npt,pts,dmin,dmax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nbins,i,j
      parameter(nbins=200)
      integer bins(nbins),maxbin,maxi
      double precision pts(npt),dmin,dmax,bspace,mode,x(nbins),y(nbins),
     .  minx,maxx,miny,maxy

      do 10 i=1,nbins
         bins(i)=0
 10   continue
      maxbin=-99

      modecal=300.0

      do 20 j=1,npt
C     calculating the bspacing between each bin
         bspace=(dmax-dmin)/dble(nbins)
C     calculating the bin number
         i=(pts(j)-dmin)/bspace
         if((i.gt.0).and.(i.le.nbins)) then
c            write(6,*) i
            bins(i)=bins(i)+1
            if(bins(i).gt.maxbin) then 
               maxbin=bins(i)
               maxi=i
               mode=dble(i)*bspace+dmin
            endif
         endif
 20   continue

      modecal=mode

c      write(6,*) "mode,nbin",modecal,maxbin
      return
      end