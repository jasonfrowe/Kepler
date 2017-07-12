subroutine pixelfix(np,xp,yp,xpos,ypos)
use precision
implicit none
integer :: np,i,j,currentpixel,npix,npixmax,ijump,pixel
real(double) :: fpix,df,pixthres,jump,outcut,dfp,dfold
real(double), dimension(np) :: xp,yp,xpos,ypos
real(double), allocatable, dimension(:) :: lx,ly

pixthres=0.3
outcut=2.0 !removing outliers.

npixmax=100
allocate(lx(npixmax),ly(npixmax))

!we assume data is sorted in time.
currentpixel=int(xpos(1))
ijump=1
jump=0.0d0
dfold=0.0d0
do i=2,np-1
   fpix=xpos(i)-int(xpos(i)+0.5d0)
   df=yp(i)-yp(i-1)
   dfp=yp(i)-yp(i+1)
   pixel=int(xpos(i))
   if(abs(fpix).le.pixthres)then
!      write(0,*) xpos(i),fpix,df
!      read(5,*)
      if((abs(df).gt.abs(jump)).and.(abs(df).gt.abs(dfp*outcut)) )then  !     &
       !.and.(abs(df).gt.abs(dfold*outcut)))then
         ijump=i
         jump=df
      endif
   endif
   if((pixel.gt.currentpixel).and.(fpix.gt.pixthres))then
      write(0,*) xp(i),currentpixel,ijump,jump
!      read(5,*)
!     reset counters

      do j=ijump,np
         yp(j)=yp(j)-jump
      enddo

      jump=0.0d0
      currentpixel=pixel
   endif
   dfold=df
enddo
500 format(108(1X,1PE17.10))
501 format(108(1X,1PE10.3))


return
end subroutine pixelfix
