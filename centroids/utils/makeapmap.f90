subroutine makeapmap(nstar,xcoo,ycoo,naxes,Image,bpix,neimax,napmap,    &
 apmap)
use precision
implicit none
integer :: nstar,neimax
integer, dimension(2) :: naxes
integer, dimension(:,:) :: napmap
integer, dimension(:,:,:) :: apmap
real(double) :: bpix
real(double), dimension(:) :: xcoo,ycoo
real(double), dimension(:,:) :: Image
!local vars
integer :: i,j,k,xc,yc,dr,nnew
real(double) :: ap,newap

!initialize map to zero
apmap=0
napmap=0

do i=1,nstar
   !central pixel
   xc=int(xcoo(i))
   yc=int(ycoo(i))
!   write(0,*) "xc,yc: ",xc,yc
   napmap(xc,yc)=napmap(xc,yc)+1
   if(napmap(xc,yc).gt.neimax)then
      write(0,*) "napmap greater than neimax"
      stop
   endif
   apmap(xc,yc,napmap(xc,yc))=i !assign central pixel to star index
   ap=Image(xc,yc) !update photometry
   dr=1 !radius to expand search
   nnew=1 !number of new pixels added to aperture
   do while(nnew.gt.0)
      nnew=0
      newap=0.0d0
      do j=max(1,xc-dr),min(naxes(1),xc+dr)
         do k=max(1,yc-dr),min(naxes(1),yc+dr)
            if((j.eq.xc-dr).or.(j.eq.xc+dr).or.(k.eq.yc-dr).or.(k.eq.yc+dr))then
               if(Image(j,k).lt.bpix)then
                  if(Image(j,k).gt.sqrt(ap))then !is the extra flux more than noise?
                     nnew=nnew+1
                     newap=newap+Image(j,k)
                     napmap(j,k)=napmap(j,k)+1
                     if(napmap(xc,yc).gt.neimax)then
                        write(0,*) "napmap greater than neimax"
                        stop
                     endif
                     apmap(j,k,napmap(j,k))=i
!                     write(0,*) "new: ",j,k,Image(j,k)
                  endif
               endif
            endif
         enddo
      enddo
      ap=ap+newap !update aperture flux
      dr=dr+1
!      write(0,*) "nnew,nextdr: ",nnew,dr
!      read(5,*)
   enddo
enddo

return
end subroutine makeapmap
