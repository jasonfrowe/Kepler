subroutine exportfit(nplanet,solout,serr,err,y,yerr,yserr)
use precision
implicit none
integer :: nplanet,i,j,k,nunit,filestatus
real(double) :: yscale
real(double), dimension(:) :: solout,err,y,yerr
real(double), dimension(:,:) :: serr,yserr
character(80) :: exportfile
character(3) :: titles(9)
character(2) :: nbtitles(7)
titles(1)='RHO'
titles(2)='NL1'
titles(3)='NL2'
titles(4)='NL3'
titles(5)='NL4'
titles(6)='DIL'
titles(7)="VOF"
titles(8)="ZPT"
titles(9)="RDR"
nbtitles(1)='MA'
nbtitles(2)='XP'
nbtitles(3)='YP'
nbtitles(4)='ZP'
nbtitles(5)='XV'
nbtitles(6)='YV'
nbtitles(7)='ZV'

nunit=10
exportfile="newfit.dat"
open(unit=nunit,file=exportfile,iostat=filestatus)

if(filestatus>0)then
   write(0,*) "Cannot open ",exportfile
   stop
endif

do i=1,8
   write(nunit,500) titles(i),solout(i),serr(i,1),serr(i,2),err(i)
enddo
500 format(A3,5(1X,1PE17.10))

!Write out star parameters
if(mserr(1,2).gt.0.0d0) mserr(1,2)=mserr(1,2)*MU/Mearth
write(nunit,501) nbtitles(1),1,m(1)*MU/Mearth,mserr(1,1),mserr(1,2),merr(1)*MU/Mearth
do i=1,6
   if(i.lt.4)then
      yscale=1.0d0*LU/AU
   else
      yscale=1.0d0/TU*LU
   endif
   if(yserr(i,2).gt.0.0d0) yserr(i,2)=yserr(i,2)*yscale
   write(nunit,501) nbtitles(i+1),1,y(i)*yscale,yserr(i,1),yserr(i,2),yerr(1)*yscale
enddo
write(nunit,501) titles(9),1,solout(9),serr(9,1),serr(9,2),err(9)

501 format(A2,I1,5(1X,1PE17.10))

do i=1,nplanet
   if(mserr(i+1,2).gt.0.0d0) mserr(i+1,2)=mserr(i+1,2)*MU/Mearth
   write(nunit,501) nbtitles(1),i+1,m(i+1)*MU/Mearth,mserr(i+1,1),mserr(i+1,2), &
      merr(i+1)*MU/Mearth
   do k=1,6
      if(k.lt.4)then
         yscale=1.0d0*LU/AU
      else
         yscale=1.0d0/TU*LU
      endif
      if(yserr(6*(i)+k,2).gt.0.0d0) yserr(6*(i)+k,2)=yserr(6*(i)+k,2)*yscale
      write(nunit,501) nbtitles(k+1),i+1,y(6*(i)+k)*yscale,yserr(6*(i)+k,1), &
         yserr(6*(i)+k,2),yerr(6*(i)+k)*yscale
   enddo
   j=8+i
   write(nunit,501) titles(9),i+1,solout(j+1),serr(j+1,1),serr(j+1,2),err(j+1)
enddo

close(nunit)

end subroutine exportfit
