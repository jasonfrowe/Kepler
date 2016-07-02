subroutine fitterv2(naxes,Image,nstar,Ic,xcoo,ycoo,ngsol,gsol,sat,bpix,   &
 starmap,numnei,nei,id)
use precision
implicit none
!import vars
integer :: nstar,ngsol
integer, dimension(2) :: naxes
integer, dimension(:) :: numnei,id
integer, dimension(:,:) :: starmap,nei
real(double) :: bpix
real(double) :: sat
real(double), dimension(:) :: Ic,xcoo,ycoo,gsol
real(double), dimension(:,:) :: Image
!local vars
integer :: npix,i,j,nfit,n,m,iprint,isave(44)
integer, allocatable, dimension(:) :: nbd,iwa
real(double) :: pgtol,factr,f,dsave(29)
real(double), allocatable, dimension(:) :: sol,l,u,g,wa
logical :: lsave(4)
character(60) :: task,csave

interface
   function chisq(naxes,Image,sat,npix,ngsol,sol,starmap,numnei,nei)
      use precision
      implicit none
      integer :: npix,ngsol
      integer, dimension(2) :: naxes
      integer, dimension(:) :: numnei
      integer, dimension(:,:) :: starmap,nei
      real(double) :: chisq,sat
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: Image
   end function chisq
end interface
interface
   subroutine gradient(nfit,f,g,naxes,Image,sat,npix,ngsol,sol,starmap, &
    numnei,nei)
      use precision
      implicit none
      integer :: npix,ngsol,nfit
      integer, dimension(2) :: naxes
      integer, dimension(:) :: numnei
      integer, dimension(:,:) :: starmap,nei
      real(double) :: chisq,sat,f
      real(double), dimension(:) :: sol,g
      real(double), dimension(:,:) :: Image
   end subroutine gradient
end interface

npix=0
do i=1,naxes(1)
   do j=1,naxes(2)
      if(Image(i,j).lt.sat)then
         npix=npix+1
      endif
   enddo
enddo
write(0,*) "Number of pixels to fit in Image: ",npix

!number of parameters to fit
nfit=ngsol+3*nstar
allocate(sol(nfit))

!initial Guesses
sol(1)=0.0
sol(2)=1.1
sol(3)=1.1
sol(4)=0.0
do i=1,nstar
   sol(ngsol+3*(i-1)+1)=Image(int(xcoo(i)),int(ycoo(i)))*0.6
   sol(ngsol+3*(i-1)+2)=xcoo(i)
   sol(ngsol+3*(i-1)+3)=ycoo(i)
enddo

n=nfit !number of variables
m=5 !corrections used in limited memory matrix
allocate(l(nfit),u(nfit),nbd(nfit))
nbd(1)=0 !no bound
nbd(2)=1 !lower bound
l(2)=1.0e-8
nbd(3)=1 !lower bound
l(3)=1.0e-8
nbd(4)=0 !no bound
do i=1,nstar
   nbd(ngsol+3*(i-1)+1)=1
     l(ngsol+3*(i-1)+1)=1.0e-8
   nbd(ngsol+3*(i-1)+2)=0
   nbd(ngsol+3*(i-1)+3)=0
enddo

allocate(g(nfit))
factr=1.0d+7
pgtol=1.0d-5
allocate ( wa(2*m*n + 5*n + 11*m*m + 8*m) )
allocate ( iwa(3*n) )
iprint=-1  !diagonistic info to print to screen (set negative to quiet)

task = 'START'

do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
               task.eq.'START')

   call setulb ( n, m, sol, l, u, nbd, f, g, factr, pgtol, &
                       wa, iwa, task, iprint,&
                       csave, lsave, isave, dsave )

   write(0,'(A6,A6)') "task: ",task

   if (task(1:2) .eq. 'FG') then

      f=chisq(naxes,Image,sat,npix,ngsol,sol,starmap,numnei,nei)

      write(0,*) "Starting gradient"
      call gradient(nfit,f,g,naxes,Image,sat,npix,ngsol,sol,starmap,numnei,nei)

      write(0,*) "f_r: ",f/npix
      write(0,*) "sol: ",(sol(i),i=1,7)
   endif

   call exportfit(nstar,nfit,sol,ngsol,id)

!   write(0,*) "Pause..."
!   read(5,*)

enddo

return
end subroutine fitterv2

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine exportfit(nstar,nfit,sol,ngsol,id)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: nfit,nstar,ngsol
integer, dimension(nstar) :: id
real(double), dimension(nfit) :: sol
!local vars
integer i,j

open(unit=10,file="solution.txt")
write(10,500) (sol(j),j=1,4)
500 format(4(1PE17.10,1X))
do i=1,nstar
   write(10,501) id(i),sol(ngsol+3*(i-1)+1),sol(ngsol+3*(i-1)+2),       &
    sol(ngsol+3*(i-1)+3)
   501 format(I9,1X,3(1PE17.10,1X))
enddo
close(10)

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine gradient(nfit,f,g,naxes,Image,sat,npix,ngsol,sol,starmap,numnei,nei)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: npix,ngsol,nfit
integer, dimension(2) :: naxes
integer, dimension(:) :: numnei
integer, dimension(:,:) :: starmap,nei
real(double) :: sat,f
real(double), dimension(:) :: sol,g
real(double), dimension(:,:) :: Image
!local vars
integer :: i
real(double) ftest,eps
real(double), allocatable, dimension(:) :: sol2

interface
   function chisq(naxes,Image,sat,npix,ngsol,sol,starmap,numnei,nei)
      use precision
      implicit none
      integer :: npix,ngsol
      integer, dimension(2) :: naxes
      integer, dimension(:) :: numnei
      integer, dimension(:,:) :: starmap,nei
      real(double) :: sat,chisq
      real(double), dimension(:) :: sol
      real(double), dimension(:,:) :: Image
   end function chisq
end interface

allocate(sol2(nfit))

eps=1.0d-8
do i=1,nfit
   sol2=sol
   sol2(i)=sol2(i)+eps
   ftest=chisq(naxes,Image,sat,npix,ngsol,sol2,starmap,numnei,nei)
   g(i)=(ftest-f)/eps
!   write(0,*) i,g(i)
enddo


return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function chisq(naxes,Image,sat,npix,ngsol,sol,starmap,numnei,nei)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: npix,ngsol
integer, dimension(2) :: naxes
integer, dimension(:) :: numnei
integer, dimension(:,:) :: starmap,nei
real(double) :: chisq,sat
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: Image
!local vars
integer npix2,i,j,k,ii,kk
real(double) :: di,dj,model,gpixmod
real(double), dimension(7) :: lsol

lsol(1)=sol(1)
lsol(2)=sol(2)
lsol(3)=sol(3)
lsol(4)=sol(4)

chisq=0.0d0
npix2=0
do i=1,naxes(1)
   di=dble(i)
   do j=1,naxes(2)
      dj=dble(j)
      if(Image(i,j).lt.sat)then
         npix2=npix2+1
         if(npix2.gt.npix)then
            write(0,*) "Npix2 greater than npix.. "
            write(0,*) "number of pixels not counted correctly"
            stop
         endif
         k=starmap(i,j)
         lsol(5)=sol(ngsol+3*(k-1)+1) !amplitude
         lsol(6)=sol(ngsol+3*(k-1)+2)-di !dx
         lsol(7)=sol(ngsol+3*(k-1)+3)-dj !dy
         model=gpixmod(lsol)
         do ii=1,numnei(k) !deal with neighbours
            kk=nei(k,ii)
            lsol(5)=sol(ngsol+3*(kk-1)+1) !amplitude
            lsol(6)=sol(ngsol+3*(kk-1)+2)-di !dx
            lsol(7)=sol(ngsol+3*(kk-1)+3)-dj !dy
            model=model+gpixmod(lsol)
         enddo
         chisq=chisq+(model-Image(i,j))**2.0/abs(Image(i,j))
!         write(0,*) model,Image(i,j),chisq
!         read(5,*)
      endif
   enddo
enddo

return
end function chisq

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function gpixmod(sol)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
real(double) :: gpixmod
real(double), dimension(3) :: temp
real(double), dimension(7) :: sol

temp(1)=-sol(6)*sol(6)/(sol(2)*sol(2))
temp(2)=-sol(7)*sol(7)/(sol(3)*sol(3))
temp(3)=2.0d0*sol(4)*sol(6)*sol(7)/(sol(2)*sol(3))

gpixmod=sol(1)+sol(5)*exp(temp(1)+temp(2)+temp(3))

return
end function gpixmod
