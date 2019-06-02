subroutine exportfit(nbodies,sol,serr)
use precision
implicit none
integer :: nunit,filestatus,nbodies,np,i,j
real(double), dimension(:) :: sol
real(double), dimension(:,:) :: serr
character(80) :: exportfile
character(3) :: titles(14)
titles(1)='RHO'
titles(2)='NL1'
titles(3)='NL2'
titles(4)='NL3'
titles(5)='NL4'
titles(6)='DIL'
titles(7)="ZPT"
titles(8)="EP"
titles(9)="PE"
titles(10)="BB"
titles(11)="RD"
titles(12)="MA"
titles(13)="EC"
titles(14)="ES"

nunit=10
exportfile="newfit.dat"
open(unit=nunit,file=exportfile,iostat=filestatus)

if(filestatus>0)then
   write(0,*) "Cannot open ",exportfile
   stop
endif

do i=1,7
    if(i.eq.1) sol(i)=abs(sol(i)) !mean stellar density must be positive
    write(nunit,500) titles(i),sol(i),serr(i,1),serr(i,2),0.0d0
enddo
500 format(A3,5(1X,1PE17.10))

do j=1,nbodies
    do i=1,7
        np=7+7*(j-1)
        if((i.eq.2).or.(i.eq.3).or.(i.eq.4).or.(i.eq.5)) sol(np+i)=abs(sol(np+i))
        write(nunit,501) titles(7+i),j,sol(np+i),serr(np+i,1),serr(np+i,2),0.0d0
    enddo
enddo

501 format(A2,I1,5(1X,1PE17.10))

close(nunit)

end subroutine exportfit
