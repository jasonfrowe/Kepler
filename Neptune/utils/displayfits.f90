subroutine displayfits(nxmax,nymax,parray,bpix,tavg,refarray)
use precision
implicit none
integer :: nxmax,nymax,npt,i,j,ncol,dumi,ndiff,k,nreplace,it1,it2
integer, dimension(4) :: nr
integer, allocatable, dimension(:) :: p
integer, allocatable, dimension(:,:) :: ia
real :: r,g,b,dumr,x2,y2,xr
real, dimension(4) :: rj
real(double) :: bpix,maxp,minp,z1,z2,med,std,stdev2,tavg,lmin,lmax, &
   ratio,lmed,lstd,minlp,replacethres
real(double), dimension(:,:) :: parray,refarray
real(double), allocatable, dimension(:) :: a
real(double), allocatable, dimension(:,:) :: lparray
character(80) :: tchar

allocate(lparray(nxmax,nymax)) !used for making a log-scale plot

ncol=64 !number of colours for display

nr=1 !range of pixels for plotting
npt=0
do i=1,nxmax
   do j=1,nymax
      if(parray(i,j).lt.bpix)then
         npt=npt+1
         if(npt.eq.1)then
            maxp=parray(i,j)
            minp=parray(i,j)
            nr(1)=i
            nr(2)=i
            nr(3)=j
            nr(4)=j
         else
            maxp=max(maxp,parray(i,j))
            minp=min(minp,parray(i,j))
            nr(1)=min(nr(1),i)
            nr(2)=max(nr(2),i)
            nr(3)=min(nr(3),j)
            nr(4)=max(nr(4),j)
         endif
      endif
   enddo
enddo
!write(0,*) "Min/Max:",minp,maxp

if(npt.le.1)then
   return
endif

!write(0,*) "NR: ",(nr(i),i=1,4)
!ndiff=(nr(2)-nr(1))-(nr(4)-nr(3))
ndiff=9.0*(nr(2)-nr(1))-16.0*(nr(4)-nr(3))
!write(0,*) "ndiff: ",ndiff
if(ndiff.gt.0)then
!   nr(4)=nr(4)+ndiff/2
!   nr(3)=nr(3)-ndiff/2

!   ratio=9.0d0*(nr(2)-nr(1))/(16.0d0*(nr(4)-nr(3)))
!   ratio=1.0/ratio
!   nr(4)=nr(4)+(nr(4)-nr(3))*(1.0d0-ratio)/2.0d0
!   nr(3)=nr(3)-(nr(4)-nr(3))*(1.0d0-ratio)/2.0d0

    it2=(nr(4)+nr(3))/2.0+int(dble(nr(2)-nr(1))*9.0/16.0)/2.0
    it1=(nr(4)+nr(3))/2.0-int(dble(nr(2)-nr(1))*9.0/16.0)/2.0
    nr(4)=it2
    nr(3)=it1

else
!   nr(2)=nr(2)-ndiff/2
!   nr(1)=nr(1)+ndiff/2

!   ratio=16.0*(nr(4)-nr(3))/(9.0d0*(nr(2)-nr(1)))
!   ratio=1.0/ratio
!   it2=nr(2)+(nr(2)-nr(1))*(1.0d0-ratio)*2.0d0
!   it1=nr(1)-(nr(2)-nr(1))*(1.0d0-ratio)*2.0d0
!   nr(2)=it2
!   nr(1)=it1
    it2=(nr(2)+nr(1))/2.0+int(dble(nr(4)-nr(3))*16.0/9.0)/2.0
    it1=(nr(2)+nr(1))/2.0-int(dble(nr(4)-nr(3))*16.0/9.0)/2.0
    nr(2)=it2
    nr(1)=it1
endif
rj=real(nr)

!write(0,*) "hello"

allocate(a(npt),p(npt))
k=0
do i=nr(1),nr(2),5
   do j=nr(3),nr(4),5
      if((parray(i,j).lt.bpix).and.(abs(parray(i,j).lt.50.0)))then
         k=k+1
         a(k)=parray(i,j)
!         write(0,*) k,a(k)
         if(k.ge.1000) goto 10
      endif
   enddo
enddo
10 continue
!write(0,*) "K:",k
if(k.ge.3)then
   call rqsort(k,a,p)
   med=a(p(k/2))
   std=stdev2(k,a,med)
else
   med=0.0
   std=1.0
endif
deallocate(a,p)

replacethres=-10.0
nreplace=0
! uncomment for replacing pixels
!do i=467,478
!   do j=nr(3),nr(4)
!      if(parray(i,j).lt.replacethres)then
!         if(parray(i+1,j).lt.bpix)then
!            parray(i,j)=parray(i+1,j)
!            nreplace=nreplace+1
!         else
!            parray(i,j)=0.0d0
!         endif
!      endif
!      if(parray(i,j).lt.replacethres)then
!         parray(i,j)=refarray(i,j)
!         nreplace=nreplace+1
!      endif
!   enddo
!enddo
!write(0,*) "Nreplace: ",nreplace

!read(5,*)

allocate(ia(nxmax,nymax))

minlp=-70.0
lmin=1000.0
lmax=-1000.0
do i=nr(1),nr(2)
   do j=nr(3),nr(4)
      if(parray(i,j).lt.bpix)then
         if(parray(i,j)-minlp+1.0.le.0.0d0)then
            lparray(i,j)=0.0
         else
            lparray(i,j)=log10(parray(i,j)-minlp+1.0)
         endif
         lmin=min(lparray(i,j),lmin)
         lmax=max(lparray(i,j),lmax)
      endif
   enddo
enddo

!allocate(a(npt),p(npt))
!k=0
!do i=nr(1),nr(2)
!   do j=nr(3),nr(4)
!      if(parray(i,j).lt.bpix)then
!         k=k+1
!         a(k)=lparray(i,j)
!         if(k.ge.1000) goto 11
!      endif
!   enddo
!enddo
!11 call rqsort(k,a,p)
!lmed=a(p(k/2))
!lstd=stdev2(k,a,med)
!deallocate(a,p)

!z1=log10(minp-minp+0.01)
!z2=log10(maxp-minp+0.01)
!write(0,*) "Z1,Z2:",z1,z2
!write(0,*) "med,std:",med,std
!z1=log10(med-10.0-minp+0.01)
!z2=log10(med+300000.0-minp+0.01)
z1=log10(-15-minlp+1.0)
z2=log10(300000.0-minlp+1.0)
!z1=log10(10.0**lmed-10)
!z2=log10(10.0**lmed+200.0)
!write(0,*) "lmin,lmax:",lmin,lmax
!write(0,*) "Z1,Z2:",z1,z2
do i=nr(1),nr(2)
   do j=nr(3),nr(4)
      if(parray(i,j).lt.bpix)then
!         IA(i,j)=int((lparray(i,j)-z1)/(z2-z1)*(NCOL-1))+16
         IA(i,j)=int((lparray(i,j)-z1)/(z2-z1)*dble(NCOL-1))+16
         if(lparray(i,j).le.z1) then
!            write(0,*) "underflow: ",i,j,lparray(i,j)
            ia(i,j)=16
         endif
      else
         ia(i,j)=15
      endif
      if(lparray(i,j).gt.z2) ia(i,j)=ncol+15
   enddo
enddo
!write(0,*) "IA:",ia(425,738),lparray(425,738)

!set up pgplot window
call pgscr(0,1.0,1.0,1.0)
call pgscr(15,0.0,0.3,0.2)

call pgsch(1.5) !make the font a bit bigger
call pgslw(3)  !make the lines a bit thicker

call pgscr(1,0.0,0.0,0.0)
call pgsci(1)
call pgvport(0.0,1.00,0.0,1.0)
call pgwindow(0.0,1.0,0.0,1.0)
call PGRECT (0.0, 1.0, 0.0, 1.0)

call pgvport(0.10,0.90,0.15,0.95) !make room around the edges for labels
call pgsci(0)
call pgwindow(rj(1),rj(2),rj(3),rj(4)) !plot scale
call pgbox("BCNTS1",0.0,0,"BCNTS1",0.0,0)
call pglabel("X (pixels)","Y (pixels)","")
call pgsci(1)

!open(unit=31,file="/iraf/iraf/unix/sun/heat.lut",status='old')
open(unit=31,file="heat.lut",status='old')
read(31,*) dumi
do i=1,ncol
   read(31,*) r,g,b
   read(31,*) dumr,dumr,dumr
   read(31,*) dumr,dumr,dumr
   read(31,*) dumr,dumr,dumr
   CALL PGSCR(I+15, R, G, B)
!   CALL PGSCR(I+15, sqrt(R), sqrt(G), sqrt(B))
!   CALL PGSCR(I+15, log10(R*9.0+1.0), log10(G*9.0+1.0), log10(B*9.0+1.0))
!   R = REAL(I-1)/REAL(NCOL-1)*0.8 + 0.2
!   G = MAX(0.0, 2.0*REAL(I-1-NCOL/2)/REAL(NCOL-1))
!   B = 0.2 + 0.4*REAL(NCOL-I)/REAL(NCOL)
!   CALL PGSCR(ncol-I+16, R, G, B)
enddo
close(31)

!do i=1,ncol
!   Red = REAL(I-1)/REAL(NCOL-1)*0.8 + 0.2
!   Green = MAX(0.0, 2.0*REAL(I-1-NCOL/2)/REAL(NCOL-1))
!   Blue = 0.2 + 0.4*REAL(NCOL-I)/REAL(NCOL)
!   R = 1.0 - real(i)/real(ncol)
!   G = 1.0 - real(i)/real(ncol)*0.5
!   B = 1.0 - real(i)/real(ncol)*0.3
!   R = real(i)/real(ncol)
!   G = real(i)/real(ncol)*0.5
!   B = real(i)/real(ncol)*0.3
!   CALL PGSCR(I+15, R, G, B)
!enddo

xr=real(max(nr(2)-nr(1),nr(4)-nr(3)))
x2=real(nr(2)-nr(1))/xr
y2=real(nr(4)-nr(3))/xr
!write(0,*) "x2,y2: ",x2,y2
!call pgpixl(ia,nxmax,nymax,nr(1),nr(2),nr(3),nr(4),0.0,x2,0.0,y2)
call pgpixl(ia,nxmax,nymax,nr(1),nr(2),nr(3),nr(4),rj(1),rj(2),rj(3),rj(4))
call pgsci(0)
call pgbox("BCNTS1",0.0,0,"BCNTS1",0.0,0)

!call pgsci(2)
!call pgsfs(2)
!call pgcirc(210.0,775.0,3.0)
!call pgsfs(1)
!call pgsci(1)

write(tchar,500) tavg
500 format(F4.1)
call PGPTXT(rj(2)-0.10*(rj(2)-rj(1)), rj(4)-0.08*(rj(4)-rj(3)), 0.0, 0.0, tchar)
call pgsch(0.9)
call PGPTXT(rj(2)-0.12*(rj(2)-rj(1)), rj(3)-0.12*(rj(4)-rj(3)), 0.0, 0.0, "Kepler/K2 Team")
call PGPTXT(rj(2)-0.12*(rj(2)-rj(1)), rj(3)-0.15*(rj(4)-rj(3)), 0.0, 0.0, "J.Rowe (UdeM/iREx)")
call pgsch(1.5)
call pgsci(1)

deallocate(ia,lparray)

!update refarray

!do i=nr(1),nr(2)
!   do j=nr(3),nr(4)
!      if(parray(i,j).gt.replacethres)then
!         refarray(i,j)=parray(i,j)
!      endif
!   enddo
!enddo


write(11,501) tavg,minp,maxp,med,std,z1,z2
501 format(7(1PE17.10,1X))

return
end subroutine displayfits
