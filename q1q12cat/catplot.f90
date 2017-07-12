program catplot
use precision
implicit none
integer :: iargc,nmax,nrec,nunit,filestatus,i,j,npt,fam,num,bmax,nplot, &
nPC,nFP
integer, allocatable, dimension(:) :: id,ils
real :: Rs,Gc,Pi,dF,bb,Per
real, allocatable, dimension(:) :: rj,xp,yp
real, allocatable, dimension(:,:) :: records
character(80) catfile,xlabel,ylabel

interface
   subroutine plotter(npt,records,rj,id,ils,xlabel,ylabel,ntype,nfpplot)
      use precision
      integer, intent(in) :: npt,ntype,nfpplot
      integer, dimension(:), intent(in) :: id,ils
      real, dimension(:), intent(in) :: rj
      real, dimension(:,:), intent(in) :: records
      character(80) :: xlabel,ylabel
   end subroutine plotter
end interface

if(iargc().lt.1) then
   write(0,*) "Usage: catplot Q1Q12catalogue"
   stop
endif

call getarg(1,catfile)

nunit=10
open(unit=nunit,file=catfile,iostat=filestatus,status='old')
if(filestatus>0)then
   write(0,*) "Cannot open ",catfile
   stop
endif

nmax=8000
nrec=37
allocate(records(nmax,nrec+3)) !add three extra records to allow custom variables

nPC=0
nFP=0
i=1
do
   if(i.gt.nmax)then !ran out of array space
      write(0,*) "not enough array space Increase nmax"
      stop
   endif
   read(nunit,*,iostat=filestatus) (records(i,j),j=1,nrec)
   if(filestatus == 0) then
!      if(records(i,29).gt.10.0)
      if(records(i,35).lt.7.1) records(i,37)=1.0d0
!      if(records(i,29).gt.20.0) records(i,37)=1.0d0
      if(records(i,37).lt.0.5)then
         nPC=nPC+1
      else
         nFP=nFP+1
      endif
      i=i+1
      cycle
   elseif(filestatus == -1) then  !EOF
      exit
   else
      write(0,*) records(i-1,2)
      write(0,*) "File Error! line:",i
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
npt=i-1
write(0,*) "Records read: ",npt
write(0,*) "Number of PCs: ",nPC
write(0,*) "Number of FPs: ",nFP
close(nunit)

call pgopen('?')
call pgpage()
call PGPAP ( 8.0 ,1.0)
call pgsch(1.2)
call pgslw(2)
call pgvport(0.15,0.85,0.15,0.85)

!1=with err,2 withouterr
!-1 all, 1 PCs, 0 FPs

allocate(id(6),rj(4),ils(2))

ils(1)=1 !logscale
ils(2)=1 !logscale

id(1)=32 !Srad
id(2)=33
id(3)=34
id(4)=29 !Rp
id(5)=30
id(6)=31

rj(1)=0.1
rj(2)=100000.0
rj(2)=10.0
rj(3)=0.2
rj(3)=0.4
rj(4)=100.0
rj(4)=20.0
rj=log10(rj)

xlabel="S (S\d\(2284)\u)"
ylabel="Radius (R\d\(2284)\u)"

call plotter(npt,records,rj,id,ils,xlabel,ylabel,1,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,1,0)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,1)
!goto 901

!ccccc
id(1)=17 !Period
id(2)=18
id(3)=18
id(4)=29 !Rp
id(5)=30
id(6)=31

rj(1)=0.1
rj(2)=1000.0
rj(3)=0.2
rj(4)=100.0
rj=log10(rj)

xlabel="Period (days)"
ylabel="Radius (R\d\(2284)\u)"

call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,1,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,0)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,1)


!ccccc
id(1)=27 !Tdur
id(2)=28
id(3)=28
id(4)=25 !Tdep
id(5)=26
id(6)=26

rj(1)=0.5
rj(2)=40.0
rj(3)=10.0
rj(4)=100000.0
rj=log10(rj)

xlabel="T\ddur\u (hours)"
ylabel="T\ddep\u (ppm)"

call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,1,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,0)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,1)

!ccccc
id(1)=27 !Tdur
id(2)=28
id(3)=28
id(4)=12 !rhostar
id(5)=13
id(6)=14

rj(1)=0.5
rj(2)=40.0
rj(3)=0.001
rj(4)=100.0
rj=log10(rj)

xlabel="T\ddur\u (hours)"
ylabel="\gr\dc\u (g/cm\u3\d)"

call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,1,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,-1)
nplot=1000
allocate(xp(nplot),yp(nplot))
Pi=acos(-1.e0)!define Pi
Rs=6.9599e+10
Gc=6.67259e-8
dF=84.0e-6 !transit depths
bb=0.78 !impact
Per=370.0*24.0*60.0*60.0
do i=1,nplot
   xp(i)=10**rj(1)+(10**rj(2)-10**rj(1))/real(nplot)*real(i)
   yp(i)=(4.0*pi*pi/(Per**2.0*Gc))*  &
    (( (1.0+sqrt(dF))**2.0 - bb*bb*(1-(sin(xp(i)*60.0*60.0*pi/Per))**2.0)) / &
    (sin(xp(i)*60.0*60.0*Pi/Per))**2.0 )**(3.0/2.0)
!   write(0,*) xp(i),yp(i)
!   read(5,*)
enddo
call pgsci(3)
xp=log10(xp)
yp=log10(yp)
call pgline(nplot,xp,yp)
call pgsci(1)

call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,0)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,1)

!ccccc
id(1)=19 !b
id(2)=20
id(3)=21
id(4)=22 !RpRs
id(5)=23
id(6)=24

rj(1)=0.0
rj(2)=2.0
rj(3)=log10(0.001)
rj(4)=log10(10.0)
ils(1)=0 !linearscale
ils(2)=1 !logscale

xlabel="b"
ylabel="Radius (Rp/Rs)"

call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,1,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,0)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,1)

!ccccc
id(1)=3 !Teff
id(2)=4
id(3)=4
id(4)=12 !rhostar
id(5)=13
id(6)=14

rj(2)=log10(3000.0)
rj(1)=log10(10000.0)
rj(4)=log10(0.001)
rj(3)=log10(100.0)
ils(1)=1 !logscale
ils(2)=1 !logscale

xlabel="Teff"
ylabel="\gr\dc\u (g/cm\u3\d)"

call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,1,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,0)
call pgsci(2)
call isoplot(1)
call isoplot(2)
call pgsci(3)
call isoplot(3)
call isoplot(4)
call pgsci(4)
call isoplot(5)
call isoplot(6)
call pgsci(1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,1)

!ccccc
Pi=acos(-1.e0)!define Pi
Rs=4.0/3.0*pi*6.9599e+10
Gc=6.67259e-8
do i=1,npt
!   write(0,*) 10**records(i,5)/(records(i,9)*Rs*Gc)
   records(i,38)=10**records(i,5)/(records(i,9)*Rs*Gc)
   records(i,39)=10**(records(i,5)+records(i,6))/ &
      ((records(i,9)+records(i,11))*Rs*Gc)
   records(i,40)=10**(records(i,5)-records(i,6))/ &
      ((records(i,9)+records(i,10))*Rs*Gc)
!   write(0,*) records(i,38),records(i,39),records(i,40)
!   read(5,*)
enddo

id(1)=38 !rhostar-KIC
id(2)=39
id(3)=40
id(4)=12 !rhostar
id(5)=13
id(6)=14

rj(1)=log10(0.1)
rj(2)=log10(30.0)
rj(3)=log10(0.1)
rj(4)=log10(30.0)
ils(1)=1 !logscale
ils(2)=1 !logscale

xlabel="\gr\dKIC\u (g/cm\u3\d)"
ylabel="\gr\dc\u (g/cm\u3\d)"

call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,1,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,0)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,1)

!ccccc
id(1)=3 !Teff
id(2)=4
id(3)=4
id(4)=27 !Tdur
id(5)=28
id(6)=28

rj(1)=log10(3000.0)
rj(2)=log10(10000.0)
rj(3)=log10(0.5)
rj(4)=log10(40.0)
ils(1)=1 !logscale
ils(2)=1 !logscale

xlabel="Teff"
ylabel="T\ddur\u (hours)"

call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,1,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,0)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,1)

!ccccc
id(1)=17 !Period
id(2)=18
id(3)=18
id(4)=27 !Tdur
id(5)=28
id(6)=28

rj(1)=log10(0.1)
rj(2)=log10(1000.0)
rj(3)=log10(0.5)
rj(4)=log10(40.0)
ils(1)=1 !logscale
ils(2)=1 !logscale

xlabel="Period (days)"
ylabel="T\ddur\u (hours)"

call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,1,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,0)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,1)

!ccccc
id(1)=3 !Teff
id(2)=4
id(3)=4
id(4)=29 !Rp
id(5)=30
id(6)=31

rj(1)=log10(3000.0)
rj(2)=log10(7000.0)
rj(3)=log10(0.2)
rj(4)=log10(100.0)
ils(1)=1 !logscale
ils(2)=1 !logscale

xlabel="Teff (K)"
ylabel="Radius (R\d\(2284)\u)"

call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,1,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,0)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,1)

!ccccc
id(1)=9 !Teff
id(2)=10
id(3)=11
id(4)=29 !Rp
id(5)=30
id(6)=31

rj(1)=log10(0.1)
rj(2)=log10(100.0)
rj(3)=log10(0.2)
rj(4)=log10(100.0)
ils(1)=1 !logscale
ils(2)=1 !logscale

xlabel="Rstar"
ylabel="Radius (R\d\(2284)\u)"

call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,1,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,0)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,1)

!ccccc
Pi=acos(-1.e0)!define Pi
Rs=4.0/3.0*pi*6.9599e+10
Gc=6.67259e-8
do i=1,npt
!   write(0,*) 10**records(i,5)/(records(i,9)*Rs*Gc)
   records(i,38)=10**records(i,5)/(records(i,9)*Rs*Gc)
   records(i,39)=10**(records(i,5)+records(i,6))/ &
      ((records(i,9)+records(i,11))*Rs*Gc)
   records(i,40)=10**(records(i,5)-records(i,6))/ &
      ((records(i,9)+records(i,10))*Rs*Gc)
!   write(0,*) records(i,38),records(i,39),records(i,40)
!   read(5,*)
enddo

id(1)=3 !Teff
id(2)=4
id(3)=4
id(4)=38 !rhostar_KIC
id(5)=39
id(6)=40

rj(2)=log10(3000.0)
rj(1)=log10(10000.0)
rj(4)=log10(0.001)
rj(3)=log10(100.0)
ils(1)=1 !logscale
ils(2)=1 !logscale

xlabel="Teff"
ylabel="\gr\dc\u (g/cm\u3\d)"

call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,1,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,-1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,0)
call pgsci(2)
call isoplot(1)
call isoplot(2)
call pgsci(3)
call isoplot(3)
call isoplot(4)
call pgsci(4)
call isoplot(5)
call isoplot(6)
call pgsci(1)
call pgpage()
call plotter(npt,records,rj,id,ils,xlabel,ylabel,2,1)

901 call pgclos()

end program catplot

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine isoplot(ntype)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
integer nunit,nisomax,niso,i,filestatus,dumi,ntype
real :: mass,logg,G,Msun,Rsun,Pi,fourthreepi
real, allocatable, dimension(:) :: lrhostar,lTeff
character(80) :: isofile

Pi=acos(-1.d0)   !Pi
fourthreepi=4.0/3.0*Pi
Msun=1.9891d30 !g  mass of Sun
Rsun=696265.0d0*1000.0d0 !cm  radius of Sun
G=6.674d-11 !N m^2 kg^-2  Gravitation constamt

select case(ntype)
   case(1)
      isofile="isochrone.1GYR.dat"
   case(2)
      isofile="isochrone.14GYR.dat"
   case(3)
      isofile="isochrone.1GYRp5.dat"
   case(4)
      isofile="isochrone.14GYRp5.dat"
   case(5)
      isofile="isochrone.1GYRm2.dat"
   case(6)
      isofile="isochrone.14GYRm2.dat"
end select

nunit=10
open(unit=nunit,file=isofile,iostat=filestatus,status='old')
if(filestatus>0)then
   write(0,*) "Cannot open ",catfile
   stop
endif

nisomax=500
allocate(lrhostar(nisomax),lTeff(nisomax))

i=1
do
   if(i.gt.nisomax)then !ran out of array space
      write(0,*) "not enough array space Increase nisomax"
      stop
   endif
   read(nunit,*,iostat=filestatus) dumi,mass,lTeff(i),logg
   if(filestatus == 0) then
      R=sqrt(G*mass*Msun/(10**logg))*10.0
      lrhostar(i)=log10(10.0**logg/(fourthreepi*G*R)/100.0/1000.0) !/g/cc
!      write(0,*) lTeff(i),rhostar(i),R/Rsun
      i=i+1
      cycle
   elseif(filestatus == -1) then  !EOF
      exit
   else
      write(0,*) "File Error! line:",i
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
niso=i-1
write(0,*) "Records read: ",niso
close(nunit)

call pgline(niso,lTeff,lrhostar)

return
end subroutine isoplot

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine plotter(npt,records,rj,id,ils,xlabel,ylabel,ntype,nfpplot)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: i,j,nplot,npt,ntype,nfpplot
integer, dimension(:) :: id,ils
integer, allocatable, dimension(:) :: icol,isym,ifp
real :: eps
real, dimension(:) :: rj
real, allocatable, dimension(:) :: xp,xm,yp,ym,x,y
real, dimension(:,:) :: records
character(80) :: xlabel,ylabel,cboxx,cboxy

allocate(xp(npt),xm(npt),yp(npt),ym(npt),x(npt),y(npt),icol(npt), &
 isym(npt),ifp(npt))

eps=1.0e-10
j=1
do i=1,npt

   x(j)=records(i,id(1))
   xp(j)=records(i,id(1))+records(i,id(2))
   xm(j)=records(i,id(1))-abs(records(i,id(3)))
   if(ils(1).eq.1)then
      x(j)=log10(x(j)+eps)
      xp(j)=log10(xp(j)+eps)
      xm(j)=log10(xm(j)+eps)
   endif
   y(j)=records(i,id(4))
   yp(j)=records(i,id(4))+records(i,id(5))
   ym(j)=records(i,id(4))-abs(records(i,id(6)))
   if(ils(2).eq.1)then
      y(j)=log10(y(j)+eps)
      yp(j)=log10(yp(j)+eps)
      ym(j)=log10(ym(j)+eps)
   endif

   ifp(j)=int(records(i,37)+0.5)
   if(records(i,37).lt.0.5)then
      isym(j)=4
      icol(j)=1
   else
      isym(j)=6
      icol(j)=2
   endif

!   if((records(i,17).gt.300.0).and.(records(i,17).lt.500.0).and.(records(i,29).lt.20.0))then
!      icol(j)=1
!      isym(j)=5
!   elseif(records(i,29).lt.10.0)then
!      icol(j)=2
!      isym(j)=6
!   else
!      icol(j)=4
!      isym(j)=7
!   endif

!   write(0,*) x(j),xp(j),xm(j)
!   read(5,*)
   if( (x(j).ge.min(rj(1),rj(2))) .and. (x(j).le.max(rj(1),rj(2))) .and.&
    (y(j).ge.min(rj(3),rj(4))) .and. (y(j).le.max(rj(3),rj(4))) .and.   &
    (xp(j).gt.x(j)) .and. (xm(j).lt.x(j)) &
    .and. (yp(j).gt.y(j)) .and. (ym(j).lt.y(j))  ) then
      j=j+1
   endif
enddo
nplot=j-1

call pgwindow(rj(1),rj(2),rj(3),rj(4))

if(ils(1).eq.1)then
   cboxx='BCLNTS1'
else
   cboxx='BCNTS1'
endif
if(ils(2).eq.1)then
   cboxy='BCLNTS1'
else
   cboxy='BCNTS1'
endif
call pgbox(cboxx,0.0,0,cboxy,0.0,0)
call pglabel(xlabel,ylabel,"")
call pgsch(0.8)
call pgbbuf()
do i=1,nplot

if(nfpplot.eq.0)then
   if(ntype.eq.1)then
      if(ifp(i).eq.0)then
         call pgsci(icol(i))
         call pgpt(1,x(i),y(i),isym(i))
         call pgerrx(1,xm(i),xp(i),y(i),1.0)
         call pgerry(1,x(i),ym(i),yp(i),1.0)
      endif
   elseif(ntype.eq.2)then
      if(ifp(i).eq.0)then
         call pgsci(icol(i))
         call pgpt(1,x(i),y(i),isym(i))
      endif
   endif
elseif(nfpplot.eq.1)then
   if(ntype.eq.1)then
      if(ifp(i).eq.1)then
         call pgsci(icol(i))
         call pgpt(1,x(i),y(i),isym(i))
         call pgerrx(1,xm(i),xp(i),y(i),1.0)
         call pgerry(1,x(i),ym(i),yp(i),1.0)
      endif
   elseif(ntype.eq.2)then
      if(ifp(i).eq.1)then
         call pgsci(icol(i))
         call pgpt(1,x(i),y(i),isym(i))
      endif
   endif
else
   if(ntype.eq.1)then
      call pgsci(icol(i))
      call pgpt(1,x(i),y(i),isym(i))
      call pgerrx(1,xm(i),xp(i),y(i),1.0)
      call pgerry(1,x(i),ym(i),yp(i),1.0)
   elseif(ntype.eq.2)then
      call pgsci(icol(i))
      call pgpt(1,x(i),y(i),isym(i))
   endif
endif

enddo
call pgsci(1)
call pgsch(1.2)
call pgebuf()

return
end
