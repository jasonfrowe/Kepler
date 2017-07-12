program hzfigure1
use precision
implicit none
integer :: subx,suby,ix,iy,nkoi,i,j,k,nmax,npt,nunit,nfit,nplanet,ii,   &
   col,nplanetmax,npt2,bins,nplot
integer, allocatable, dimension(:) :: koi,kid,ntt,dtype,flag,dtype2
integer, allocatable, dimension(:,:) :: pplot
real, dimension(4) :: rbb
real, allocatable, dimension(:) :: px,py,pz
real(double) :: ztime,period,toff,ph1,ph2,tdurmax,xwidth,ratio,t1,t2
real(double), allocatable, dimension(:) :: time,flux,ferr,itime,sol,err,&
   sol2,tmodel,phase,phase2,flux2,ferr2,time2,itime2
real(double), allocatable, dimension(:,:) :: serr,tobs,omc
character(80) :: photfile,n0file,label

interface
   subroutine readdata2(obsfile,nmax,npt,time,flux,ferr,itime,ztime)
      use precision
      implicit none
      integer :: nmax,npt
      real(double) :: ztime
      real(double), dimension(:) :: time,flux,ferr,itime
      character(80) :: obsfile
   end subroutine readdata2
end interface

nkoi=18 !number of KOIs
allocate(koi(nkoi),kid(nkoi),pplot(nkoi,9))
pplot=0 !default is to not plot. (1=plot planet transit)

koi(1)=172
kid(1)=8692861
pplot(1,2)=1

koi(2)=438
kid(2)=12302530
pplot(2,2)=1

koi(3)=463
kid(3)=8845205
pplot(3,1)=1

koi(4)=812
kid(4)=4139816
pplot(4,3)=1

koi(5)=854
kid(5)=6435936
pplot(5,1)=1

koi(6)=2418
kid(6)=10027247
pplot(6,1)=1

koi(7)=2626
kid(7)=11768142
pplot(7,1)=1

koi(8)=2650
kid(8)=8890150
pplot(8,1)=1

koi(9)=3010
kid(9)=3642335
pplot(9,1)=1

koi(10)=3282
kid(10)=12066569
pplot(10,1)=1

koi(11)=3497
kid(11)=8424002
pplot(11,1)=1

koi(12)=4036
kid(12)=11415243
pplot(12,1)=1

koi(13)=4054
kid(13)=6428794
pplot(13,1)=1

koi(14)=4356
kid(14)=8459663
pplot(14,1)=1

koi(15)=4450
kid(15)=7429240
pplot(15,1)=1

koi(16)=4550
kid(16)=5977470
pplot(16,1)=1

koi(17)=5236
kid(17)=6067545
pplot(17,1)=1

koi(18)=5856
kid(18)=11037818
pplot(18,1)=1

!koi(19)=7235
!kid(19)=9821428
!pplot(19,1)=1

!arrays for storing photometry
nmax=2000000 !maximum number of observations that can be read in
allocate(time(nmax),flux(nmax),ferr(nmax),itime(nmax),dtype(nmax),      &
   phase(nmax),phase2(nmax),flux2(nmax),ferr2(nmax),flag(nmax))
allocate(px(nmax),py(nmax),pz(nmax)) !arrays for plotting.

!arrays for holding transit solution
nfit=108
allocate(sol(nfit),serr(nfit,2),err(nfit),sol2(nfit))
nunit=10 !unit number used for file input

subx=6 !number of sub-plots
suby=3
ix=2   !index of sub-plot to start at.
iy=2

!Set global scale for all the transit models.
xwidth=15.0d0 !range in x-axes in hours
tdurmax=xwidth/24.0d0
rbb(1)=-real(xwidth)
rbb(2)=+real(xwidth)
rbb(3)=0.995
rbb(4)=1.0025

nplot=1000 !number of samples for plotting transit model
allocate(time2(nplot),itime2(nplot),dtype2(nplot))

call pgopen('?') !open PGPlot device
!call pgopen('test.png/png')
call pgask(.true.) !don't ask for new page.. just do it.
call pgpap(8.0,0.5) !page size
call pgslw(2) !thicker lines
call pgsubp(subx+2,suby+2)  !break up plot into grid

!loop over all KOIs
do i=1,nkoi

   !Read in photometry
   write(photfile,500) "koi",koi(i),".n/klc",kid(i),".dc.dat"
   500 format(A3,i0,A6,i8.8,A7)
   ztime=54900.0
   call readdata2(photfile,nmax,npt,time,flux,ferr,itime,ztime)
   dtype=0 !set all data to 'photometry' class
   write(0,*) "Number of data points read: ",npt

   !Read in best-fit solution for transit model
   write(n0file,501) "koi",koi(i),".n/n0.scld.dat"
   501 format(A3,i0,A14)
   open(unit=nunit,file=n0file,status='old',err=902)
!   write(0,*) "reading in input solution"
   call getfitpars(nunit,nfit,nplanet,sol,serr,err)
!   write(0,*) "done reading input solution"
   close(nunit) !release unit number as we are done with file
   goto 903 !goofy goto to use F77 code
   902 write(0,*) "Cannot open ",n0file
   stop
   903 continue

!   write(0,*) "Nplanet: ",nplanet
   allocate(ntt(nplanet),tobs(nplanet,npt),omc(nplanet,npt))
   ntt=0 !ignoring TTVs for this exercise. (it's only a paper figure)

   !array to hold transit model
   allocate(tmodel(npt))

   !for each KOI, loop over all planets in the system
   do ii=1,nplanet
      if(pplot(i,ii).eq.1)then !only a subset of planets are plottedd

         !calculate offset to pick off parameters of individual planets
         col=10*(ii-1)

         !write(6,*) "ix,iy: ",ix,iy
         !Set up window for plotting (size, scale)
         call pgpanl(ix,iy)
         call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))
         call pgvport(0.00,1.00,0.0,1.0)


!        Create model with just 1 planet
         do j=1,nfit
            sol2(j)=sol(j)
         enddo
         sol2(12+col)=0.0d0 !set r/R*=0
         period=sol2(10+col)
         nplanetmax=nplanet
         call transitmodel(nfit,nplanet,nplanetmax,sol2,nmax,npt,time,  &
            itime,ntt,tobs,omc,tmodel,dtype)

!        Determine range of phase that will be plotted.
         toff=0.5-(sol(9+col)/sol(10+col)-int(sol(9+col)/sol(10+col)))
         if(toff.lt.0.0)toff=toff+1.0
         ph1=0.5-1.0d0*tdurmax/sol(10+col)
         if(ph1.lt.0.25)ph1=0.25
         ph2=0.5+1.0d0*tdurmax/sol(10+col)
         if(ph2.gt.0.75)ph2=0.75

         !range for plotting data.
         rbb(1)=real(ph1)
         rbb(2)=real(ph2)

         !phase data to orbital period and offset.
         ratio=1.0d0
         call phasept(npt,time,phase,period,toff,ratio)

         !get data points that will be plotted
         k=0
         do j=1,npt
            if((phase(j).gt.ph1).and.(phase(j).lt.ph2))then
               k=k+1
               px(k)=real(phase(j))
               py(k)=real(flux(j)-tmodel(j))+1.0
               !write(0,*) px(k),py(k)
            endif
         enddo

         !plot the observations.
         !write(6,*) "rbb:",rbb
         call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))
         call pgpt(k,px,py,17)

         !bin and plot observations
         k=0
         do j=1,npt
            if((real(flux(j)).gt.rbb(3)).and.(real(flux(j)).lt.rbb(4)))then
               k=k+1
               phase2(k)=phase(j)
               flux2(k)=flux(j)-tmodel(j)+1.0d0
               ferr2(k)=ferr(j)
            endif
         enddo
         npt2=k

         !use 30 minute bins
         bins=int(sol(col+8+2)*1440.0d0/30.0d0)
!         write(6,*) "bc: ",sol(col+8+2),bins,col+8+2
         call binp(npt2,phase2,flux2,ferr2,bins,flag)

         !get data that is in plotting range, convert to real and plot.
         j=0
         do k=1,bins
            if(flag(k).eq.0)then
               if((phase2(k).gt.ph1).and.(phase2(k).lt.ph2))then
                  j=j+1
                  px(j)=real(phase2(k))
                  py(j)=real(flux2(k))
                  pz(j)=real(ferr2(k))
               endif
            endif
         enddo
         npt2=j
         call pgsci(4)
         call pgpt(npt2,px,py,17)
         call pgerrb(6,npt2,px,py,pz,1.0)
         call pgsci(1)

         !plot transit moddel
         do j=1,8
            sol2(j)=sol(j)
         enddo
         do j=9,18
            sol2(j)=sol(j+col)
         enddo

         t1=sol2(9)-1.5*tdurmax
         t2=sol2(9)+1.5*tdurmax
         do k=1,nplot
            time2(k)=t1+(t2-t1)/dble(nplot)*dble(k-1)
            itime2(k)=itime(1)
            dtype2(k)=0
         enddo
         call transitmodel(nfit,1,nplanetmax,sol2,nmax,nplot,time2,     &
            itime2,ntt,tobs,omc,tmodel,dtype2)

         call phasept(nplot,time2,phase,period,toff,ratio)
         do k=1,nplot
            tmodel(k)=tmodel(k)-sol(8) !correct zero point
            px(k)=real(phase(k))
            py(k)=real(tmodel(k))
         enddo
         call sort2(nplot,px,py)
         call pgsci(2)
         call pgline(nplot,px,py)
         call pgsci(1)

         !add axes
         rbb(1)=-real(xwidth) !scale for labeling axes.
         rbb(2)=+real(xwidth)
         call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))

         !increase font size
         call pgsch(3.0)
         call pgslw(1)

         !add label
         write(label,502) "KOI-",koi(i),".",ii
         502 format(A4,I0,A1,I2.2)
         call pgptxt((rbb(1)+rbb(2))/2.0,rbb(3)+0.1*(rbb(4)-rbb(3)),    &
            0.0,0.5,label)

         if((ix.eq.2).and.(iy.ne.suby+1))then
            call pgbox('BCTS',0.0,0,'BCNTSV1',0.0,0)
         elseif((ix.ne.2).and.(iy.eq.suby+1))then
            call pgbox('BCNTS1',0.0,0,'BCTS',0.0,0)
         elseif((ix.eq.2).and.(iy.eq.suby+1))then
            call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
         else
            call pgbox('BCTS',0.0,0,'BCTS',0.0,0)
         endif

         if((ix.eq.2).and.(iy.eq.3))                                    &
            call pgptxt(rbb(1)-0.28*(rbb(2)-rbb(1)),(rbb(4)+rbb(3))/2.0,&
            90.0,0.5,'Relative Flux')
         if((ix.eq.5).and.(iy.eq.4))                                    &
            call pglabel('Time from mid-transit (hours)','','')

         call pgsch(1.0)
         call pgslw(2)

         !Update panel that we are plotting to.
         ix=ix+1
         if(ix.gt.subx+1)then
            iy=iy+1
            ix=2
         endif
      endif

   enddo
   deallocate(ntt,tobs,omc)
   deallocate(tmodel)

enddo

call pgclos()

end program hzfigure1

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine readdata2(obsfile,nmax,npt,time,flux,ferr,itime,ztime)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: nmax,npt
real(double) :: ztime
real(double), dimension(:) :: time,flux,ferr,itime
character(80) :: obsfile
!local vars
integer :: nunit,filestatus,i
real(double) :: t,f,e,sec2day,it

sec2day=86400.0d0

nunit=10
open(unit=nunit,file=obsfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",obsfile
   stop
endif



i=0
do
   if(i.gt.nmax)then
      write(0,*) "Increase nmax to match data points"
      write(0,*) "nmax: ",nmax
      stop
   endif
   read(nunit,*,iostat=filestatus) t,f,e!,it
   it=0.0
   if(filestatus == 0) then
      i=i+1
      time(i)=t-ztime+0.5d0
      flux(i)=f+1.0
      ferr(i)=e
      if (it.lt.0.5) then
         itime(i)=1765.5/sec2day
      else
         itime(i)=58.85/sec2day !short cadence
      endif
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!! Line:",i+1
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit) !close file
npt=i

return
end subroutine readdata2


