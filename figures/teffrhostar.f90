program teffrhostar
implicit none
integer :: nmax,i,nunit,filestatus,npt,nplot,j,k
real :: dumr,Rsun,G,Msun,R,M,mold
real, parameter :: pi = 3.1415926535897932
real, allocatable, dimension(:) :: lteff,logg,mass,bb,rhostar,lrhostar, &
   x,y,mp,lx,ly
character(80) :: filename,dumc,text

Rsun=6.96342e10 !cm
G   =6.67408e-8 !dyne cm2 g-2
Msun=1.989e33   ! g

do j=1,4
   if (j.eq.1) then
      filename="dartmouth_feh0_0.25Gyr.dat"
      nmax=400
      allocate(lteff(nmax),logg(nmax),mass(nmax))
   elseif(j.eq.2)then
      filename="dartmouth_feh0_1Gyr.dat"
   elseif(j.eq.3)then
      filename="dartmouth_feh0_6Gyr.dat"
   elseif(j.eq.4)then
      filename="dartmouth_feh0_12Gyr.dat"
   endif

   nunit=10
   open(unit=nunit,file=filename,iostat=filestatus,status='old')
   if(filestatus>0)then !trap missing file errors
      write(0,*) "Cannot open ",filename
      stop
   endif

   !first two lines are comments
   read(nunit,*) dumc
   read(nunit,*) dumc

   i=1
   do
      if(i.gt.nmax)then
         write(0,*) "Increase nmax to match data points"
         write(0,*) "nmax: ",nmax,i
         stop
      endif
      read(nunit,*,iostat=filestatus) dumr,mass(i),lteff(i),logg(i)
      if(filestatus == 0) then
      i=i+1
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
   npt=i-1
   write(0,*) "npt: ",npt

   if(j.eq.1) allocate(rhostar(npt),x(npt),y(npt),mp(npt))

   nplot=0
   do i=1,npt
      M=mass(i)*Msun
      R=sqrt(G*M/10**(logg(i)))
      rhostar(i)=M/(4.0/3.0*pi*R*R*R)
      if(rhostar(i).gt.0)then
         nplot=nplot+1
         x(nplot)=lteff(i)
         y(nplot)=log10(rhostar(i))
         mp(nplot)=mass(i)
      endif
!      write(0,*) 10.0**lteff(i),mass(i),R/Rsun,rhostar(i)
   enddo

!set up PGPLOT
!open PGPLOT device
   if(j.eq.1)then
      call pgopen('?')  !'?' lets the user choose the device.
      call PGPAP (8.0,1.0) !use a square 8" across
      call pgpage() !create a fresh page
      call pgslw(3) !thicker lines

      allocate(bb(4),lx(2),ly(2))
      bb(1)=log10(15000.0)!maxval(1.02*x(1:nplot))
      bb(2)=log10(2500.0)!minval(0.98*x(1:nplot))
      bb(3)=log10(200.0)!maxval(1.15*y(1:nplot))
      bb(4)=log10(1.0e-6)!minval(1.05*y(1:nplot))
      call pgwindow(bb(1),bb(2),bb(3),bb(4))

      CALL PGBOX('BCLNTS1',0.0,0,'BCLNTS1',0.0,0)
      call PGSLS(4)
      CALL PGBOX('GL',0.0,0,'GL',0.0,0)
      lx(1)=log10(5000.0)
      lx(2)=log10(5000.0)
      ly(1)=bb(3)
      ly(2)=bb(4)
      call pgline(2,lx,ly)
      call PGSLS(1)
      call pglabel("Teff (K)","rhostar (g/cm\u-3\d)","")
   endif

   call pgsci(j)
   call pgpt(nplot,x,y,17)
   call pgline(nplot,x,y)


   if(j.eq.4)then
      call pgsch(0.6)

      k=1
      write(text,'(F5.2)') mp(k)
      call PGPTXT (X(k)-0.02*(bb(1)-bb(2)),                             &
                   Y(k)-0.01*(bb(3)-bb(4)), 0.0, 0.0, TEXT)
      lx(1)=X(k)-0.02*(bb(1)-bb(2))
      ly(1)=Y(k)-0.01*(bb(3)-bb(4))
      lx(2)=X(k)
      ly(2)=Y(k)
      call pgline(2,lx,ly)

      k=nplot
      write(text,'(F5.2)') mp(k)
      call PGPTXT (X(k)-0.02*(bb(1)-bb(2)),                             &
                   Y(k)-0.01*(bb(3)-bb(4)), 0.0, 0.0, TEXT)
      lx(1)=X(k)-0.02*(bb(1)-bb(2))
      ly(1)=Y(k)-0.01*(bb(3)-bb(4))
      lx(2)=X(k)
      ly(2)=Y(k)
      call pgline(2,lx,ly)

      mold=mp(1)
      do k=2,nplot
         if(mp(k)-mold.gt.0.2)then
            mold=mp(k)
            write(text,'(F5.2)') mp(k)
            call PGPTXT (X(k)-0.02*(bb(1)-bb(2)),                        &
                         Y(k)-0.01*(bb(3)-bb(4)), 0.0, 0.0, TEXT)
            lx(1)=X(k)-0.02*(bb(1)-bb(2))
            ly(1)=Y(k)-0.01*(bb(3)-bb(4))
            lx(2)=X(k)
            ly(2)=Y(k)
            call pgline(2,lx,ly)
         endif
      enddo
      call pgsch(1.0)
   elseif(j.eq.1)then
      call pgsch(0.6)
      mold=mp(1)
      do k=2,nplot
         if(mp(k)-mold.gt.0.2)then
            mold=mp(k)
            write(text,'(F5.2)') mp(k)
            if(mp(k).lt.2.9)then
               call PGPTXT (X(k)+0.02*(bb(1)-bb(2)),                     &
                            Y(k)+0.04*(bb(3)-bb(4)), 45.0, 1.0, TEXT)
               lx(1)=X(k)+0.02*(bb(1)-bb(2))
               ly(1)=Y(k)+0.04*(bb(3)-bb(4))
               lx(2)=X(k)
               ly(2)=Y(k)
               call pgline(2,lx,ly)
            else
               call PGPTXT (X(k)+0.06*(bb(1)-bb(2)),                     &
                            Y(k)-0.02*(bb(3)-bb(4)), 0.0, 1.0, TEXT)
               lx(1)=X(k)+0.06*(bb(1)-bb(2))
               ly(1)=Y(k)-0.02*(bb(3)-bb(4))
               lx(2)=X(k)
               ly(2)=Y(k)
               call pgline(2,lx,ly)
            endif
         endif
      enddo
      call pgsch(1.0)
   endif

   if(j.eq.1)then
      TEXT="0.25 Gyr"
   elseif(j.eq.2)then
      TEXT="1.00 Gyr"
   elseif(j.eq.3)then
      TEXT="6.00 Gyr"
   else
      TEXT="12.0 Gyr"
   endif

   call PGPTXT(bb(1)-(bb(1)-bb(2))*0.05,bb(4)+(bb(3)-bb(4))*0.04*real(j),  &
      0.0,0.0,TEXT)
!   write(0,*) bb(4)+(bb(3)-bb(4))*0.03*real(j)

enddo

call pgclos()

end program teffrhostar
