      program histanim
      implicit none
      integer nmax,i,nunit,npt,j,jj,nframe
      parameter(nmax=27)
      real year(nmax),nother(nmax),nkep(nmax),nmulti(nmax),dumr,mbin,
     .   xrange
      character*80 filename
      character*4 cyear

      filename="planethistogram.dat"

      nunit=10
      open(unit=nunit,file=filename,status='old',err=901)
      i=1
 10   read(nunit,*,end=11) year(i),dumr,nother(i),nkep(i),nmulti(i)
         i=i+1
      goto 10
 11   continue
      close(nunit)
      npt=i-1

      call pgopen('?')
c      call pgask(.false.)
      call pgpage()
      call PGPAP (8.0 ,0.6)
      call pgsch(1.5)
      call pgslw(3)

      call pgvport(0.15,0.85,0.25,0.95)

      nframe=10
      do 15 jj=3*nframe,npt*nframe-nframe+1
      j=jj/nframe
      call pgpage()

      mbin=0.0
      do 16 i=1,j
         mbin=max(mbin,nmulti(i))
 16   continue
      if(j.lt.npt)then
         if(nmulti(j+1).gt.mbin)then
            write(0,*) (nmulti(j+1)-mbin)/nframe*real(jj-j*nframe)
            write(0,*) j*nframe,jj
            mbin=mbin+(nmulti(j+1)-mbin)/nframe*real(jj-j*nframe)
         endif
      endif
      mbin=mbin+mbin*0.15
c      mbin=800.0/(real(nframe)*real(npt-3))*jj

      xrange=(2014.0-1988.0)/(real(nframe)*real(npt-3))*jj+1988.0
      call pgwindow(1988.0,xrange,0.0,mbin)
      call pgbox('BCTS',0.0,0,'BCNTS1',0.0,0)
      call pglabel("","Number of New Planets","")
      call pgptxt((1988.0+xrange)/2.0,-mbin*0.25,0.0,0.5,
     .   "Discovery Year")
c      call pgbin(npt,year,nmulti,.true.)

      do 17 i=1,j
         write(cyear,500) int(year(i))
 500     format(I4)
         call pgptxt(year(i)+0.4/real(npt*nframe)*real(jj),
     .      -mbin*0.05,90.0,1.0,cyear)
 17   continue

c      call pgsch(0.8)
c      if(year(j).ge.2000.0) call pgptxt(2000.0+0.4/real(npt*nframe)*
c     .   real(jj),nmulti(12),90.0,0.0,"HD209458")
c      call pgsch(1.5)

      call pgsci(3)
      do 12 i=1,npt-1
         call pgrect(year(i)-0.5,year(i)+0.5,0.0,nmulti(i))
 12   continue
      call pgsci(2)
      do 13 i=1,npt-1
         call pgrect(year(i)-0.5,year(i)+0.5,0.0,nkep(i))
 13   continue
      call pgsci(4)
      do 14 i=1,npt-1
         call pgrect(year(i)-0.5,year(i)+0.5,0.0,nother(i))
 14   continue
      call pgsci(1)

 15   continue

      call pgclos()

      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 999  end
