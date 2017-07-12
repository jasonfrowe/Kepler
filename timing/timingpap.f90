program timingpap
use precision
implicit none
integer npt,iargc,filestatus,nfile,nunit,ntt,nttmax,i,j
integer, allocatable, dimension(:) :: p
real(double) :: chisq,rchisq,tprec,t,f,ferr,esum,tsum,fsum,rnpt
real(double), allocatable, dimension(:) :: time,tt,tterr,rtime,rtt,rtterr,temp
character(80) :: listfile,ttfile

interface
   subroutine readttfile(ttfile,ntt,rtime,rtt,rtterr)
      use precision
      implicit none
      integer :: ntt
      real(double), dimension(:) :: rtime,rtt,rtterr
      character(80) :: ttfile
   end subroutine readttfile
end interface
interface
   subroutine chisqtest(ntt,rtt,rtterr,chisq)
      use precision
      implicit none
      integer :: ntt
      real(double) :: chisq
      real(double), dimension(:) :: rtt,rtterr
   end subroutine chisqtest
end interface

if(iargc().lt.1)then
   write(0,*) "Usage: timingpap ttfiles.list"
   stop
endif

call getarg(1,listfile)  !get filename for photometry
nunit=10 !unit number for file list
open(unit=nunit,file=listfile,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",listfile
   stop
endif

npt=2200000  !number of expecting timing measurements
allocate(time(npt),tt(npt),tterr(npt))

nttmax=6000 !maximum number of timings in any one file
allocate(rtime(nttmax),rtt(nttmax),rtterr(nttmax))

nfile=0
npt=0
do
   read(nunit,500,iostat=filestatus) ttfile
   500 format(A80)
   if(filestatus == 0) then
      nfile=nfile+1 !count number of files

      call readttfile(ttfile,ntt,rtime,rtt,rtterr) !read in timing data
!     need to check the S/N of the detection.. ignore FAs
      call chisqtest(ntt,rtt,rtterr,chisq)
      rchisq=chisq/dble(ntt)
!     test the timing file
      if((rchisq.gt.0.15).and.(rchisq.lt.3.5))then
         do i=1,ntt
            npt=npt+1
            time(npt)=rtime(i)
            tt(npt)=rtt(i)
            tterr(npt)=rtterr(i)
         enddo
      endif
!      write(6,*) chisq
!      write(6,501) ttfile,ntt,chisq/dble(ntt)
      501 format(A30,I5,1X,F10.3)
!      read(5,*)

      cycle
   elseif(filestatus == -1) then
      exit  !successively break from data read loop.
   else
      write(0,*) "File Error!!",listfile
      write(0,900) "iostat: ",filestatus
      900 format(A8,I3)
      stop
   endif
enddo
close(nunit)
write(0,*) "nfiles: ",nfile
write(0,*) "Timings: ",npt

allocate(p(npt))
call rqsort(npt,time,p)
allocate(temp(npt))
do i=1,npt
   temp(i)=time(p(i))
enddo
time(1:npt)=temp(1:npt)
do i=1,npt
   temp(i)=tt(p(i))
enddo
tt(1:npt)=temp(1:npt)
do i=1,npt
   temp(i)=tterr(p(i))
enddo
tterr(1:npt)=temp(1:npt)


tprec=10.0/(24.0*60.0*60.0) !5 second timing precision requirement

i=1
j=2
do while(j.lt.npt)
   ferr=99.9e30
   do while(ferr.gt.tprec)
      !calculate Weighted Sums
      tsum=Sum(time(i:j)/tterr(i:j)**2.0d0)
      fsum=Sum(tt(i:j)/tterr(i:j)**2.0d0)
!      rnpt=dble(j-i+1)
      esum=Sum(1.0d0/tterr(i:j)**2.0d0)
      !calcalate Weighted Errors
      t=tsum/esum
      f=fsum/esum
!      ferr=sqrt(rnpt)/esum
      ferr=1.0d0/sqrt(esum)
!      write(0,*) t,f,ferr
!      read(5,*)
      j=j+1
      if(j.gt.npt) ferr=0.0d0
   enddo
   write(0,'(I7,1X,I7,3(1X,F12.7))') i,j,t,f,ferr
!   read(5,*)
   if(ferr.gt.0.0d0) write(6,'(3(1X,F13.8))') t,f,ferr
   i=j+1
   j=j+2
enddo


!do i=1,npt
!   write(6,*) time(i),tt(i),tterr(i)
!enddo

end program timingpap
