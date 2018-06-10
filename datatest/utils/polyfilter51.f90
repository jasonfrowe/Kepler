subroutine polyfilter(npt,time,flux,ferr,tflag,boxbin,nfitp)
use precision
implicit none
!import vars
integer :: npt,nfitp
integer, dimension(:) :: tflag
real(double) :: boxbin
real(double), dimension(:) :: time,flux,ferr
!local vars
integer :: i,nt,i1,i2,npt2,nmax,j,k
integer, allocatable, dimension(:) :: ts,ts2,tfsort
real(double) :: bbd2,meddt,off,tzero
real(double), allocatable, dimension(:) :: dt,x,y,z,offset,tsort,fsort,fesort

!interface
!	subroutine polydetrend(npt,x,y,z,nfitp,tzero,off)
!		use precision
!		implicit none
!		integer :: npt,nfitp
!		real(double) :: off,tzero
!		real(double), dimension(:) :: x,y,z
!	end subroutine polydetrend
!end interface

bbd2=boxbin/2.0d0

!make sure data is sorted by time
allocate(ts(npt))
call rqsort(npt,time,ts)

allocate(tsort(npt),fsort(npt),fesort(npt),tfsort(npt))
!$OMP PARALLEL DO
do i=1,npt
	tsort(i)=time(ts(i))
	fsort(i)=flux(ts(i))
	fesort(i)=ferr(ts(i))
	tfsort(i)=tflag(ts(i))
enddo
!$OMP END PARALLEL DO
deallocate(ts)


allocate(dt(npt-1))
!$OMP PARALLEL DO
do i=1,npt-1
	dt(i)=tsort(i+1)-tsort(i)
enddo
!$OMP END PARALLEL DO

allocate(ts2(npt-1))
call rqsort(npt-1,dt,ts2)
meddt=dt(ts2((npt-1)/2))
deallocate(ts2,dt)

nt=int(bbd2/meddt)

if(nt.lt.nfitp+1)then
	write(0,*) "Warning:  boxbin is likely too small or nfitp too large"
endif

nmax=2*nt+1
allocate(x(nmax),y(nmax),z(nmax),offset(npt))

!$OMP PARALLEL DO FIRSTPRIVATE(tzero,i1,i2,npt2,j,k,x,y,z,off)
do i=1,npt
	!write(0,*) i
	tzero=tsort(i)
	i1=max(1,i-nt) !start of array
	i2=min(npt,i+nt) !end of array
	npt2=i2-i1+1
	k=0
	do j=i1,i2
		if(tfsort(j).eq.0)then !only use data outside of transit.
			k=k+1
			x(k)=tsort(j)
			y(k)=fsort(j)
			z(k)=fesort(j)
		endif
	enddo
	if(npt2.gt.nfitp+1)then
!       write(0,*) "into polydetrend"
        call polydetrend(k,x,y,z,nfitp,tzero,off)
!       write(0,*) "out of polydetrend"
        offset(i)=off
    elseif(npt2.eq.0)then
        offset(i)=0.0d0 !if no data, then no correction.. duh..
    else
        offset(i)=Sum(y(1:npt2))/dble(npt2) !fall back to simple mean if too few points
    endif
enddo
!$OMP END PARALLEL DO

write(0,*) "npt: ",npt
!$OMP PARALLEL DO
do i=1,npt
	time(i)=tsort(i)
	flux(i)=fsort(i)-offset(i)
	ferr(i)=ferr(i)
enddo
!$OMP END PARALLEL DO

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine polydetrend(npt,time,mag,merr,nfit,tzero,off)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
!import vars
integer :: npt,nfit
real(double) :: off,tzero
real(double), dimension(npt) :: time,mag,merr
!local vars
integer :: i,j,ii,maxiter
integer, allocatable, dimension(:) :: ia
real(double) :: meanT,chisq,T,std,stdev,mx,dchi,ochi,sigcut
real(double), allocatable, dimension(:) :: ans,work,merr2
real(double), allocatable, dimension(:,:) :: covar

sigcut=3.0

allocate(ia(nfit),covar(nfit,nfit),ans(nfit),work(npt),merr2(npt))
ia=1 !fit all variables

meanT=Sum(time(1:npt))/dble(npt) !mean Time,
time=time-meanT !remove mean.

call lfit(time,mag,merr,npt,ans,ia,nfit,covar,nfit,chisq)

ii=0
dchi=1.0 !initiate delta chisq
ochi=chisq
maxiter=10
do while((ii.lt.maxiter).and.(dchi.gt.0.1))

	do i=1,npt
    	T=0.
        do j = 1, nfit
            T = ans(j)*((time(i))**dble(j-1)) +  T
 		enddo
        work(i)=mag(i)-T
 	enddo

	std=stdev(npt,work,mx) !standard deviation

	j=0
    do i=1,npt
    	if(abs(work(i)-mx).lt.sigcut*std)then
            j=j+1
            time(j)=time(i)
            mag(j)=mag(i)
            merr(j)=merr(i)
            merr2(j)=merr(i)+abs(work(i))
         endif
 	enddo
 	!if(npt.eq.j) exit !break loop, because no data was clipped.
    npt=j
    call lfit(time,mag,merr2,npt,ans,ia,nfit,covar,nfit,chisq)
    dchi=abs(chisq-ochi)
    ochi=chisq
    ii=ii+1  !count number of iterations to avoid infinite loops

enddo

T=0.
do j=1,nfit
	T = ans(j)*((tzero-meanT)**dble(j-1)) +  T
enddo
off=T

return
end

