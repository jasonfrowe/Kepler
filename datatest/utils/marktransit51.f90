subroutine marktransit(np,nplanet,npt,time,tflag,nfit,sol,ntt,tobs,omc)
use precision
implicit none
!import vars
integer :: np,npt,nfit,nplanet
integer, dimension(:) :: tflag,ntt
real(double), dimension(:) :: time,sol
real(double), dimension(:,:) :: tobs,omc
!local vars
integer :: col,i
real(double) :: tdur,transitdur,ph1,ph2,toff,epo,period,ttcor
real(double), allocatable, dimension(:) :: phase,tcor

tdur=transitdur(np,nfit,sol)/86400.0d0+0.03 !extra 0.03 to include 30-min cadence

col=10*(np-1)
epo=sol(9+col)
period=sol(10+col)
      
ph1=0.75-0.5d0*tdur/period
if(ph1.lt.0.5)ph1=0.5
ph2=0.75+0.5d0*tdur/period
if(ph2.gt.1.0)ph2=1.0
      
toff=0.75-(epo/period-int(epo/period))
if(toff.lt.0.0)toff=toff+1.0

allocate(tcor(npt))
!$OMP PARALLEL DO PRIVATE(ttcor)
do i=1,npt
    !routine from ttcor.
	call lininterp(tobs,omc,nplanet,npt,np,ntt,time(i),ttcor)
    tcor(i)=time(i)-ttcor
enddo
!$OMP END PARALLEL DO

allocate(phase(npt))
call phasept(npt,tcor,phase,period,toff)


!$OMP PARALLEL DO
do i=1,npt
   	if((phase(i).ge.ph1).and.(phase(i).le.ph2))then           
		tflag(i)=1
    endif
enddo
!$OMP END PARALLEL DO

return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function transitdur(np,nfit,sol)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use precision
implicit none
integer :: np,nfit
real(double) :: sol(nfit),b,Psec,G,aConst,Pi,adrs,cincl,temp(4),bb,rdr, &
 transitdur
      
Pi=acos(-1.d0)   !Pi
G=6.674d-11 !N m^2 kg^-2  Gravitation constant
aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)
      
bb=sol(11+10*(np-1))
b=sqrt(bb)
Psec=sol(10+10*(np-1))*8.64d4 !sec ; period of planet
adrs=1000.0*sol(1)*G*(Psec)**2/(3.0d0*Pi)
adrs=adrs**(1.0d0/3.0d0) !a/R*
cincl=b/adrs !cos(i)
rdr=sol(12+10*(np-1))

temp(1)=Psec/Pi
temp(2)=1.0d0/adrs
temp(3)=(1+rdr)**2.0-bb
temp(4)=1-cincl*cincl
! Transit duration in days
transitdur=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))
      
return
end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine phasept(npt,time,phase,period,toff)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     npt - number of data points
!     time - times of observations
!     mag - the observations
!     phase - returned phase of observation
!     period - fixed period for data
use precision
implicit none
integer :: npt
real(double) :: time(npt),phase(npt),period,toff

integer i
real(double) temp

!$OMP PARALLEL DO PRIVATE(temp)
do i=1,npt
	temp=time(i)
!   Get the phase
    phase(i)=temp/period-int(temp/period)
!   apply optional phase offset to make plot pretty
    phase(i)=phase(i)+toff
!   make sure phase is between 0 and 1
	if(phase(i).lt.0.0) phase(i)=phase(i)+1.0
	if(phase(i).gt.1.0) phase(i)=phase(i)-1.0
enddo
!$OMP END PARALLEL DO

return
end
