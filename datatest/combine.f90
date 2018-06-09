program combine
implicit none
integer nmax,nptSC,nptLC,nunit,filestatus,iargc,i,flag,j,npt
integer, allocatable, dimension(:) :: itype
real(8) :: dt,dtthres
real(8), allocatable, dimension(:) :: timeSC,fluxSC,fluxerrSC,timeLC,fluxLC,fluxerrLC, &
    time,flux,fluxerr
character(80) :: SCname,LCname,cline

if(iargc().lt.2)then
    write(0,*) "Usage: combine SCname LCname > Combined"
    stop
endif

nmax=2000000

call getarg(1,SCname)
call getarg(2,LCname)

nunit=10 !unit number for file list
open(unit=nunit,file=SCname,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",SCname
   stop
endif

allocate(timeSC(nmax),fluxSC(nmax),fluxerrSC(nmax))
i=1
!read in data from file line by line
do
    read(nunit,*,iostat=filestatus) timeSC(i),fluxSC(i),fluxerrSC(i)
	if(filestatus == 0) then
		if(i.gt.nmax)then
            write(0,*) "Error : SCname has too many entries"
            write(0,*) "i     : ",i
            write(0,*) "nmax  : ",nmax
            stop
    	endif
        !timeSC(i)=timeSC(i)+0.5d0-54900.0 !this is for sharing with external
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
nptSC=i-1
write(0,*) "nptSC: ",nptSC
close(nunit)

open(unit=nunit,file=LCname,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",LCname
   stop
endif

allocate(timeLC(nmax),fluxLC(nmax),fluxerrLC(nmax))
i=1
!read in data from file line by line
do
    read(nunit,*,iostat=filestatus) timeLC(i),fluxLC(i),fluxerrLC(i)
    if(filestatus == 0) then
        if(i.gt.nmax)then
            write(0,*) "Error : LCname has too many entries"
            write(0,*) "i     : ",i
            write(0,*) "nmax  : ",nmax
            stop
        endif
        !timeLC(i)=timeLC(i)+0.5d0-54900.0 !this is for sharing with external
        i=i+1
    elseif(filestatus == -1) then
        exit  !successively break from data read loop.
    else
        write(0,*) "File Error!! Line:",i+1
        write(0,900) "iostat: ",filestatus
        stop
    endif
enddo
nptLC=i-1
write(0,*) "nptLC: ",nptLC
close(nunit)


allocate(time(nmax),flux(nmax),fluxerr(nmax),itype(nmax))
npt=0
dtthres=2.1d-2
do i=1,nptLC
    flag=0
	do j=1,nptSC
		if(abs(timeLC(i)-timeSC(j)).lt.dtthres)then
            flag=1
		endif
	end do
	if(flag.eq.0)then
        npt=npt+1
        time(npt)=timeLC(i)
        flux(npt)=fluxLC(i)
        fluxerr(npt)=fluxerrLC(i)
        itype(npt)=0
	endif
enddo
write(0,*) "npt: ",npt

time(npt+1:npt+nptSC)=timeSC(1:nptSC)
flux(npt+1:npt+nptSC)=fluxSC(1:nptSC)
fluxerr(npt+1:npt+nptSC)=fluxerrSC(1:nptSC)
itype(npt+1:npt+nptSC)=1

npt=npt+nptSC

do i=1,npt
    write(6,500) time(i),flux(i),fluxerr(i),itype(i)
end do

500 format(F17.11,1X,F17.11,1X,F17.11,1X,I1,1X,F17.11)


end