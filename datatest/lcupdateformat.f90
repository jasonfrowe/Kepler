program lcupdateformat
use precision
implicit none
integer :: iargc,nunit,filestatus,itype,i
real(double) :: time,flux,ferr
character(180) :: filename

itype=0 !mark all data as Long Cadence

if(iargc().lt.1)then
    write(0,*) "Usage: lcupdateformat filename > filename_new"
    stop
endif

call getarg(1,filename)

nunit=10 !unit number for file list
open(unit=nunit,file=filename,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",filename
   stop
endif

i=0 !line counter
do
    read(nunit,*,iostat=filestatus) time,flux,ferr
	if(filestatus == 0) then
		i=i+1
		write(6,500) time,flux,ferr,itype
	elseif(filestatus == -1) then
        exit  !successively break from data read loop.
	else
        write(0,*) "File Error!! Line:",i+1
        write(0,900) "iostat: ",filestatus
		900 format(A8,I3)
        stop
	endif
enddo

close(nunit)

500 format(3(F17.11,1X),I1)

end

