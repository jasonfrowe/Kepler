program timingfigure
implicit none
real :: dt,a,d,p,pfunc,dfunc,tx,ty
real, allocatable, dimension(:) :: bb,px,py
allocate(bb(4))

bb(1)=2.0
bb(2)=1000.0
bb(3)=2.0
bb(4)=1000.0
bb=log10(bb)

call pgopen('?')
call PGPAP (8.0 ,1.0) !use a square 8" across
call pgslw(3) !thicker lines

call pgwindow(bb(1),bb(2),bb(3),bb(4))
CALL PGBOX('BCLNTS1',0.0,0,'BCLNTS1',0.0,0)
call pglabel("Distance (pc)","Period (days)","")

allocate(px(2),py(2))

!first line
call pgsls(1)
call pgsci(1)
dt=1.0  !s
a =50.0 !AU
d=10.0**bb(1)
px(1)=log10(d)
py(1)=log10(pfunc(dt,a,d))
d=10.0**bb(2)
px(2)=log10(d)
py(2)=log10(pfunc(dt,a,d))
call pgline(2,px,py)
tx=(log10(dfunc(dt,a,10.0**bb(3)))+bb(2))/2.0d0-(bb(2)-bb(1))*0.01
ty=(bb(3)+py(2))/2.0d0+(bb(4)-bb(3))*0.01
call PGPTXT (tx,ty,45.0,0.5,"1 s, 50 AU")

!second line
call pgsls(2)
call pgsci(1)
dt=10.0  !s
a =50.0 !AU
d=10.0**bb(1)
px(1)=log10(d)
py(1)=log10(pfunc(dt,a,d))
d=10.0**bb(2)
px(2)=log10(d)
py(2)=log10(pfunc(dt,a,d))
call pgline(2,px,py)
tx=(log10(dfunc(dt,a,10.0**bb(4)))+bb(1))/2.0d0-(bb(2)-bb(1))*0.01
ty=(bb(4)+py(1))/2.0d0+(bb(4)-bb(3))*0.01
call PGPTXT (tx,ty,45.0,0.5,"10 s, 50 AU")

!third line
call pgsls(1)
call pgsci(2)
dt=1.0  !s
a =1.0 !AU
d=10.0**bb(1)
px(1)=log10(d)
py(1)=log10(pfunc(dt,a,d))
d=10.0**bb(2)
px(2)=log10(d)
py(2)=log10(pfunc(dt,a,d))
call pgline(2,px,py)
tx=(log10(dfunc(dt,a,10.0**bb(4)))+bb(1))/2.0d0+(bb(2)-bb(1))*0.025
ty=(bb(4)+py(1))/2.0d0-(bb(4)-bb(3))*0.025
call PGPTXT (tx,ty,45.0,0.5,"1 s, 1 AU")

!third line
call pgsls(4)
call pgsci(1)
dt=60.0  !s
a =50.0 !AU
d=10.0**bb(1)
px(1)=log10(d)
py(1)=log10(pfunc(dt,a,d))
d=10.0**bb(2)
px(2)=log10(d)
py(2)=log10(pfunc(dt,a,d))
call pgline(2,px,py)
tx=(log10(dfunc(dt,a,10.0**bb(4)))+bb(1))/2.0d0-(bb(2)-bb(1))*0.01
ty=(bb(4)+py(1))/2.0d0+(bb(4)-bb(3))*0.01
call PGPTXT (tx,ty,45.0,0.5,"60 s, 50 AU")

!forth line
call pgsls(2)
call pgsci(2)
dt=10.0  !s
a =1.0 !AU
d=10.0**bb(1)
px(1)=log10(d)
py(1)=log10(pfunc(dt,a,d))
d=10.0**bb(2)
px(2)=log10(d)
py(2)=log10(pfunc(dt,a,d))
call pgline(2,px,py)
tx=(log10(dfunc(dt,a,10.0**bb(4)))+bb(1))/2.0d0-(bb(2)-bb(1))*0.01
ty=(bb(4)+py(1))/2.0d0+(bb(4)-bb(3))*0.01
call PGPTXT (tx,ty,45.0,0.5,"10 s, 1 AU")

call pgsci(1)
call pgclos()

end program timingfigure

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
real function dfunc(dt,a,p)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
implicit none
real :: dt,p,a,C
C=0.067

dfunc=p*a*C/dt

return
end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
real function pfunc(dt,a,d)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
implicit none
real :: dt,d,a,C
C=0.067

pfunc=dt*d/(C*a)

return
end



