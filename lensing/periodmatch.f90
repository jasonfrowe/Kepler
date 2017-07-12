program periodmatch
implicit none
integer :: ntransit,nmatch
real, dimension(3) :: parsin,parsout

!T0 (days)
parsin(1)=66.432951
parsout(1)=66.433828
!Period (days)
parsin(2)=2.23
parsout(2)=4.45999
!duration
parsin(3)=0.135488
parsout(3)=0.0605906

call eventmatch(parsin,parsout,ntransit,nmatch)
write(0,*) "ntransit, nmatch: ",ntransit,nmatch

end program periodmatch

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine eventmatch(parsin,parsout,ntransit,nmatch)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
implicit none
integer :: itransitmatch,ntransit,nmatch,Transitnumber,                 &
   Transitnumber_old
real :: time,Pmin,T0min,Durmin,Pmax,T0max,Durmax,Tstart,Tmission,       &
   Tfinish,Tstep,Phimax,Phimin,Phimaxstart,Phimaxend,Phiminstart,       &
   Phiminend
real, dimension(3) :: parsin,parsout

!figure out which event has the longest period
if(parsin(2)<parsout(2))then
   Pmin=parsin(2)  !save as min and max values
   T0min=parsin(1)
   Durmin=parsin(3)
   Pmax=parsout(2)
   T0max=parsout(1)
   Durmax=parsout(3)
else
   Pmax=parsin(2)
   T0max=parsin(1)
   Durmax=parsin(3)
   Pmin=parsout(2)
   T0min=parsout(1)
   Durmin=parsout(3)
endif

Tstart=T0max-Pmax/2.0 !Time at beginning of observations
Tmission=4*365.25 !length of observations (4 years)
Tfinish=Tstart+Tmission !Time at end of observations
Tstep=30.0/60.0/24.0 !30 minute observations (time-steps)

!use duration to get phases for transit matches
Phimaxstart=(-Durmax)/Pmax+0.5
Phimaxend  =( Durmax)/Pmax+0.5
!write(0,*) phimaxstart,phimaxend
Phiminstart=(-Durmin)/Pmin+0.5
Phiminend  =( Durmin)/Pmin+0.5
!write(0,*) phiminstart,phiminend

!initialize counters
ntransit=0 !number of transits
nmatch=0   !number of matches

Transitnumber_old=1 !keep track of transit number so we know when to reset counters
time=Tstart !time step
itransitmatch=0 !0=no match, 1=we found a match for the event being tested
do while(time<=Tfinish)  !loop over 4 years of observations
   Transitnumber=int((time-T0max+Pmax/2.0)/Pmax)+1 !current transit number
   Phimax=(time-T0max+Pmax/2.0)/Pmax-int((time-T0max+Pmax/2.0)/Pmax) !Phase of event long
   Phimin=(time-T0min+Pmin/2.0)/Pmin-int((time-T0min+Pmin/2.0)/Pmin) !Phase of event short

!  see if we have a transit for max event
   if((Phimax.gt.Phimaxstart).and.(Phimax.lt.Phimaxend))then
   !if we do, check the min event
      if((Phimin.gt.Phiminstart).and.(Phimin.lt.Phiminend))then
         !if we are also in transit for the min event, then we have a match
         itransitmatch=1
      endif
   endif

!when we move to next transit event, we check to see if we had a match previously
   if(Transitnumber.ne.Transitnumber_old)then
      ntransit=ntransit+1  !count number of transits
      if(itransitmatch==1)then
         nmatch=nmatch+1 !if we matched an event, then record it.
      endif
      itransitmatch=0 !reset flag to indicate if we found a match
      Transitnumber_old=Transitnumber !update tracking of current transit number
   endif

!   write(0,*) time,Transitnumber,Phimax,ntransit,nmatch
!   read(5,*)
   time=time+Tstep
enddo

!we should also check if the last event had a chance to produce a match.
if(Phimax.gt.Phimaxend)then !we made it through the last transit (partials are missed)
   ntransit=ntransit+1  !count the transit
   if(itransitmatch==1)then
         nmatch=nmatch+1 !if we matched the last event, count it.
   endif
endif

return
end
