      subroutine marktransit(npt,phase,time,flux,tflag,nfit,sol)
      implicit none
      integer npt,nfit,tflag(npt),i
      double precision time(npt),flux(npt),sol(nfit),ph1,ph2,transitdur,
     .  toff,tdur,phase(npt),period
      
      tdur=transitdur(nfit,sol)/86400.0d0+0.03
      
      ph1=0.75-0.5d0*tdur/sol(5)
      if(ph1.lt.0.5)ph1=0.5
      ph2=0.75+0.5d0*tdur/sol(5)
      if(ph2.gt.1.0)ph2=1.0
      
      toff=0.75-(sol(7)/sol(5)-int(sol(7)/sol(5)))
      if(toff.lt.0.0)toff=toff+1.0
c      write(0,*) "Toff:",toff
      
      period=sol(5)
      call phasept(npt,time,phase,period,toff)
      
      do 10 i=1,npt
        if((phase(i).ge.ph1).and.(phase(i).le.ph2))then
            tflag(i)=1
c            write(0,*) time(i)
        endif
 10   continue
      
      
      return
      end