      subroutine marktransit(np,npt,phase,time,tflag,nfit,sol,ntt,tobs,
     .  omc)
      implicit none
      integer npt,nplanetmax,nmax,nfit
      parameter(nmax=2000000,nplanetmax=10)
      integer tflag(npt),i,np,ntt(nplanetmax),col
      double precision time(npt),sol(nfit),ph1,ph2,transitdur,
     .  toff,tdur,phase(npt),period,tobs(nplanetmax,nmax),
     .  omc(nplanetmax,nmax),tcor(nmax),ttcor,epo,tdurfac
      
      tdur=transitdur(np,nfit,sol)/86400.0d0+0.03
c      write(0,*) 'tdur: ',tdur

      col=10*(np-1)
      epo=sol(9+col)
      period=sol(10+col)

C     tdurfac - 0.5: remove exactly 1 transit duration
C     tdurfac - 1.0: remove +/- 1 transit duration centred on transit
      tdurfac=1.0d0
      ph1=0.75-tdurfac*tdur/period
      if(ph1.lt.0.5)ph1=0.5
      ph2=0.75+tdurfac*tdur/period
      if(ph2.gt.1.0)ph2=1.0
c      write(0,*) "ph1,ph2",ph1,ph2
      
      toff=0.75-(epo/period-int(epo/period))
      if(toff.lt.0.0)toff=toff+1.0
c      write(0,*) "Toff:",toff

      do 24 i=1,npt
         call lininterp(tobs,omc,nplanetmax,nmax,np,ntt,time(i),ttcor)
         tcor(i)=time(i)-ttcor
c         write(0,*) tcor(i),ttcor
 24   continue

      call phasept(npt,tcor,phase,period,toff)
      
      do 10 i=1,npt
c        if(phase(i).lt.0.0d0) phase(i)=phase(i)+1.0d0
c        write(0,*) phase(i),ph1,ph2
        if((phase(i).ge.ph1).and.(phase(i).le.ph2))then           
            tflag(i)=1
c            write(0,*) time(i)
c            read(5,*)
        endif
 10   continue
      
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function transitdur(np,nfit,sol)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer np,nfit
      double precision sol(nfit),b,Psec,G,aConst,Pi,adrs,cincl,temp(4),
     .  bb,rdr
      
      Pi=acos(-1.d0)   !Pi
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)
      
      b=sol(11+10*(np-1))
      bb=b*b
      Psec=sol(10+10*(np-1))*8.64d4 !sec ; period of planet
      adrs=1000.0*sol(1)*G*(Psec)**2/(3.0d0*Pi)
      adrs=adrs**(1.0d0/3.0d0) !a/R*
      cincl=b/adrs !cos(i)
      rdr=sol(12+10*(np-1))
c      write(0,*) bb,adrs,rdr
        
      temp(1)=Psec/Pi
      temp(2)=1.0d0/adrs
      temp(3)=(1+rdr)**2.0-bb
      temp(4)=1-cincl*cincl
C     Transit duration in days
      transitdur=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))
      
      return
      end
