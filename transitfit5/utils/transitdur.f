CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function transitdur(nfit,sol,np)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,col,np
      double precision sol(nfit),Pi,Msun,Rsun,G,aConst,b,Psec,adrs,
     .  cincl,temp(4),tpi,rdr,bb,esw,ecw,eccn,drs_ratio,sinincl

      Pi=acos(-1.d0)   !Pi
      tpi=2.0d0*pi
      Msun=1.9891d30 !kg  mass of Sun
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)

      col=8+(np-1)*10
      bb=sol(col+3)
      b=sqrt(bb)
      Psec=sol(col+2)*8.64d4 !sec ; period of planet
      adrs=1000.0*sol(1)*G*(Psec)**2/(3.0d0*Pi)
      adrs=adrs**(1.0d0/3.0d0) !a/R*
      rdr=sol(col+4)
      esw=sol(col+5)
      ecw=sol(col+6)

      eccn=sqrt(esw*esw+ecw*ecw)
      if(eccn.ge.1.0d0) eccn=0.999d0

      drs_ratio=(1.0d0-eccn*eccn)/(1.0d0+esw) ! r/a at transit
      cincl=b/(adrs*drs_ratio) ! cos(i)

      if(abs(cincl).ge.1.0d0) then
          sinincl=0.0d0
      else
          sinincl=sqrt(1.0d0-cincl*cincl)
      endif

      if((sinincl.gt.0.0d0).and.(adrs*drs_ratio.gt.0.0d0)) then
          temp(1)=Psec * ((1.0d0-eccn*eccn)**1.5d0)/
     .            (Pi*(1.0d0+esw)**2.0d0)
          temp(2)=1.0d0/(adrs*drs_ratio)
          temp(3)=(1.0d0+rdr)**2.0d0 - bb
          if(temp(3).ge.0.0d0) then
              transitdur=temp(1)*asin(temp(2)*sqrt(temp(3))/sinincl)
          else
              transitdur=0.0d0
          endif
      else
          transitdur=0.0d0
      endif

      return
      end
