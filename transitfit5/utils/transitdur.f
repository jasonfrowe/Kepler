CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function transitdur(nfit,sol,np)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Total transit duration T14 (seconds) of planet np, for an eccentric
C     orbit (Kipping 2010, MNRAS 407, 301, eq. 15).  sol(col+3) is b (not b^2) and
C     sol(col+5), sol(col+6) are sqrt(e)cos(omega), sqrt(e)sin(omega),
C     with omega the argument of periastron of the RV literature
C     (transit at true anomaly f = pi/2 - omega), as in transitmodel.f.
      implicit none
      integer nfit,col,np
      double precision sol(nfit),Pi,G,b,Psec,adrs,cincl,temp(3),rdr,bb,
     .  esw,ecw,eccn,esinw,drs_ratio,sinincl

      Pi=acos(-1.d0)   !Pi
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant

      col=8+(np-1)*10
      b=abs(sol(col+3))
      bb=b*b
      Psec=sol(col+2)*8.64d4 !sec ; period of planet
      adrs=1000.0d0*sol(1)*G*(Psec)**2/(3.0d0*Pi)
      adrs=adrs**(1.0d0/3.0d0) !a/R*
      rdr=abs(sol(col+4))
      ecw=sol(col+5)
      esw=sol(col+6)

      eccn=ecw*ecw+esw*esw
      if(eccn.ge.1.0d0) eccn=0.99d0
      esinw=sqrt(eccn)*esw    ! e sin(omega)

      drs_ratio=(1.0d0-eccn*eccn)/(1.0d0+esinw) ! r/a at transit
      cincl=b/(adrs*drs_ratio) ! cos(i)

      if(abs(cincl).ge.1.0d0) then
          sinincl=0.0d0
      else
          sinincl=sqrt(1.0d0-cincl*cincl)
      endif

      transitdur=0.0d0
      if((sinincl.gt.0.0d0).and.(adrs*drs_ratio.gt.0.0d0)) then
          temp(1)=Psec/Pi*drs_ratio**2/sqrt(1.0d0-eccn*eccn)
          temp(2)=1.0d0/(adrs*drs_ratio*sinincl)
          temp(3)=(1.0d0+rdr)**2-bb
          if(temp(3).ge.0.0d0) then
              transitdur=temp(1)*asin(min(1.0d0,temp(2)*sqrt(temp(3))))
          endif
      endif

      return
      end
