CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine oddeven(npt,aT,aM,aE,dtype,nfit,sol,serr,npars,
     .  ipars,p,pvary,Dpvary,dchistop,alpha,covar,ia,inclmin,err,phase,
     .  time2,mag2,merr2,sol2,serr2,doe,err2,dtype2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nfit,i,npt2,n1,n2
      integer ia(nfit),dtype(npt),dtype2(npt)
      real*8 aT(npt),aM(npt),aE(npt),sol(nfit),serr(nfit,2),npars,
     .  ipars(nfit),p(nfit),pvary(nfit),Dpvary(nfit),alpha(nfit,nfit),
     .  covar(nfit,nfit),dchistop,inclmin,err(nfit),epoch,phi,
     .  phase(npt),per,epoff,time2(npt),mag2(npt),merr2(npt),
     .  opra,epra,oerr,eerr,sol2(nfit),serr2(nfit,2),doe,tdur,
     .  transitdur,phi2,temp(2),ph,err2(nfit)
      include "physcons.f"
     
      doe=0.0d0
      
      tdur=transitdur(nfit,sol)/86400.0d0
     
      epoch=sol(7)
      per=sol(5)
      phi=pi-2.0*pi*(epoch/per-floor(epoch/per))
      temp(1)=pi-2.0*pi*((epoch+tdur)/per-floor((epoch+tdur)/per))
      temp(2)=pi-2.0*pi*((epoch+tdur)/per-floor((epoch+tdur)/per))
      phi2=min(abs(phi-temp(1)),abs(phi-temp(2)))
      epoff=(phi-0.5)*per
      
      n1=0
      n2=0
      npt2=0
      do 10 i=1,npt
        phase(i)=floor((aT(i)-epoff)/per)
        ph=(aT(i)-epoff)/per-phase(i)
        if(mod(phase(i)+1.0d0,2.0d0).eq.0)then
            if(abs(0.5d0-ph).lt.phi2)n1=n1+1
            npt2=npt2+1
            time2(npt2)=aT(i)
            mag2(npt2)=aM(i)
            merr2(npt2)=aT(i)
            dtype2(npt2)=dtype(i)
        else
            if(abs(0.5d0-ph).lt.phi2)n2=n2+1
        endif
 10   continue
c      write(6,*) phi2
c      write(6,*) "n:",n1,n2
      if((n1.eq.0).or.(n2.eq.0))return
      if(npt2.le.2) return !lets get at least 2 data points..
 
      do 11 i=1,nfit
        sol2(i)=sol(i)
        serr2(i,1)=0.0!serr(i,1)
        serr2(i,2)=0.0!serr(i,2)
 11   continue
      serr2(4,2)=-1.0
 
      dchistop=1.0d-10 !criteria for good fit.
 
      call fittransitmodel2(npt2,time2,mag2,merr2,dtype2,nfit,sol2,
     .  serr2,npars,ipars,p,pvary,Dpvary,dchistop,alpha,covar,ia,
     .  inclmin,err2)
      oerr=err2(4)
      opra=sol2(4)
      
      npt2=0
      do 12 i=1,npt
        phase(i)=floor((aT(i)-epoff)/per)
        if(mod(phase(i),2.0d0).eq.0)then
            npt2=npt2+1
            time2(npt2)=aT(i)
            mag2(npt2)=aM(i)
            merr2(npt2)=aE(i)
            dtype2(npt2)=dtype(i)
        endif
 12   continue      
      
      if(npt2.le.2) return !lets get at least 2 data points..
      
      do 13 i=1,nfit
        sol2(i)=sol(i)
        serr2(i,1)=0.0!serr(i,1)
        serr2(i,2)=0.0!serr(i,2)
 13   continue
      serr2(4,2)=-1.0
 
      dchistop=1.0d-10 !criteria for good fit.
 
      call fittransitmodel2(npt2,time2,mag2,merr2,dtype2,nfit,sol2,
     .  serr2,npars,ipars,p,pvary,Dpvary,dchistop,alpha,covar,ia,
     .  inclmin,err2)
      eerr=err2(4)
      epra=sol2(4)      
      
      write(0,*) oerr,opra
      write(0,*) eerr,epra
      doe=abs(opra-epra)/sqrt(oerr**2.0+eerr**2.0)
      write(0,*) "doe:", doe
      
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function transitdur(nfit,sol)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This is for circular orbits only.
      implicit none
      integer nfit
      double precision sol(nfit),Psec,M1,M2,R1,R2,asemi,temp(4),incl
C     Read in physical constants (Pi,G,Msun,..etc.)      
      include "physcons.f"
      
      M1=sol(1)*Msun !kg ; mass of star
      M2=sol(2)*Mjup !kg ; mass of planet
      Psec=sol(5)*24.0*60.0*60.0 !sec ; period of planet
      asemi=(Psec*Psec*G*(M1+M2)/(4.0*Pi*Pi))**(1.0/3.0) !m
      R1=sol(3)*Rsun  !radius of star
      R2=sol(4)*Rjup  !radius of planet
      incl=Pi*sol(6)/180.0d0
      
      temp(1)=Psec/Pi
      temp(2)=R1/asemi
      temp(3)=(1+(R2/R1))**2.0-((asemi/R1)*cos(incl))**2.0
      temp(4)=1-cos(incl)*cos(incl)
      
      transitdur=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))
      
      return
      end
      
