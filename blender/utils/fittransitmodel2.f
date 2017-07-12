CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function fittransitmodel2(npt,aT,aM,aE,dtype,
     .  nfit,sol,serr,npars,ipars,p,pvary,Dpvary,dchistop,alpha,covar,
     .  ia,inclmin,err)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,i,j,ipars(nfit),npars,it,itmax,npt,
     .  ia(nfit),dtype(npt),goagain
      double precision sol(nfit),serr(nfit,2),p(nfit),
     .  pvary(nfit),Dpvary(nfit),bchi,ochi,
     .  dchi,dchistop,aT(npt),aM(npt),aE(npt),alamda,chisq,
     .  alpha(nfit,nfit),covar(nfit,nfit),alamhi,inclmin,M1,M2,Psec,
     .  asemi,R1,err(nfit),R2,mintime,maxtime,epoch
      logical loop
      include "utils/physcons.f"
      
 36   goagain=0 !if we need to repeat the fitting procedure
      
      M1=sol(1)*Msun !kg ; mass of star
      M2=sol(2)*Mjup !kg ; mass of planet
      Psec=sol(5)*24.0*60.0*60.0 !sec ; period of planet
      asemi=(Psec*Psec*G*(M1+M2)/(4.0*Pi*Pi))**(1.0/3.0)
      R1=sol(3)*Rsun
      R2=sol(3)*Rjup
      inclmin=180.0*tan((R1+R2)/asemi)/pi
      inclmin=90.0-inclmin
      
      mintime=aT(1)
      maxtime=aT(1)
      do 34 i=2,npt
        mintime=min(mintime,aT(i))
        maxtime=max(maxtime,aT(i))
 34   continue
      
      j=0
      do 10 i=1,nfit  !pick out parameters that we are going to fit
        if((serr(i,2).ne.0.0).and.(i.ne.18))then !exclude DIL from fits
            j=j+1
            ipars(j)=i !store indices of fitted variables
        endif
 10   continue
      npars=j
 
      write(0,500) "To Fit: ",(ipars(i),i=1,npars)
 500  format(A8,16(1X,I2))

      do 22 i=1,npars
        p(i)=sol(ipars(i))
        pvary(i)=Dpvary(ipars(i))
        ia(i)=1 !fit parameters fed to mrqmin
 22   continue

      alamhi=100000.
      alamda=-1.0
      loop=.true.
      itmax=75  !takes about 5 iterations to get minimum.
      it=1
      do 30 while(loop)
C       Going to put parameters in order of what I want fitted.
C       Initial guesses go into element nfit of p
         call mrqmin(aT,aM,aE,dtype,npt,p,ia,npars,covar,alpha,nfit,
     .      chisq,alamda)
         do 32 i=1,npars
            if(ipars(i).eq.3) p(i)=abs(p(i))
            if(ipars(i).eq.4) p(i)=abs(p(i))
            if(ipars(i).eq.5) then
c                epoch=(sol(7)-pi)/(2.0*pi)
c                if(epoch.lt.0.0) epoch=epoch+2.0*pi
c                if(epoch.gt.2.0*pi) epoch=epoch-2.0*pi
c                epoch=epoch*p(i)
c                epoch=epoch+p(i)*(floor((mintime-epoch)/p(i))+1)
cc                write(0,*) "epoch:",epoch
                epoch=sol(7)
C               Conditions for ending Period fit
c                if(epoch+p(i).gt.maxtime)then
c                    ia(i)=0
c                    serr(5,2)=0.0d0
c                    write(6,*) "epoch:",epoch,epoch+p(i)
c                endif
            endif
            if(ipars(i).eq.6) then
c                write(0,*) "INCL:",p(i),inclmin
                if(p(i).gt.90.0) p(i)=180.0-p(i)
c                if (p(i).lt.inclmin) then
c                    p(i)=inclmin!(90.0+inclmin)/2.0
c                    ia(i)=0 !stop fitting inclination
c                    serr(6,2)=0.0d0
c                endif
                if (p(i).gt.89.95) then !stop when we get to 90 deg.
                    p(i)=90.0d0
                    ia(i)=0
                    serr(6,2)=0.0d0 !turn off inclination
                    serr(3,2)=-1.0d0 !turn on radius
                    do 35 j=1,npars
                        goagain=1 !repeat the fit if we need to redo rad
                        if(ipars(j).eq.3)then
                            ia(j)=1 !fit stellar radius
                            goagain=0 !radius is being fit
                        endif
 35                 continue
                endif
c                write(0,*) "INCL:",p(i),inclmin
            endif
C           Old ECN,WWW stuff
c            if(ipars(i).eq.14) then
c                p(i)=abs(p(i))
c            endif
c            if(ipars(i).eq.15) then
c                if(p(i).lt.0.0d0) p(i)=p(i)+360.0d0
c                if(p(i).gt.3.0d2) p(i)=p(i)-360.0d0
c            endif
 32      continue
 
         write(0,504) "Best fit for iteration # ",it
 504     format(A25,I4)
         bchi=chisq
         if(it.gt.1)then
            dchi=(ochi-bchi)!/ochi
         else
            dchi=1.0e30
         endif
         write(0,*) "alamda",alamda,dchi
         write(0,503) (p(i),i=1,npars),bchi,dchi
         ochi=bchi
 503     format(17(1PE10.3,1X))
         it=it+1
c         if(it.gt.itmax) loop=.false.
c         if((abs(dchi).lt.dchistop).and.(abs(dchi).gt.0.0)) loop=.false.
c         if(abs(dchi).lt.dchistop) then
          if((it.gt.itmax).or.
     .    ((dchi.lt.dchistop).and.(dchi.gt.0.0d0)).or.
     .    (alamda.gt.alamhi))then
            loop=.false.
            alamda=0.0
            call mrqmin(aT,aM,aE,dtype,npt,p,ia,npars,covar,alpha,nfit,
     .      chisq,alamda)
         endif
 30   enddo
 
c      bchi=y(tp(1))
      do 33 i=1,nfit
        err(i)=0.0
 33   continue
      do 31 i=1,npars
         sol(ipars(i))=p(i)
         err(ipars(i))=sqrt(covar(i,i)*bchi/dble(npt-1))
c         Dpvary(ipars(i))=pvary(i)
 31   continue
 
      if(goagain.gt.1) goto 36 !repeat the fit
 
      write(0,502) (sol(i),i=1,nfit)
 502  format(13(F7.3,1X))
      write(0,*) "Chi-Sq:",bchi
      fittransitmodel2=bchi
 
      return
      end