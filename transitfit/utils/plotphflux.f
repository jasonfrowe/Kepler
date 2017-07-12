CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotphflux(npt,time,mag,merr,period,x1,x2,pan,sym,
     .  nmax,px,py,pz,phase,toff,nerrplot,bbp,tdepth,mag2,merr2,bins,
     .  nfit,sol,err,nsensor,Thpars,flag,rad,tmodel,exptime,dtype,tdur)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Plots the phased light curve
C     npt - number of data points
C     time - times of observations
C     mag  - the observations
C     merr - observational error
C     period - fixed period of the data
C     x1,x2 - phase bounds to plot
C     panx,pany - pgplot panel to plot light curve
C     sym - plotting symbol to use
C     toff - phase offset to center plots
C     nerrplot - 0: plot errors, 1: do not plot errors
C     nmax,px,py,pz,phase - workspace variables
C     bbp : bounds used for plot
      implicit none
      integer npt,pan,i,sym,nmax,nerrplot,npoints,npt2,bins,nfit,
     .  nsensor,flag(npt),j,dtype(npt),nres
      double precision time(npt),mag(npt),period,phase(nmax),mmax,mmin,
     .  ax1,ax2,stdev,std,sigcut,mean,merr(npt),x1,x2,bbp(4),width,
     .  tdepth,mag2(npt),merr2(npt),sol(nfit),err(nfit),Thpars(3,2),
     .  sn,transn,redj,Pi,G,aConst,Msun,Mjup,M1,M2,Psec,asemi,incl,
     .  Mearth,R1,Rsun,b,rad,exptime(npt),tmodel(npt),tdur,ph1,ph2
      real px(nmax),py(nmax),pz(nmax),toff,rbb(4)
      character*80 label
      
      nres=0 !set to 1 to plot residuals.
      
      Pi=acos(-1.d0)   !Pi
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)
      Msun=1.9891d30 !kg  mass of Sun
      Mearth=5.974d24 !kg mass of Earth
      Mjup=317.833d0*Mearth !kg  mass of Jupiter  
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      
      redj=0.09205 !radius of Earth divided by Jupiter
      sigcut=3.0
      mmin= 99.9d30
      mmax=-99.9d30

      ph1=0.75-1.1d0*tdur/sol(5)
      if(ph1.lt.0.5)ph1=0.5
      ph2=0.75+1.1d0*tdur/sol(5)
      if(ph2.gt.1.0)ph2=1.0

      call phasept(npt,time,phase,period,toff)
      mean=0.0d0
      j=0
      do 5 i=1,npt
         if((phase(i).lt.ph1).or.(phase(i).gt.ph2))then
            j=j+1
            mean=mean+mag(i)
            mag2(j)=mag(i)
            merr2(j)=merr(i)
         endif
 5    continue
      mean=mean/dble(j)

      std=stdev(j,mag2,mean)
c      call pgpanl(panx,pany)
      call panel(pan)
      npoints=0
      do 10 i=1,npt
        if((phase(i).ge.x1).and.(phase(i).le.x2))then
            if((mag(i).gt.mmax).and.
     .          (abs(mag(i)-mean).lt.sigcut*std))mmax=mag(i)
            if((mag(i).lt.mmin).and.
     .          (abs(mag(i)-mean).lt.sigcut*std))mmin=mag(i)
            if(abs(mag(i)-mean).lt.sigcut*std)npoints=npoints+1
        endif
 10   continue
c      write(0,*) "npoints:",npoints

      do 6 i=1,npt
         mag2(i)=mag(i)
         merr2(i)=merr(i)
 6    continue
      
c      write(0,*) "mmin",mmin
      if(1.0d0-mmin.lt.tdepth) mmin=1.0d0-tdepth-std
c      write(0,*) "mmin",mmin

      if((nres.eq.1).and.(pan.eq.3))then
        ax2=(mmax+0.20*(mmax-mmin)+3.0*std)
      else
        ax2=(mmax+0.10*(mmax-mmin))
      endif
      ax1=(mmin-0.10*(mmax-mmin))
c      ax1=0.9994
c      ax2=1.0002
      bbp(1)=x1
      bbp(2)=x2
      bbp(3)=ax1
      bbp(4)=ax2
      call pgwindow(real(x1),real(x2),real(ax1),real(ax2))
c      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
      if(pan.le.2) then
        call pgptxt(real((x1+x2)/2.0),real(ax1-0.25*(ax2-ax1)),
     .      0.0,0.5,"Phase")
        call pgptxt(real(x1-0.1*(x2-x1)),real((ax2+ax1)/2),
     .      90.0,0.5,"Flux")
      elseif(pan.eq.3) then
        call pgptxt(real((x1+x2)/2.0),real(ax1-0.12*(ax2-ax1)),
     .      0.0,0.5,"Time from mid-transit (hours)")
c1        call pgptxt(real(x1-0.2*(x2-x1)),real((ax2+ax1)/2),
c1     .      90.0,0.5,"Flux")
      endif
      do 11 i=1,npt
        px(i)=real(phase(i))
        py(i)=real(mag(i))
        pz(i)=real(merr(i))
 11   continue
c      if(npoints.lt.1000) then 
c        sym=17
c      else
        sym=1
c      endif
      call pgpt(npt,px,py,sym)
      if(nerrplot.eq.0) call pgerrb(6,npt,px,py,pz,1.0)
      
C     Plot residuals       
      if((pan.eq.3).and.(nres.eq.1)) then
        call transitmodel(npt,time,exptime,dtype,tmodel,nfit,sol)
        do 14 i=1,npt
            py(i)=real(tmodel(i)+2.5*log10(mag(i))+ax2/2.0+0.5)
c            write(0,*) py(i)
 14     continue
        call pgsci(5)
        call pgpt(npt,px,py,sym)
        call pgsci(1)
      endif
      
      call pgsci(2)
      call pgpt1(0.75,real(mmin),18)
      call pgsci(1)
      
      npt2=npt
      call binp(npt2,phase,mag2,merr2,bins,flag)
      j=0
      do 12 i=1,npt2
        if(flag(i).eq.0)then
            j=j+1
            px(j)=real(phase(i))
            py(j)=real(mag2(i))
            pz(j)=real(merr2(i))
        endif
 12   continue
      npt2=j
c      write(0,*) "npt2:",npt2
      if(pan.eq.2)then
        call pgsci(3)
      else
        call pgsci(4)
      endif
      call pgsch(1.5)
      call pgpt(npt2,px,py,17)
      call pgsch(0.8)
      if(x1.gt.0.5) call pgerrb(6,npt2,px,py,pz,1.0)
      call pgsci(1)
      
      do 13 i=1,4
        rbb(i)=real(bbp(i))
 13   continue
      
      if(pan.eq.2)then
        call pgsch(0.7)
                
        sn=transn(nfit,sol,npt,time,mag,merr,mag2)
c        write(0,*) "S/N:",sn
        write(label,513) "S/N=",sn
 513    format(A4,F7.1)
        if(nsensor.lt.2) call pgptxt(0.75,
     .      rbb(4)+0.05*(rbb(4)-rbb(3)),0.0,0.5,label)
        
        write(label,502) "Tc(E)=",sol(7),"(BJD)+E(",sol(5),"days)"
 502    format(A6,F8.4,A8,F10.5,A5)
c        write(0,*) label
        if(nsensor.lt.2) call pgptxt(rbb(1)+0.03*(rbb(2)-rbb(1)),
     .      rbb(4)+0.18*(rbb(4)-rbb(3)),0.0,0.0,label)
        write(label,504) "       ",err(7),"         ",err(5),"     "
 504    format(A7,F8.4,A9,F9.5,A5)
        if(nsensor.lt.2) call pgptxt(rbb(1)+0.03*(rbb(2)-rbb(1)),
     .      rbb(4)+0.05*(rbb(4)-rbb(3)),
     .      0.0,0.0,label)
      
        call pgsch(0.8)
      elseif(pan.eq.3)then
        call pgsch(0.7)
        
c        write(label,501) "Rp=",sol(4)/redj,"Re"
        write(label,501) "Rp=",sol(4)/redj*rad/sol(3),"Re" 
 501    format(A3,F6.2,1X,A2)
        if(nsensor.lt.2) call pgptxt(rbb(1)+0.01*(rbb(2)-rbb(1)),
     .          rbb(4)+0.04*(rbb(4)-rbb(3)),
     .          0.0,0.0,label)
     
        M1=sol(1)*Msun
        M2=sol(2)*Mjup
        R1=sol(3)*Rsun
        Psec=sol(5)*8.64d4
        incl=sol(6)
        asemi=(M1+M2)**(1.0e0/3.0e0)*Psec**(2.0d0/3.0d0)*aConst
        b=asemi/R1*tan(Pi*(90.0-incl)/180.0)
        write(label,514) "b=",b
 514    format(A2,F4.2)
        if(nsensor.lt.2) call pgptxt(rbb(1)+0.33*(rbb(2)-rbb(1)),
     .          rbb(4)+0.04*(rbb(4)-rbb(3)),
     .          0.0,0.0,label)
        
        if(err(16).gt.0.0d0)then
            if(sol(16)/err(16).lt.10.0)then
                write(label,505) "Eclip-sig ",sol(16)/err(16)
            else
                write(label,512) "Eclip-sig ",sol(16)/err(16)
            endif
 505        format(A10,F4.1)
 512        format(A10,F5.1)
            if(nsensor.lt.2) call pgptxt(rbb(1)-0.02*(rbb(2)-rbb(1)),
     .          rbb(4)+0.10*(rbb(4)-rbb(3)),
     .          0.0,0.0,label)
        endif
        
        if(sol(16)/err(16).ge.2.0d0)then !only plot if reasonable detect
            if(Thpars(1,1).lt.10.0)then
                write(label,506) "Ag= ",Thpars(1,1),"(",Thpars(1,2),")"
 506            format(A4,F4.2,A1,F4.2,A1)
            else
                write(label,509) "Ag=",Thpars(1,1),"(",Thpars(1,2),")"
 509            format(A3,F6.2,A1,F6.2,A1)
            endif
            if(nsensor.lt.2) call pgptxt(rbb(1)+0.28*(rbb(2)-rbb(1)),
     .          rbb(4)+0.10*(rbb(4)-rbb(3)),
     .          0.0,0.0,label)
            if(Thpars(3,1).lt.1.0d4)then
                write(label,508) "Teff=",int(Thpars(3,1)+0.5),"(",
     .              int(Thpars(3,2)+0.5),") K"
 508            format(A5,I4,A1,I4,A3)
            else
                write(label,511) "Teff=",int(Thpars(3,1)+0.5),"(",
     .              int(Thpars(3,2)+0.5),") K"
 511            format(A5,I5,A1,I5,A3)
            endif
            if(nsensor.lt.2) call pgptxt(rbb(1)+0.66*(rbb(2)-rbb(1)),
     .          rbb(4)+0.10*(rbb(4)-rbb(3)),
     .          0.0,0.0,label)
        endif        
        if(Thpars(2,1).lt.1.0d4)then
            write(label,507) "Teq=",int(Thpars(2,1)+0.5),"(",
     .          int(Thpars(2,2)+0.5),") K"
 507        format(A4,I4,A1,I4,A3)
        else
            write(label,510) "Teq=",int(Thpars(2,1)+0.5),"(",
     .          int(Thpars(2,2)+0.5),") K"
 510        format(A4,I5,A1,I5,A3)
        endif
        if(nsensor.lt.2) call pgptxt(rbb(1)+0.66*(rbb(2)-rbb(1)),
     .      rbb(4)+0.04*(rbb(4)-rbb(3)),
     .      0.0,0.0,label)
        call pgsch(0.8)
      endif
 
      if(x1.ge.0.5)then
        width=(x2-x1)*period*24.0d0/2.0d0
        call pgwindow(real(-width),real(width),real(ax1),real(ax2))
      else
        call pgwindow(real(x1),real(x2),real(ax1),real(ax2))
      endif
      call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0) 
c      call pgbox('BCNTS1',0.0,0,'BCTSV1',0.0,0) 

      return
      end
