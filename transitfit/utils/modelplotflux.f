CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine modelplotflux(nfit,sol,nmax,phase,bb,bbp,bbp2,bbp3,
     .  avgitime,avgvtime,toff)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
C     nfit: the number parameters for the model
      integer nfit,nplot,i,nmax
      parameter(nplot=10000) !number of points for plotting
C     sol: parameters for the transitmodel
C     time: time for each model point
C     tmodel: the model at each time point
C     nmax,px,py,phase: workspace variables
C     bb,bbp: plotting boundies from time-series and phase plot
C     avgitime: average integration time
      integer dtype(nplot)
      real px(nplot),py(nplot)
      double precision sol(nfit),time(nplot),exptime(nplot),dnpt,
     .  tmodel(nplot),phase(nmax),bb(4),bbp(4),bbp2(4),bbp3(4),avgitime,
     .  toff,avgvtime
 
      dnpt=dble(nplot-1)
      do 10 i=1,nplot  !generate 1 phase of model data
C     Psec contains the orbital Period of the planet in seconds
        time(i)=bb(1)+(bb(2)-bb(1))*dble(i-1)/dnpt
        exptime(i)=avgitime
        dtype(i)=0
 10   continue
      call pgsci(2)
 
c      goto 20 !only do panel 3 
      
c      call transitmodel(nplot,time,exptime,dtype,tmodel,nfit,sol)
      
c      call panel(1)
c      call pgwindow(real(bb(1)),real(bb(2)),real(bb(3)),real(bb(4)))
c      do 11 i=1,nplot
c        px(i)=real(time(i))
c        py(i)=real(tmodel(i)*1000.0)
cc        write(6,*) px(i),py(i)
c 11   continue

c      call pgsci(2)
c      call pgline(nplot,px,py)
c      call pgsci(1)
      
c      open(unit=13,file="tmodel.dat")
c      do 16 i=1,nplot
c        write(13,*) time(i),10**(tmodel(i)/-2.5)
c 16   continue
c      write(13,*) "# PHASE"
      
      call panel(2)
      
      do 13 i=1,nplot  !generate 1 phase of model data
C     Psec contains the orbital Period of the planet in seconds
c        time(i)=bb(1)+(sol(5)-bb(1))*dble(i-1)/dnpt
        time(i)=bb(1)+sol(5)*dble(i-1)/dnpt
        exptime(i)=avgitime
        dtype(i)=0
 13   continue
  
      call transitmodel(nplot,time,exptime,dtype,tmodel,nfit,sol)
      
      call pgwindow(real(bbp(1)),real(bbp(2)),real(bbp(3)),real(bbp(4)))
      call phasept(nplot,time,phase,sol(5),toff)
      call sort2(nplot,phase,tmodel)
      do 12 i=1,nplot
        px(i)=real(phase(i))
        py(i)=real(10.0**(tmodel(i)/-2.5))
 12   continue
      call pgline(nplot,px,py,-1)
      
c      do 17 i=1,nplot
c        write(13,*) phase(i),10**(tmodel(i)/-2.5)
c 17   continue
c      close(13)

c      goto 999

 20   call panel(3)
      
      do 18 i=1,nplot  !generate 1 phase of model data
C     Psec contains the orbital Period of the planet in seconds
c        time(i)=bb(1)+(sol(5)-bb(1))*dble(i-1)/dnpt
        time(i)=sol(7)-sol(5)*(0.75-bbp2(1))+
     .      sol(5)*dble(i-1)/dnpt*(bbp2(2)-bbp2(1))
c        write(6,*) time(i),sol(7)
c        read(5,*) 
        exptime(i)=avgitime
        dtype(i)=0
 18   continue
  
      call transitmodel(nplot,time,exptime,dtype,tmodel,nfit,sol)
      call phasept(nplot,time,phase,sol(5),toff)
      call sort2(nplot,phase,tmodel)
      do 19 i=1,nplot
        px(i)=real(phase(i))
        py(i)=real(10.0**(tmodel(i)/-2.5))
 19   continue
      
      call pgwindow(real(bbp2(1)),real(bbp2(2)),real(bbp2(3)),
     .  real(bbp2(4)))
      call pgsch(3.0)
      call pgline(nplot,px,py,-1)
      call pgsch(1.0)
      
c      do 15 i=1,nplot  !generate 1 phase of model data
c        time(i)=bb(1)+sol(5)*dble(i-1)/dnpt
c        exptime(i)=avgvtime
c        dtype(i)=1
c 15   continue
      
c      call transitmodel(nplot,time,exptime,dtype,tmodel,nfit,sol)
c      
c      call panel(4)
c      call pgwindow(real(bbp3(1)),real(bbp3(2)),real(bbp3(3)),
c     .  real(bbp3(4)))
c      call phasept(nplot,time,phase,sol(5),toff)
c      call sort2(nplot,phase,tmodel)
c
cc      open(unit=20,file="rvphase.dat")
c      do 14 i=1,nplot
c        px(i)=real(phase(i))
c        py(i)=real(tmodel(i))
cc        write(20,*) px(i),py(i)
c 14   continue
cc      close(20)
c      call pgline(nplot,px,py,-1)
      


 999  call pgsci(1)
       
      
      return
      end
      
