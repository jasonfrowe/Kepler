CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine explore(nexpl,npt,aT,aM,aE,aIT,dtype,tmodel,nfit,sol,
     .          sol2,serr,Dpvary)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nexpl,npt,nfit,nplot,i,nans
      integer dtype(npt)
      parameter(nplot=100)
      real px(nplot),py(nplot)
      double precision aT(npt),aM(npt),aE(npt),aIT(npt),
     .  tmodel(npt),sol(nfit),Dpvary(nfit),x1,x2,chi(nplot),sol2(nfit),
     .  dnplotm1,chimin,chimax,dnptm1,vary,solmin,dsol,serr(nfit,2),
     .  likelihood
      logical loop
      include "titles.f"
     
C     Select plotting panel
      call panel(4)
      
      if(Dpvary(nexpl).eq.0.0d0)then
        write(0,*) "Enter exploration length"
        read(5,*) vary
      else
        vary=Dpvary(nexpl)
      endif
      loop=.true.
      do while(loop)
      
C       Get the range to explore chi-sq      
        x1=sol(nexpl)-vary
        x2=sol(nexpl)+vary

C     make a copy of array      
        do 5 i=1,nfit
            sol2(i)=sol(i)
 5      continue
      
        dnptm1=dble(npt-1)
        dnplotm1=dble(nplot-1)
        sol2(nexpl)=x1
c        chi(1)=chisquared(npt,tmodel,time,mag,merr,itime,nfit,sol2)
c     .      /dnptm1
        chi(1)=likelihood(npt,tmodel,aT,aM,aE,aIT,dtype,nfit,sol2,serr)
        chimin=chi(1)
        chimax=chi(1)
        px(1)=real(x1)
        py(1)=real(chi(1))
        do 10 i=2,nplot
            sol2(nexpl)=x1+(x2-x1)*dble(i-1)/dnplotm1
c            chi(i)=chisquared(npt,tmodel,time,mag,merr,itime,nfit,sol2)
c     .          /dnptm1
             chi(i)=likelihood(npt,tmodel,aT,aM,aE,aIT,dtype,nfit,sol2,
     .          serr)
c        write(0,*) sol2(nexpl),chi(i)
c        read(5,*)
            if(chi(i).lt.chimin) then
                chimin=chi(i)
                solmin=sol2(nexpl)
            endif
            chimax=max(chi(i),chimax)
            px(i)=real(sol2(nexpl))
            py(i)=real(chi(i))
 10     continue
        dsol=sol(nexpl)-solmin
        write(0,500) "Chimin: ",solmin,chimin
        write(0,500) "Oldsol: ",sol(nexpl),dsol
 500    format(A8,3(1X,1PE12.5))
 
        call pgwindow(real(x1),real(x2),real(chimin),real(chimax))
        call pgbox('BCNTS1',0.0,0,'BCNTS',0.0,0)
        call pgptxt(real((x1+x2)/2.0),real(chimin+0.12*(chimin-chimax)),
     .      0.0,0.5,titles(nexpl))
        call pgptxt(real(x1-0.1*(x2-x1)),real((chimin+chimax)/2),
     .      90.0,0.5,"Chi-Sq")
        call pgline(nplot,px,py)
      
        write(0,*) "Accept new Dpvary (0), Explore (1), Exit (2)"
        read(5,*) nans
        if(nans.eq.0) then
            Dpvary(nexpl)=-dsol
            loop=.false.
        elseif(nans.eq.1)then
            call pgsci(0)
            call pgbox('BCNTS1',0.0,0,'BCNTS',0.0,0)
            call pgptxt(real((x1+x2)/2.0),
     .          real(chimin+0.12*(chimin-chimax)),0.0,0.5,titles(nexpl))
            call pgptxt(real(x1-0.1*(x2-x1)),real((chimin+chimax)/2),
     .          90.0,0.5,"Chi-Sq")
            call pgline(nplot,px,py)
            call pgsci(1)
            write(0,*) "Enter exploration length"
            read(5,*) vary
        else
            loop=.false.
        endif
        if(vary.eq.0.0d0) loop=.false.
      enddo
      
      
     
      return
      end