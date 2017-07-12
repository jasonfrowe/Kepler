CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotsetup(rbb,nfit,nplanet,sol,which)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer CI1,CI2,P,col,nfit,i,nplanet,which
      real rbb(3,4),sign,contra,bright
      double precision per,adrs,eccn,aph,Pi,G,sol(nfit)

      Pi=acos(-1.d0)!define Pi
      G=6.674d-11 

      if((which.eq.0).or.(which.eq.1))then
        call pgpanl(1,1)
        call pgeras()
        rbb(1,1)=-1.1
        rbb(1,2)= 1.1
        rbb(1,3)=-1.1
        rbb(1,4)= 1.1
        call pgvport(0.15,0.85,0.2,0.9)
        call pgwindow(rbb(1,1),rbb(1,2),rbb(1,3),rbb(1,4))
        call pgslw(3)
        call pgsch(1.5)
        CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
        call pglabel('r/R*','r/R*','')
      endif

      CI1 = 16 
      CI2 = 86
      call PGSCIR ( CI1 , CI2 ) 
      P      = 3
      SIGN   = 1.0
      CONTRA = 1.0!0.56
      BRIGHT = 0.5  
      call Palett ( P , SIGN*CONTRA , BRIGHT )

      if((which.eq.0).or.(which.eq.2))then      
        call pgpanl(2,1)
        call pgeras()
        col=11*(1-1)+15
        per=sol(col+2) !period in days
        do 10 i=2,nplanet
            col=11*(i-1)+15
            per=max(per,sol(col+2)) !period in days
 10     continue
        adrs=1000.0*sol(1)*G*(Per*86400.0d0)**2/(3.0d0*Pi)
        adrs=adrs**(1.0d0/3.0d0) !a/R*
        eccn=sqrt(sol(col+5)*sol(col+5)+sol(col+6)*sol(col+6))
        if(eccn.ge.1.0) eccn=0.99
        aph=adrs*(1.0-eccn*eccn)/(1.0-eccn)
c      write(0,*) 'adrs:',adrs,col
        rbb(2,1)=real(-aph-aph*0.1)
        rbb(2,2)=real( aph+aph*0.1)
        rbb(2,3)=real(-aph-aph*0.1)
        rbb(2,4)=real( aph+aph*0.1)  
        call pgvport(0.15,0.85,0.2,0.9)
        call pgwindow(rbb(2,1),rbb(2,2),rbb(2,3),rbb(2,4))
        CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
        call pglabel('r/R*','r/R*','')
        call pgsci(ci2)
        call pgcirc(0,0,1.0)
        call pgsci(1)
      endif
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine plotphase(nfit,np,sol,npt,time,flux,ferr,dtype,tmodel,
     .  rbb)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,np,npt,col,i,dtype(npt)
      real rbb(3,4)
      double precision sol(nfit),time(npt),flux(npt),tmodel(npt),
     .  transitdur,tdur,toff,ph1,ph2,phase,mmin,mmax,ferr(npt)
     
c      write(0,*) "npt:",npt
      mmin=flux(1)
      mmax=flux(1)
      do 10 i=2,npt
        if(dtype(i).eq.0)then
            mmin=min(mmin,flux(i))
            mmax=max(mmax,flux(i))
        endif
c        write(0,*) mmin,mmax,flux(i)
 10   continue
        
      rbb(3,3)=real(mmin-(mmax-mmin)*0.1)
      rbb(3,4)=real(mmax+(mmax-mmin)*0.1)  
        
        
      col=11*(np-1)+15
      
      tdur=transitdur(np,nfit,sol)/(24.0*60.0*60.0)
c      write(0,*) "Tdur: ",tdur

      toff=0.5-(sol(1+col)/sol(2+col)-int(sol(1+col)/sol(2+col)))
      if(toff.lt.0.0)toff=toff+1.0
      ph1=0.5-1.0d0*tdur/sol(2+col)
      if(ph1.lt.0.25)ph1=0.25
      ph2=0.5+1.0d0*tdur/sol(2+col)
      if(ph2.gt.0.75)ph2=0.75
        
      rbb(3,1)=real(ph1)
      rbb(3,2)=real(ph2)     
c      write(0,*) toff,ph1,ph2
      
      call pgpanl(3,1)
      call pgeras()
      call pgvport(0.15,0.85,0.2,0.9)
      call pgwindow(rbb(3,1),rbb(3,2),rbb(3,3),rbb(3,4))
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel('Phase','Relative Flux','')
      
      call pgbbuf()
      do 11 i=1,npt

        if(dtype(i).eq.0)then
            phase=time(i)/sol(2+col)-int(time(i)/sol(2+col))
            phase=phase+toff
            if(phase.lt.0.0) phase=phase+1.0
            if(phase.gt.1.0) phase=phase-1.0
            if((phase.gt.ph1).and.(phase.lt.ph2))then
                call pgpt1(real(phase),real(flux(i)),-1)
                call pgsci(2)
                call pgpt1(real(phase),real(tmodel(i)),-1)
                call pgsci(1)
c            write(6,*) time(i)+54900.0-0.5,flux(i)-1.0d0,ferr(i)
            endif
        endif
        
 11   continue     
      call pgebuf()
      
      return
      end
      
c     --------------------------------------------
      SUBROUTINE PALETT ( TYPE , CONTRA , BRIGHT )
c     --------------------------------------------
c
c     Sets a "palette" of colors in the range of defined color indices 
c     This subroutine is distributed with PGPLOT in one of the demos
c
      INTEGER TYPE
      REAL CONTRA, BRIGHT
C
      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL CL(5), CR(5), CG(5), CB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)
C
      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
C
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
C
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
C
      DATA CL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA CR /0.0, 0.0, 0.0, 0.3, 1.0/
      DATA CG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA CB /0.0, 0.5, 1.0, 1.0, 1.0/
C
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
C
      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
     :         0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
     :         0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
     :         0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
     :         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
C
      IF (TYPE.EQ.1) THEN
C        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.2) THEN
C        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.3) THEN
C        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.4) THEN
C        -- freeze
         CALL PGCTAB(CL, CR, CG, CB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.5) THEN
C        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
      END IF
      END
 