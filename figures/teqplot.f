      program teqplot
      implicit none
      integer npt,nmax,i
      parameter(nmax=5000)
      integer flag(nmax),sym,npri(nmax)
      real dumr,rad(nmax),Teq(nmax),RjdRs,px(2),py(2),Rs,Rl,Ts,Tl,
     .  koi(nmax)
      character*80 filename,dumc
      
      RjdRs=10.863 !Size of Jupiter relative to Earth
      Rs=0.4 !smaller radius to plot
      Rl=20.0  !largest radius to plot
      Ts=0.0   !Temperature range to plot
      Tl=2000.0
      
      filename="Teq.dat"
      open(unit=10,file=filename,status='old',err=901)
      read(10,*) dumc
      read(10,*) dumc
      
      i=0
 10   read(10,*,end=11) koi(i),dumr,dumr,dumr,rad(i),dumr,dumr,dumr,
     .  Teq(i),flag(i),npri(i)
        if (Teq(i).le.0.0) goto 10
c        write(6,*) Teq(i)
        rad(i)=log10(rad(i)*RjdRs)
c        Teq(i)=log10(Teq(i))
        i=i+1
        goto 10
 11   continue
      npt=i-1
      
      write(0,*) "Points read: ",i
      
      call pgopen('?')
      call PGPAP ( 7.0 ,1.0)  
c      call PGENV ( -0.5 , 0.5 , -0.5 , 0.5 , 1 , -2 )  
      call pgwindow(Ts,Tl,log10(Rs),log10(Rl))
      call pgslw(3)
      CALL PGBOX('BCNTS1',0.0,0,'BCLNTS1',0.0,0)
      call pglabel("T\deq\u (K)","Radius (R\d\(2284)\u)","")
      call pgbbuf()
      call pgsch(2.5)
      do 12 i=1,npt
c        sym=17          !sequestered data
c        call pgsci(4)
c        if(npri(i).eq.4) sym=18
c        if(flag(i).eq.1)then  !single transits
c            sym=21
cc            if(npri(i).eq.4) sym=0
c            call pgsci(1)
c        elseif(flag(i).eq.2)then  !false-positive
c            sym=1
c            call pgsci(1)
c            goto 12
c        elseif(flag(i).eq.3)then  !released data
c            sym=17
c            if(npri(i).eq.4) sym=18
c            call pgsci(2)
c        elseif(flag(i).eq.5)then !new candidates
c            sym=17
c            call pgsci(3)
        if((koi(i).lt.964.0).and.(flag(i).eq.0))then
            sym=17
            call pgsci(2)
            call pgpt1(teq(i),rad(i),sym)
        elseif((koi(i).gt.964.0).and.(flag(i).eq.5))then
            sym=17
            call pgsci(3)
            call pgpt1(teq(i),rad(i),sym)
        elseif((koi(i).gt.1611.0).and.(flag(i).eq.6))then
            sym=17
            call pgsci(4)
            call pgpt1(teq(i),rad(i),sym)
        endif        
 12   continue
      call pgsch(1)
      call pgsci(1)
      call pgebuf()
      call pgsch(1.5) !plot the Earth
      call pgpt1(255.0,log10(1.0),12)
c      call pgsch(1.0)
      py(1)=log10(Rs)
      py(2)=log10(Rl)
      px(1)=273.0
      px(2)=273.0
      call pgline(2,px,py)
      px(1)=373.0
      px(2)=373.0
      call pgline(2,px,py)
      px(1)=Ts
      px(2)=Tl
      py(1)=0.0 !Earth sized
      py(2)=0.0
      call pgline(2,px,py)
      py(1)=log10(1.25) !Earth sized
      py(2)=log10(1.25)
      call pgline(2,px,py)
      call pgclos()
      
      close(10)
      
      goto 999
 901  write(0,*) "Cannot open ",filename
      goto 999
 999  end