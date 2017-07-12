      program plot975
C     Plots for KOI975 paper
      implicit none
      integer nmodelmax,nmodel,nmax,i,nhx
      parameter(nmodelmax=1000,nmax=100000)
      real px,py
      double precision massi,Z,tage(nmodelmax),tTeff(nmodelmax),
     .  tlogL(nmodelmax),trad(nmodelmax),trho(nmodelmax),
     .  tdrhodt(nmodelmax),dumr,hxL(nmax),hxTeff(nmax),hxchisq(nmax),
     .  chimin,hxmass(nmax),hxrad(nmax)
C     Yonsei Yale Model Variables
      integer inda,indz,indm,nline
      parameter (nline=150)
      parameter (inda=3)
      parameter (indz=11)
      parameter (indm=35)
      integer a_one,z_one,m_one
      real*8 avalue(inda)
      real*8 zvalue(indz)
      real*8 xmass(indm)
      common /grid/avalue,zvalue,xmass
      common /single/a_one,z_one,m_one
 
      open(unit=10,file='koi975.all.hmx',status='old',err=901)

      i=1
 10   read(10,*,end=11) hxmass(i),hxrad(i),dumr,dumr,hxTeff(i),hxL(i),
     .  dumr,dumr,hxchisq(i)
        if(hxmass(i).gt.4.0d0) goto 10
        i=i+1
        goto 10
 11   continue
      nhx=i-1

      close(10)
      
      chimin=hxchisq(1)
      do 13 i=2,nhx
        chimin=min(hxchisq(i),chimin)
 13   continue

      call pgopen('?')
      call pgpage()
      call PGPAP ( 6.0 ,1.0)  
      call pgvport(0.15,0.85,0.2,0.9)
c      call PGENV ( -0.5 , 0.5 , -0.5 , 0.5 , 1 , -2 )  
      call pgwindow(15000.0,2000.0,-2.0,4.0)
      call pgslw(3)
      call pgsch(1.5)
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("T\deff\u (K)",
     .  "log(L/L\d\(2281)\u)","")
 
c      call pgsci(1)
      massi=0.6
      Z=0.02
      do while(massi.lt.4.2)
        call gettrack(massi,Z,nmodelmax,nmodel,tage,tTeff,tlogL,
     .      trad,trho,tdrhodt)
        call pgbbuf()
        do 17 i=1,nmodel
c            if(tage(i).lt.12.0d0)then
                px=tTeff(i)
                py=tlogL(i)
                call pgpt1(px,py,-1)
c            endif
 17     continue
        call pgebuf()
        massi=massi+0.2
      enddo      
      call pgsci(1)
     
c      call pgsci(2)
c      call pgbbuf()
c      do 12 i=1,nhx
c        px=real(hxTeff(i))
c        py=real(hxL(i))
c        call pgpt1(px,py,-1)
c 12   continue
c      call pgebuf()
c      call pgsci(1)
      
      call pgsci(3)
      call pgbbuf()
      do 14 i=1,nhx
        if(hxchisq(i).lt.chimin+3)then
            px=real(hxTeff(i))
            py=real(hxL(i))
            call pgpt1(px,py,17)
        endif
 14   continue
      call pgebuf()
      call pgsci(1)
     
c      call pgsci(4)
c      call pgbbuf()
c      do 15 i=1,nhx
c        if(hxchisq(i).lt.chimin+2)then
c            px=real(hxTeff(i))
c            py=real(hxL(i))
c            call pgpt1(px,py,17)
c        endif
c 15   continue
c      call pgebuf()
c      call pgsci(1)
      
c      call pgsci(5)
c      call pgbbuf()
c      do 16 i=1,nhx
c        if(hxchisq(i).lt.chimin+1)then
c            px=real(hxTeff(i))
c            py=real(hxL(i))
c            call pgpt1(px,py,17)
c        endif
c 16   continue
c      call pgebuf()
c      call pgsci(1)
     
      call pgpage()
      
      call pgwindow(0.4,4.0,0.5,3.0)
      call pgslw(3)
      call pgsch(1.5)
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Mass (M\d\(2281)\u)",
     .  "Radius (R\d\(2281)\u)","")
          
      massi=0.6
      Z=0.02
      do while(massi.lt.4.2)
        call gettrack(massi,Z,nmodelmax,nmodel,tage,tTeff,tlogL,
     .      trad,trho,tdrhodt)
        call pgbbuf()
        do 18 i=1,nmodel
c            if(tage(i).lt.12.0d0)then
                px=massi
                py=trad(i)
                call pgpt1(px,py,-1)
c            endif
 18     continue
        call pgebuf()
        massi=massi+0.2
      enddo      
      call pgsci(1)
      
      call pgsci(3)
      call pgbbuf()
      do 19 i=1,nhx
        if(hxchisq(i).lt.chimin+3)then
            px=real(hxmass(i))
            py=real(hxrad(i))
            call pgpt1(px,py,17)
        endif
 19   continue
      call pgebuf()
      call pgsci(1)     
     
      call pgclos()

  
      goto 999
 901  write(0,*) "Cannot open koi975.all.hmx"
      goto 999
 999  end