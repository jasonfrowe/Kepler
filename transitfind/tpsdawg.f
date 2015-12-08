      program tpsdawg
      implicit none
      integer nmax,iargc,nunit,i,npt,j,nplot,nPmatch,nPmatch2,nPmatch3
      parameter(nmax=10000)
      integer koi(nmax),z(nmax)
      real q3Per(nmax),q3Epo(nmax),q3MES(nmax),q3MDS(nmax),
     .  q6Per(nmax),q6Epo(nmax),q6MES(nmax),q6MDS(nmax),tpspars(5),
     .  q8Per(nmax),q8Epo(nmax),q8MES(nmax),q8MDS(nmax),x(nmax),y(nmax),
     .  dd(4),pcut,s2,ecut,px(2),py(2)
      character*80 filename      
      
      s2=sqrt(2.0)
      pcut=1.0e-3
      ecut=0.02
      
      do 12 i=1,nmax
        q3MES(i)=0.0d0
        q6MES(i)=0.0d0
        q8MES(i)=0.0d0
 12   continue
      
      if(iargc().lt.1) goto 901 !check number of command line arguements
      call getarg(1,filename) !get filename for input data
      nunit=10
      
      open(unit=nunit,file=filename,status='old',err=902)
      
      i=1
 10   read(nunit,*,end=11) koi(i)
        i=i+1
      goto 10
 11   continue
      npt=i-1
      write(6,*) "Number of KOIs: ",npt
      
      close(nunit)
      
      filename="tps_q1q3.dat"
      open(unit=nunit,file=filename,status='old',err=902)
      
 13   read(nunit,*,end=15) (tpspars(i),i=1,5)
        do 14 i=1,npt
            if(koi(i).eq.int(tpspars(1)))then
                if(q3MES(i).lt.tpspars(4))then
                    q3MES(i)=tpspars(4)
                    q3Per(i)=tpspars(2)
                    q3MDS(i)=tpspars(3)
                    q3Epo(i)=tpspars(5)
                endif
             endif
 14     continue
      goto 13
 15   continue
      close(nunit)
      
      j=0
      do 16 i=1,npt
        if((q3MES(i).gt.7.1).and.(q3MDS(i).gt.s2)) j=j+1
 16   continue
      write(6,*) "Q1Q3 KOIs",j,real(j)/real(npt)

      filename="tps_q1q6.dat"
      open(unit=nunit,file=filename,status='old',err=902)
      
 17   read(nunit,*,end=19) (tpspars(i),i=1,5)
        do 18 i=1,npt
            if(koi(i).eq.int(tpspars(1)))then
                if(q6MES(i).lt.tpspars(4))then
                    q6MES(i)=tpspars(4)
                    q6Per(i)=tpspars(2)
                    q6MDS(i)=tpspars(3)
                    q6Epo(i)=tpspars(5)
                endif
             endif
 18     continue
      goto 17
 19   continue
      close(nunit)
      
      j=0
      do 20 i=1,npt
        if((q6MES(i).gt.7.1).and.(q6MDS(i).gt.s2)) j=j+1
 20   continue
      write(6,*) "Q1Q6 KOIs",j,real(j)/real(npt)      
      
      filename="tps_q1q8.dat"
      open(unit=nunit,file=filename,status='old',err=902)
      
 21   read(nunit,*,end=23) (tpspars(i),i=1,5)
        do 22 i=1,npt
            if(koi(i).eq.int(tpspars(1)))then
                if(q8MES(i).lt.tpspars(4))then
                    q8MES(i)=tpspars(4)
                    q8Per(i)=tpspars(2)
                    q8MDS(i)=tpspars(3)
                    q8Epo(i)=tpspars(5)
                endif
             endif
 22     continue
      goto 21
 23   continue
      close(nunit)
      
      j=0
      do 24 i=1,npt
        if((q8MES(i).gt.7.1).and.(q8MDS(i).gt.s2))j=j+1
 24   continue
      write(6,*) "Q1Q8 KOIs",j,real(j)/real(npt)
      
      call pgopen('?')
      call pgpage()
      call PGPAP ( 6.0 ,1.0)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      dd(1)= 9.9e30
      dd(2)=-9.9e30
      dd(3)= 9.9e30
      dd(4)=-9.9e30
      j=0
      nPmatch=0
      nPmatch2=0
      nPmatch3=0
      do 25 i=1,npt
        if((q3MES(i).gt.0.0).and.(q6MES(i).gt.0.0))then
c        if((q3MES(i).gt.7.1).and.(q6MES(i).gt.7.1).and.
c     .    (q3MDS(i).gt.s2).and.(q6MDS(i).gt.s2))then
            j=j+1
            x(j)=real(j)
            y(j)=q3Per(i)/q6Per(i)-int(q3Per(i)/q6Per(i))
            if(y(j).lt.0.0) y(j)=y(j)+1.0
            if((y(j).lt.Pcut).or.(1.0-y(j).lt.Pcut)) nPmatch=nPmatch+1
            if(abs(y(j)-0.5).lt.Pcut) nPmatch2=nPmatch2+1
            if(abs(y(j)-1.0/3.0).lt.Pcut) nPmatch3=nPmatch3+1
            if((q3MES(i).gt.7.1).and.(q6MES(i).gt.7.1).and.
     .       (q3MDS(i).gt.s2).and.(q6MDS(i).gt.s2))then
                z(j)=1
            else
                z(j)=2
            endif
            dd(1)=min(dd(1),x(j))
            dd(2)=max(dd(2),x(j))
            dd(3)=min(dd(3),y(j))
            dd(4)=max(dd(4),y(j))
        endif
 25   continue
      nplot=j
      
      write(0,*) "Q1Q3, Q1Q6 KOI     Match: ",nplot
      write(0,*) "Q1Q3, Q1Q6 Period  Match: ",nPmatch,
     .  real(nPmatch)/real(nplot)
      write(0,*) "Q1Q3, Q1Q6 Period Double: ",nPmatch2,
     .  real(nPmatch2)/real(nplot)
      write(0,*) "Q1Q3, Q1Q6 Period Triple: ",nPmatch3,
     .  real(nPmatch3)/real(nplot)
      
      call pgwindow(dd(1),dd(2),
     .  dd(3)-0.1*(dd(4)-dd(3)),dd(4)+0.1*(dd(4)-dd(3)))
      call pgslw(3)
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Candidate #","P3/P6-int(P3/P6)","")
c      call pgpt(nplot,x,y,17)

      call pgbbuf()
      do 26 i=1,nplot
        call pgsci(z(i))
        call pgpt1(x(i),y(i),17)
 26   continue
      call pgebuf()
      call pgsci(1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call pgpage()

      dd(1)= 9.9e30
      dd(2)=-9.9e30
      dd(3)= 9.9e30
      dd(4)=-9.9e30
      j=0
      nPmatch=0
      nPmatch2=0
      nPmatch3=0
      do 27 i=1,npt
        if((q3MES(i).gt.0.0).and.(q8MES(i).gt.0.0))then
c        if((q3MES(i).gt.7.1).and.(q8MES(i).gt.7.1).and.
c     .    (q3MDS(i).gt.s2).and.(q8MDS(i).gt.s2))then
            j=j+1
            x(j)=real(j)
            y(j)=q3Per(i)/q8Per(i)-int(q3Per(i)/q8Per(i))
            if(y(j).lt.0.0) y(j)=y(j)+1.0
            if((y(j).lt.Pcut).or.(1.0-y(j).lt.Pcut)) nPmatch=nPmatch+1
            if(abs(y(j)-0.5).lt.Pcut) nPmatch2=nPmatch2+1
            if(abs(y(j)-1.0/3.0).lt.Pcut) nPmatch3=nPmatch3+1
            if((q3MES(i).gt.7.1).and.(q8MES(i).gt.7.1).and.
     .       (q3MDS(i).gt.s2).and.(q8MDS(i).gt.s2))then
                z(j)=1
            else
                z(j)=2
            endif
            dd(1)=min(dd(1),x(j))
            dd(2)=max(dd(2),x(j))
            dd(3)=min(dd(3),y(j))
            dd(4)=max(dd(4),y(j))
        endif
 27   continue
      nplot=j
      
      write(0,*) "Q1Q3, Q1Q8 KOI     Match: ",nplot
      write(0,*) "Q1Q3, Q1Q8 Period  Match: ",nPmatch,
     .  real(nPmatch)/real(nplot)
      write(0,*) "Q1Q3, Q1Q8 Period Double: ",nPmatch2,
     .  real(nPmatch2)/real(nplot)
      write(0,*) "Q1Q3, Q1Q8 Period Triple: ",nPmatch3,
     .  real(nPmatch3)/real(nplot)
      
      call pgwindow(dd(1),dd(2),
     .  dd(3)-0.1*(dd(4)-dd(3)),dd(4)+0.1*(dd(4)-dd(3)))
      call pgslw(3)
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Candidate #","P3/P8-int(P3/P8)","")
c      call pgpt(nplot,x,y,17)

      call pgbbuf()
      do 28 i=1,nplot
        call pgsci(z(i))
        call pgpt1(x(i),y(i),17)
 28   continue
      call pgebuf()
      call pgsci(1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call pgpage()

      dd(1)= 9.9e30
      dd(2)=-9.9e30
      dd(3)= 9.9e30
      dd(4)=-9.9e30
      j=0
      nPmatch=0
      nPmatch2=0
      nPmatch3=0
      do 29 i=1,npt
        if((q6MES(i).gt.0.0).and.(q8MES(i).gt.0.0))then
c        if((q6MES(i).gt.7.1).and.(q8MES(i).gt.7.1).and.
c     .    (q6MDS(i).gt.s2).and.(q8MDS(i).gt.s2))then
            j=j+1
            x(j)=real(j)
            y(j)=q6Per(i)/q8Per(i)-int(q6Per(i)/q8Per(i))
            if(y(j).lt.0.0) y(j)=y(j)+1.0
            if((y(j).lt.Pcut).or.(1.0-y(j).lt.Pcut)) nPmatch=nPmatch+1
            if(abs(y(j)-0.5).lt.Pcut) nPmatch2=nPmatch2+1
            if(abs(y(j)-1.0/3.0).lt.Pcut) nPmatch3=nPmatch3+1
            if((q6MES(i).gt.7.1).and.(q8MES(i).gt.7.1).and.
     .       (q6MDS(i).gt.s2).and.(q8MDS(i).gt.s2))then
                z(j)=1
            else
                z(j)=2
            endif
            dd(1)=min(dd(1),x(j))
            dd(2)=max(dd(2),x(j))
            dd(3)=min(dd(3),y(j))
            dd(4)=max(dd(4),y(j))
        endif
 29   continue
      nplot=j
      
      write(0,*) "Q1Q6, Q1Q8 KOI     Match: ",nplot
      write(0,*) "Q1Q6, Q1Q8 Period  Match: ",nPmatch,
     .  real(nPmatch)/real(nplot)
      write(0,*) "Q1Q6, Q1Q8 Period Double: ",nPmatch2,
     .  real(nPmatch2)/real(nplot)
      write(0,*) "Q1Q6, Q1Q8 Period Triple: ",nPmatch3,
     .  real(nPmatch3)/real(nplot)
      
      call pgwindow(dd(1),dd(2),
     .  dd(3)-0.1*(dd(4)-dd(3)),dd(4)+0.1*(dd(4)-dd(3)))
      call pgslw(3)
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Candidate #","P6/P8-int(P6/P8)","")
c      call pgpt(nplot,x,y,17)

      call pgbbuf()
      do 30 i=1,nplot
        call pgsci(z(i))
        call pgpt1(x(i),y(i),17)
 30   continue
      call pgebuf()
      call pgsci(1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call pgpage()

      dd(1)= 9.9e30
      dd(2)=-9.9e30
      dd(3)= 9.9e30
      dd(4)=-9.9e30
      j=0
      nPmatch=0
c      nPmatch2=0
c      nPmatch3=0
      do 31 i=1,npt
        if((q3MES(i).gt.0.0).and.(q6MES(i).gt.0.0))then
c        if((q6MES(i).gt.7.1).and.(q8MES(i).gt.7.1).and.
c     .    (q6MDS(i).gt.s2).and.(q8MDS(i).gt.s2))then
            j=j+1
            x(j)=real(j)
            y(j)=(q3Epo(i)-q6Epo(i))/q6Per(i)-
     .          int((q3Epo(i)-q6Epo(i))/q6Per(i))
            if((q3MES(i).gt.7.1).and.(q6MES(i).gt.7.1).and.
     .       (q3MDS(i).gt.s2).and.(q6MDS(i).gt.s2))then
                z(j)=1
                nPmatch=nPmatch+1
            else
                z(j)=2
            endif
            dd(1)=min(dd(1),x(j))
            dd(2)=max(dd(2),x(j))
            dd(3)=min(dd(3),y(j))
            dd(4)=max(dd(4),y(j))
        endif
 31   continue
      nplot=j
      
      call pgwindow(dd(1),dd(2),
     .  dd(3)-0.1*(dd(4)-dd(3)),dd(4)+0.1*(dd(4)-dd(3)))
      call pgslw(3)
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Candidate #","(E3-E6)/P6-int((E3-E6)/P6)","")
c      call pgpt(nplot,x,y,17)

      call pgbbuf()
      do 32 i=1,nplot
        call pgsci(z(i))
        call pgpt1(x(i),y(i),17)
 32   continue
      call pgebuf()
      call pgsci(1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call pgpage()

      dd(1)= 9.9e30
      dd(2)=-9.9e30
      dd(3)= 9.9e30
      dd(4)=-9.9e30
      j=0
      nPmatch=0
c      nPmatch2=0
c      nPmatch3=0
      do 33 i=1,npt
        if((q3MES(i).gt.0.0).and.(q8MES(i).gt.0.0))then
c        if((q6MES(i).gt.7.1).and.(q8MES(i).gt.7.1).and.
c     .    (q6MDS(i).gt.s2).and.(q8MDS(i).gt.s2))then
            j=j+1
            x(j)=real(j)
            y(j)=(q3Epo(i)-q8Epo(i))/q8Per(i)-
     .          int((q3Epo(i)-q8Epo(i))/q8Per(i))
            if((q3MES(i).gt.7.1).and.(q8MES(i).gt.7.1).and.
     .       (q3MDS(i).gt.s2).and.(q8MDS(i).gt.s2))then
                z(j)=1
                nPmatch=nPmatch+1
            else
                z(j)=2
            endif
            dd(1)=min(dd(1),x(j))
            dd(2)=max(dd(2),x(j))
            dd(3)=min(dd(3),y(j))
            dd(4)=max(dd(4),y(j))
        endif
 33   continue
      nplot=j
      
      call pgwindow(dd(1),dd(2),
     .  dd(3)-0.1*(dd(4)-dd(3)),dd(4)+0.1*(dd(4)-dd(3)))
      call pgslw(3)
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Candidate #","(E3-E8)/P8-int((E3-E8)/P8)","")
c      call pgpt(nplot,x,y,17)

      call pgbbuf()
      do 34 i=1,nplot
        call pgsci(z(i))
        call pgpt1(x(i),y(i),17)
 34   continue
      call pgebuf()
      call pgsci(1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call pgpage()

      dd(1)= 9.9e30
      dd(2)=-9.9e30
      dd(3)= 9.9e30
      dd(4)=-9.9e30
      j=0
      nPmatch=0
c      nPmatch2=0
c      nPmatch3=0
      do 35 i=1,npt
        if((q6MES(i).gt.0.0).and.(q8MES(i).gt.0.0))then
c        if((q6MES(i).gt.7.1).and.(q8MES(i).gt.7.1).and.
c     .    (q6MDS(i).gt.s2).and.(q8MDS(i).gt.s2))then
            j=j+1
            x(j)=real(j)
            y(j)=(q6Epo(i)-q8Epo(i))/q8Per(i)-
     .          int((q6Epo(i)-q8Epo(i))/q8Per(i))
            if((q6MES(i).gt.7.1).and.(q8MES(i).gt.7.1).and.
     .       (q6MDS(i).gt.s2).and.(q8MDS(i).gt.s2))then
                z(j)=1
                nPmatch=nPmatch+1
            else
                z(j)=2
            endif
            dd(1)=min(dd(1),x(j))
            dd(2)=max(dd(2),x(j))
            dd(3)=min(dd(3),y(j))
            dd(4)=max(dd(4),y(j))
        endif
 35   continue
      nplot=j
      
      call pgwindow(dd(1),dd(2),
     .  dd(3)-0.1*(dd(4)-dd(3)),dd(4)+0.1*(dd(4)-dd(3)))
      call pgslw(3)
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Candidate #","(E6-E8)/P8-int((E6-E8)/P8)","")
c      call pgpt(nplot,x,y,17)

      call pgbbuf()
      do 36 i=1,nplot
        call pgsci(z(i))
        call pgpt1(x(i),y(i),17)
 36   continue
      call pgebuf()
      call pgsci(1)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call pgpage()

      dd(1)= 9.9e30
      dd(2)=-9.9e30
      dd(3)= 9.9e30
      dd(4)=-9.9e30

      j=0
      do 37 i=1,npt
        if((q3MES(i).gt.0.0).and.(q6MES(i).gt.0.0).and.
     .    (q6Per(i).gt.0.0))then
            j=j+1
            x(j)=log10(q6Per(i))
c            write(0,*) q6Per(i),x(j)
            y(j)=q6MES(i)/q3MES(i)
            if((q3MES(i).gt.7.1).and.(q6MES(i).gt.7.1).and.
     .       (q3MDS(i).gt.s2).and.(q6MDS(i).gt.s2))then
                z(j)=1
                nPmatch=nPmatch+1
            else
                z(j)=2
            endif
            dd(1)=min(dd(1),x(j))
            dd(2)=max(dd(2),x(j))
            dd(3)=min(dd(3),y(j))
            dd(4)=max(dd(4),y(j))
        endif
 37   continue
      nplot=j
      dd(4)=3.0
      
      call pgwindow(dd(1),dd(2),
     .  dd(3)-0.1*(dd(4)-dd(3)),dd(4)+0.1*(dd(4)-dd(3)))
      call pgslw(3)
      CALL PGBOX('BCLNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Period(d)","Q6MES/Q3MES","")
c      call pgpt(nplot,x,y,17)

      call pgbbuf()
      do 38 i=1,nplot
        call pgsci(z(i))
        call pgpt1(x(i),y(i),17)
 38   continue
      call pgebuf()
      px(1)=dd(1)
      px(2)=dd(2)
      py(1)=1.51 !ratio of data
      py(2)=1.51
      call pgsci(4)
      call pgline(2,px,py)

      call pgsci(1)
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call pgpage()

      dd(1)= 9.9e30
      dd(2)=-9.9e30
      dd(3)= 9.9e30
      dd(4)=-9.9e30

      j=0
      do 39 i=1,npt
        if((q3MES(i).gt.0.0).and.(q8MES(i).gt.0.0).and.
     .    (q8Per(i).gt.0.0))then
            j=j+1
            x(j)=log10(q8Per(i))
c            write(0,*) q6Per(i),x(j)
            y(j)=q8MES(i)/q3MES(i)
            if((q3MES(i).gt.7.1).and.(q8MES(i).gt.7.1).and.
     .       (q3MDS(i).gt.s2).and.(q8MDS(i).gt.s2))then
                z(j)=1
                nPmatch=nPmatch+1
            else
                z(j)=2
            endif
            dd(1)=min(dd(1),x(j))
            dd(2)=max(dd(2),x(j))
            dd(3)=min(dd(3),y(j))
            dd(4)=max(dd(4),y(j))
        endif
 39   continue
      nplot=j
      dd(4)=3.0
      
      call pgwindow(dd(1),dd(2),
     .  dd(3)-0.1*(dd(4)-dd(3)),dd(4)+0.1*(dd(4)-dd(3)))
      call pgslw(3)
      CALL PGBOX('BCLNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Period(d)","Q8MES/Q3MES","")
c      call pgpt(nplot,x,y,17)

      call pgbbuf()
      do 40 i=1,nplot
        call pgsci(z(i))
        call pgpt1(x(i),y(i),17)
 40   continue
      call pgebuf()
      px(1)=dd(1)
      px(2)=dd(2)
      py(1)=1.73 !ratio of data
      py(2)=1.73
      call pgsci(4)
      call pgline(2,px,py)
      call pgsci(1)     

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call pgpage()

      dd(1)= 9.9e30
      dd(2)=-9.9e30
      dd(3)= 9.9e30
      dd(4)=-9.9e30

      j=0
      do 41 i=1,npt
        if((q6MES(i).gt.0.0).and.(q8MES(i).gt.0.0).and.
     .    (q8Per(i).gt.0.0))then
            j=j+1
            x(j)=log10(q8Per(i))
c            write(0,*) q6Per(i),x(j)
            y(j)=q8MES(i)/q6MES(i)
            if((q6MES(i).gt.7.1).and.(q8MES(i).gt.7.1).and.
     .       (q6MDS(i).gt.s2).and.(q8MDS(i).gt.s2))then
                z(j)=1
                nPmatch=nPmatch+1
            else
                z(j)=2
            endif
            dd(1)=min(dd(1),x(j))
            dd(2)=max(dd(2),x(j))
            dd(3)=min(dd(3),y(j))
            dd(4)=max(dd(4),y(j))
        endif
 41   continue
      nplot=j
      dd(4)=3.0
      
      call pgwindow(dd(1),dd(2),
     .  dd(3)-0.1*(dd(4)-dd(3)),dd(4)+0.1*(dd(4)-dd(3)))
      call pgslw(3)
      CALL PGBOX('BCLNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("Period(d)","Q8MES/Q6MES","")
c      call pgpt(nplot,x,y,17)

      call pgbbuf()
      do 42 i=1,nplot
        call pgsci(z(i))
        call pgpt1(x(i),y(i),17)
 42   continue
      call pgebuf()
      px(1)=dd(1)
      px(2)=dd(2)
      py(1)=1.15 !ratio of data
      py(2)=1.15
      call pgsci(4)
      call pgline(2,px,py)
      call pgsci(1) 
      
      call pgclos()
      
      goto 999
 901  write(0,*) "Usage: tpsdawg (eb,KOI,var).list"
      goto 999
 902  write(0,*) "Cannot open ",filename
      goto 999
 999  end