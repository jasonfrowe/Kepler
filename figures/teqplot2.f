      program teqplot2
C     For plotting koiprops files
      implicit none
      integer nmax,npt,nunit,i,nfp,nunit2,fpflag,j
      parameter(nmax=100000)
      integer fp(nmax)
      real period(nmax),rdr(nmax),rstar(nmax),asemi(nmax),teq(nmax),
     .  koi(nmax),kid,koinum(nmax),period2(nmax),rdr2(nmax),
     .  asemi2(nmax),teq2(nmax),Redrs,Rs,Rl,Ts,Tl,rad(nmax),rad2(nmax),
     .  px(2),py(2),teff(nmax),teff2(nmax),rstar2(nmax)
      character*80 fplist,file1,file2,dumc

      RedRs=109.17 !Size of Sun relative to earth
      Rs=0.4 !smaller radius to plot
      Rl=20.0  !largest radius to plot
      Ts=0.0   !Temperature range to plot
      Tl=2000.0
      
      nunit=10
      nunit2=11
      fplist="koi_FPs.20120918.csv"
      i=1
      open(unit=nunit,file=fplist,status='old',err=903)
 5    read(nunit,*,end=6) koi(i),kid,fp(i)
        if(koi(i)-int(koi(i)).eq.0.0) koi(i)=koi(i)+0.01
        if(fp(i).ge.1)then
c            write(0,*) koi(i),kid,fp(i)
            i=i+1
        endif
      goto 5
 6    continue
      nfp=i-1
      write(0,*) "nFP: ",nfp
      close(nunit)
      
      file1="koiprops_awked.dat"
      file2="koiprops_awked2.dat"
      
      open(unit=nunit,file=file1,status='old',err=901)
      open(unit=nunit2,file=file2,status='old',err=902)

      do 10 i=1,2 !read in header
        read(nunit,*) dumc
        read(nunit2,*) dumc
 10   continue
 
      i=1
 12   read(nunit,*,end=11,err=904) koinum(i),period(i),rdr(i),rstar(i),
     .   asemi(i),teq(i),teff(i)
c         write(6,*) i,koinum(i)
      read(nunit2,*,end=11) koinum(i),period2(i),rdr2(i),rstar2(i),
     .  asemi2(i),teq2(i),teff2(i)
c        if((rstar(i).lt.0.1).or.(rstar2(i).lt.0.1)) goto 12
        rad(i)=rstar(i)*rdr(i)*RedRs
        rad2(i)=rstar2(i)*rdr2(i)*RedRs
        if(koinum(i)-int(koinum(i)).eq.0.0) koinum(i)=koinum(i)+0.01

        fpflag=0
        do 25 j=1,nfp
            if(koinum(i).eq.koi(j))then
                fpflag=1
            endif
 25     continue
c        write(6,*) koinum,fpflag,i
        if(fpflag.eq.0) then
            if((rad(i).lt.20.00).and.(teq(i).gt.185.0)
     .       .and.(teq(i).lt.303.0)) write(6,500) koinum(i),rad(i),
     .       teq(i)
 500            format(F8.2,1X,F5.2,1X,F4.0)
            if((teq(i).lt.303).and.(rad(i).gt.0.2)) i=i+1
c            if(rad(i).gt.0.1) i=i+1
        endif
      goto 12
 11   continue
      npt=i-1
      write(6,*) "number of points: ",npt
      close(nunit)
      close(nunit2)
      
      call pgopen('?')
      call PGPAP ( 7.0 ,1.0)  
c      call PGENV ( -0.5 , 0.5 , -0.5 , 0.5 , 1 , -2 )  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call pgpage()
      call pgwindow(Ts,Tl,log10(Rs),log10(Rl))
      call pgslw(3)
      CALL PGBOX('BCNTS1',0.0,0,'BCLNTS1',0.0,0)
      call pglabel("T\deq\u (K)","Radius (R\d\(2284)\u)","New KIC")
      call pgbbuf()
      call pgsch(1.5)
      
      do 13 i=1,npt
        if(koinum(i).lt.965.0)then
            call pgsci(2)
        elseif(koinum(i).lt.1609)then
            call pgsci(3)
        elseif(koinum(i).lt.2669)then
            call pgsci(4)
        else
            call pgsci(6)
        endif
        if((teq(i).gt.Ts).and.(teq(i).lt.Tl).and.
     .   (log10(rad(i)).gt.log10(Rs)).and.
     .   (log10(rad(i)).lt.log10(Rl)))then
            if(rstar(i).gt.0.9) call pgpt1(teq(i),log10(rad(i)),17)
        endif
 13   continue
      call pgsci(1)

      call pgsch(1.5) !plot the Earth
      call pgpt1(255.0,log10(1.0),12)
c      call pgsch(1.0)
      py(1)=log10(Rs)
      py(2)=log10(Rl)
      px(1)=185.0
      px(2)=185.0
      call pgline(2,px,py)
      px(1)=303.0
      px(2)=303.0
      call pgline(2,px,py)
      px(1)=Ts
      px(2)=Tl
      py(1)=0.0 !Earth sized
      py(2)=0.0
      call pgline(2,px,py)
      py(1)=log10(1.25) !Earth sized
      py(2)=log10(1.25)
      call pgline(2,px,py)
      call pgsch(1.0)     

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call pgpage()
      call pgwindow(Ts,Tl,log10(Rs),log10(Rl))
      call pgslw(3)
      CALL PGBOX('BCNTS1',0.0,0,'BCLNTS1',0.0,0)
      call pglabel("T\deq\u (K)","Radius (R\d\(2284)\u)","Old KIC")
      call pgbbuf()
      call pgsch(1.5)
      
      do 14 i=1,npt
        if(koinum(i).lt.965.0)then
            call pgsci(2)
        elseif(koinum(i).lt.1609)then
            call pgsci(3)
        elseif(koinum(i).lt.2669)then
            call pgsci(4)
        else
            call pgsci(6)
        endif
        if((teq(i).gt.Ts).and.(teq(i).lt.Tl).and.
     .   (log10(rad(i)).gt.log10(Rs)).and.
     .   (log10(rad(i)).lt.log10(Rl)))then
            call pgpt1(teq2(i),log10(rad2(i)),17)
        endif
 14   continue
      call pgsci(1)

      call pgsch(1.5) !plot the Earth
      call pgpt1(255.0,log10(1.0),12)
c      call pgsch(1.0)
      py(1)=log10(Rs)
      py(2)=log10(Rl)
      px(1)=185.0
      px(2)=185.0
      call pgline(2,px,py)
      px(1)=303.0
      px(2)=303.0
      call pgline(2,px,py)
      px(1)=Ts
      px(2)=Tl
      py(1)=0.0 !Earth sized
      py(2)=0.0
      call pgline(2,px,py)
      py(1)=log10(1.25) !Earth sized
      py(2)=log10(1.25)
      call pgline(2,px,py)
      call pgsch(1.0)
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call pgpage()
      call pgwindow(Ts,Tl,log10(Rs),log10(Rl))
      call pgslw(3)
      CALL PGBOX('BCNTS1',0.0,0,'BCLNTS1',0.0,0)
      call pglabel("T\deq\u (K)","Radius (R\d\(2284)\u)","")
      call pgbbuf()
      call pgsch(0.5)
      
      do 15 i=1,npt
        if(teff(i).lt.4500.0)then
            call pgsci(2)
        else
            call pgsci(1)
        endif
        px(1)=teq(i)
        px(2)=teq2(i)
        py(1)=log10(rad(i))
        py(2)=log10(rad2(i))
c        call pgline(2,px,py)

        if((px(1).gt.Ts).and.(px(1).lt.Tl).and.
     .   (px(2).gt.Ts).and.(px(2).lt.Tl).and.
     .   (py(1).gt.log10(Rs)).and.(py(1).lt.log10(Rl)).and.
     .   (py(2).gt.log10(Rs)).and.(py(2).lt.log10(Rl)))then

            if((px(1).eq.px(2)).and.(py(1).eq.py(2)))then
               call pgpt1(teq2(i),log10(rad2(i)),17)
            else
               call pgarro(px(2),py(2),px(1),py(1))
            endif

         endif
c        call pgpt1(teq2(i),log10(rad2(i)),17)
 15   continue

      call pgsch(1.5) !plot the Earth
      call pgpt1(255.0,log10(1.0),12)
c      call pgsch(1.0)
      py(1)=log10(Rs)
      py(2)=log10(Rl)
      px(1)=185.0
      px(2)=185.0
      call pgline(2,px,py)
      px(1)=303.0
      px(2)=303.0
      call pgline(2,px,py)
      px(1)=Ts
      px(2)=Tl
      py(1)=0.0 !Earth sized
      py(2)=0.0
      call pgline(2,px,py)
      py(1)=log10(1.25) !Earth sized
      py(2)=log10(1.25)
      call pgline(2,px,py)
      call pgsch(1.0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call pgpage()
      call pgwindow(150.0,400.0,log10(0.7),log10(2.0))
      call pgslw(3)
      CALL PGBOX('BCNTS1',0.0,0,'BCLNTS1',0.0,0)
      call pglabel("T\deq\u (K)","Radius (R\d\(2284)\u)","")
      call pgbbuf()
      call pgsch(0.5)
      
      do 16 i=1,npt
        if(teff(i).lt.4500.0)then
            call pgsci(2)
        else
            call pgsci(1)
        endif
        px(1)=teq(i)
        px(2)=teq2(i)
        py(1)=log10(rad(i))
        py(2)=log10(rad2(i))
c        call pgline(2,px,py)

        if((px(1).gt.Ts).and.(px(1).lt.Tl).and.
     .   (px(2).gt.Ts).and.(px(2).lt.Tl).and.
     .   (py(1).gt.log10(Rs)).and.(py(1).lt.log10(Rl)).and.
     .   (py(2).gt.log10(Rs)).and.(py(2).lt.log10(Rl)))then

            if((px(1).eq.px(2)).and.(py(1).eq.py(2)))then
               call pgpt1(teq2(i),log10(rad2(i)),17)
            else
               call pgarro(px(2),py(2),px(1),py(1))
            endif

         endif

 16   continue
      call pgsci(1)

      call pgsch(1.5) !plot the Earth
      call pgpt1(255.0,log10(1.0),12)
c      call pgsch(1.0)
      py(1)=log10(Rs)
      py(2)=log10(Rl)
      px(1)=185.0
      px(2)=185.0
      call pgline(2,px,py)
      px(1)=303.0
      px(2)=303.0
      call pgline(2,px,py)
      px(1)=Ts
      px(2)=Tl
      py(1)=0.0 !Earth sized
      py(2)=0.0
      call pgline(2,px,py)
      py(1)=log10(1.25) !Earth sized
      py(2)=log10(1.25)
      call pgline(2,px,py)
      call pgsch(1.0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      call pgpage()
      
      call PGHIST(npt,rad,0.0,15.0,50,0)
      call pglabel("Radius (R\d\(2284)\u)","Frequency (N)","New KIC")
      call PGHIST(npt,rad2,0.0,15.0,50,0)
      call pglabel("Radius (R\d\(2284)\u)","Frequency (N)","Old KIC")
      call PGHIST(npt,teq,0.0,2000.0,50,0)
      call pglabel("T\deq\u (K)","Frequency (N)","New KIC")
      call PGHIST(npt,teq2,0.0,2000.0,50,0)
      call pglabel("T\deq\u (K)","Frequency (N)","Old KIC")
      call PGHIST(npt,asemi,0.0,1.0,50,0)
      call pglabel("a (AU)","Frequency (N)","New KIC")
      call PGHIST(npt,asemi2,0.0,1.0,50,0)
      call pglabel("a (AU)","Frequency (N)","New KIC")
      
      call pgclos()

      goto 999
 901  write(0,*) "Cannot open ",file1
      goto 999
 902  write(0,*) "Cannot open ",file2      
      goto 999
 903  write(0,*) "Cannot open ",fplist
      goto 999
 904  write(0,*) "Error: ",i,koinum(i-1)
      goto 999
 999  end 
