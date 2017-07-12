      program koi70p4plot
      implicit none
      integer nmax,nfit,npt,iargc,nunit,nplanet,i,nplanet2,nplanetmax,j,
     .  col,k,order(10),ii,nplot,iii,npt2,bins,kk,kk2
      parameter (nmax=600000,nfit=108,nplanetmax=10,nplot=1000)
      integer dtype(nmax),dtype2(nmax),flag(nmax)
      real px(nmax),py(nmax),rbb(4),pz(nmax)
      double precision sol(nfit),time(nmax),tmodel(nmax),
     .  flux(nmax),ferr(nmax),itime(nmax),Keplertime,serr(nfit,2),
     .  err(nfit),sol2(nfit),tdur(nplanetmax),transitdur,tdurmax,
     .  transitdepth,tdepth(nplanetmax),tdepthmax,toff,ph1,ph2,
     .  phase(nmax),period,mmax,mmin,width,time2(nmax),itime2(nmax),t1,
     .  t2,ratio,bx(2),mintime,flux2(nmax),
     .  ferr2(nmax),phase2(nmax),std,mean,stdev,phd1,phd2,sn
      character*80 inputsol,obsfile,label
c      data order /5,4,1,3,2,6,7,8,9,10/
c      data order /1,2,3,4,5,6,7,8,9,10/
      data order /2,1,3,4,5,6,7,8,9,10/
      
      if(iargc().lt.2) goto 901
      
C     Parse the name of the observations data file from the commandline
      call getarg(1,obsfile)
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,flux,ferr,itime,Keplertime)
      call readkeplc(nunit,nmax,npt,time,flux,ferr,itime,Keplertime)
      close(nunit)!release unit number as we are done with the file.
      
      mmin=flux(1)
      mmax=flux(1)
      mintime=time(1)
      ferr(1)=1.0e-5
      do 15 i=2,npt
        mmin=min(mmin,flux(i))
        mmax=max(mmax,flux(i))
        mintime=min(mintime,time(i))
        ferr(i)=1.0e-5
 15   continue
      
      call getarg(2,inputsol) !get filename for input solution
      nunit=10 !unit number used for file input
      open(unit=nunit,file=inputsol,status='old',err=902)
C     We start by reading in solution from input file
      write(0,*) "reading in input solution"
      call getfitpars(nunit,nfit,nplanet,sol,serr,err)
      write(0,*) "done reading input solution"
c      call getfitpars(nunit,nfit,sol,serr,err)
      close(nunit) !release unit number as we are done with file
      
      do 17 i=1,npt !central database of all data
        dtype(i)=0 !0 marks that we have photometric data
 17   continue

      call pgopen('?') !open PGPlot device
      call pgask(.true.) !don't ask for new page.. just do it.
c      call PGPAP (10.0/real(nplanet+1) ,real(nplanet+1)/2.0) !paper size
      call PGPAP(8.0/real(nplanet+1)*2.0 ,real(nplanet+1)/2.0) 
      call pgsubp(1,nplanet+1)  !break up plot into grid
 
      do 10 i=1,nplanet
        tdur(i)=transitdur(i,nfit,sol)/8.64d4
        tdepth(i)=transitdepth(i,nfit,sol,sol2)
c        write(6,*) "tdur:",tdur(i),tdepth(i)
 10   continue
      
      tdurmax=tdur(1)
      tdepthmax=tdepth(1)
      do 11 i=2,nplanet
        tdurmax=max(tdur(i),tdurmax)
        tdepthmax=max(tdepth(i),tdepthmax)
 11   continue     

      rbb(3)=real(mmin)
c      rbb(3)=0.9997
      rbb(4)=real(mmax)
c      rbb(4)=1.00025

      do 12 ii=1,nplanet
        i=order(ii)
        col=10*(i-1)

c       ratio=(sol(10+10*(4-1))/sol(10+10*(order(ii)-1)))**(1.0d0/3.0d0)
        ratio=1.0

c        if(ii.lt.nplanet)then
c            ratio=sol(10+10*(order(ii+1)-1))/sol(10+10*(order(ii)-1))
c            ratio=ratio**(1.0d0/3.0d0)
c        else
c            ratio=1.0
c        endif
      
        do 13 j=1,nfit
            sol2(j)=sol(j)
 13     continue
        sol2(12+col)=0.0d0 !set r/R*=0
        period=sol2(10+col)
        call transitmodel(nfit,nplanet,sol2,npt,time,itime,tmodel,dtype)
        
        toff=0.5-(sol(9+col)/sol(10+col)-int(sol(9+col)/sol(10+col)))
        if(toff.lt.0.0)toff=toff+1.0
        ph1=0.5-1.25d0*tdurmax/sol(10+col)
        if(ph1.lt.0.25)ph1=0.25
        ph2=0.5+1.25d0*tdurmax/sol(10+col)
        if(ph2.gt.0.75)ph2=0.75
        
        rbb(1)=real(ph1)
        rbb(2)=real(ph2)
        
        call phasept(npt,time,phase,period,toff,ratio)
        
        k=0
        kk=0
        kk2=0
        phd1=0.5-0.5*tdur(i)/sol(10+col)
        phd2=0.5+0.5*tdur(i)/sol(10+col)
        do 14 j=1,npt
            if((phase(j).gt.ph1).and.(phase(j).lt.ph2))then
                    k=k+1
                    px(k)=real(phase(j))
                    py(k)=real(flux(j)-tmodel(j))+1.0
                    
c                write(0,*) px(k),py(k)
            endif
            if((phase(j).lt.phd1).or.(phase(j).gt.phd2))then
                    kk=kk+1
                    flux2(kk)=flux(j)-tmodel(j)
            else
                kk2=kk2+1
            endif
c            write(0,*) kk,kk2
c            read(5,*)
 14     continue
        std=stdev(kk,flux2,mean)
        sn=tdepth(i)/std*sqrt(dble(kk2))
        write(0,*) "S/N:",i,tdepth(i)/std*sqrt(dble(kk2))
c        write(0,*) tdepth(i),std,kk2
        call pgpage()
        call pgvport(0.25,0.95,0.0,1.0)
c        if(ii.le.2) goto 12
        call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))
c        write(6,*) rbb(1),rbb(2),rbb(3),rbb(4)
        call pgsch(1.0)
        call pgpt(k,px,py,17)
        call pgsch(3.5)
        
        write(label,500) 'S/N=',sn
 500    format(A4,F5.1)
c        call pgptxt(rbb(2)-0.20*(rbb(2)-rbb(1)),rbb(3)
c     .      +0.10*(rbb(4)-rbb(3)),0.0,0.5,label)

CCCCCCCCCCCCCCCCCCCCCCCCCC       

        k=0
        do 6 j=1,npt
                if((flux(j).gt.rbb(3)).and.(flux(j).lt.rbb(4)))then
                    k=k+1
                    phase2(k)=phase(j)
                    flux2(k)=flux(j)-tmodel(j)+1.0d0
                    ferr2(k)=ferr(j)
                endif
 6      continue
        npt2=k
        
        bins=int(sol(10*(i-1)+8+2)*1440.0d0/30.0d0)
        call binp(npt2,phase2,flux2,ferr2,bins,flag)
        j=0
        do 21 iii=1,npt2
            if(flag(i).eq.0)then
                if((phase2(iii).gt.ph1).and.(phase2(iii).lt.ph2))then
                    j=j+1
                    px(j)=real(phase2(iii))
                    py(j)=real(flux2(iii))
                    pz(j)=real(ferr2(iii))
                endif
            endif
 21     continue
        npt2=j
c      write(0,*) "npt2:",npt2
        call pgsci(4)
c        call pgsch(1.5)
        call pgpt(npt2,px,py,17)
c        call pgsch(0.8)
        call pgerrb(6,npt2,px,py,pz,1.0)
        call pgsci(1)
CCCCCCCCCCCCCCCCC

c        do 16 j=1,nplanet
c            do 18 k=1,nfit
c                sol2(k)=sol(k)
c 18         continue
c            if(i.ne.j) sol2(12+10*(j-1))=0.0d0 !only want one transit
c            write(6,*) 'rdr',j,sol2(12+10*(j-1))
c 16     continue

        do 16 j=1,8
            sol2(j)=sol(j)
  16    continue
        do 18 j=9,18          
            sol2(j)=sol(j+col)
  18    continue
c        sol2(8)=0.0
        t1=sol2(9+col)-1.5*tdurmax
        t2=sol2(9+col)+1.5*tdurmax
        do 19 k=1,nplot
            time2(k)=t1+(t2-t1)/dble(nplot)*dble(k-1)
            itime2(k)=itime(1)
            dtype2(k)=0
 19     continue
        call transitmodel(nfit,1,sol2,nplot,time2,itime2,tmodel,dtype2)
        call phasept(nplot,time2,phase,period,toff,ratio)
        do 20 k=1,nplot
            px(k)=real(phase(k))
            py(k)=real(tmodel(k))
c            write(0,*) px(k),py(k)
 20     continue
c        read(5,*)
        call sort2(k,px,py)
c        call pgsci(i+1)
        call pgsci(2)
        call pgline(nplot,px,py)
        
        bx(1)=sol2(9)-tdur(i)/2.0d0
        bx(2)=sol2(9)+tdur(i)/2.0d0
c        write(0,*) sol2(9),tdur(i)
        call phasept(2,bx,phase,period,toff,ratio)
c        write(6,*) "bx:",phase(1)
        px(1)=real(phase(1))
        px(2)=px(1)
        py(1)=rbb(3)
        py(2)=rbb(4)
c        call pgline(2,px,py)
        px(1)=real(phase(2))
        px(2)=px(1)
c        call pgline(2,px,py)
        
        call pgsci(1)

        call pgsch(3.5)
        
        if(ratio.eq.1.0)then
            if(ii.eq.nplanet) call pgptxt((rbb(1)+rbb(2))/2.0,rbb(3)-
     .          0.20*(rbb(4)-rbb(3)),0.0,
     .          0.5,"Time from mid-transit (hours)")
        else     
           if(ii.eq.nplanet) call pgptxt((rbb(1)+rbb(2))/2.0,rbb(3)-
     .          0.20*(rbb(4)-rbb(3)),0.0,
     .          0.5,"Scaled Transit Duration")
        endif
     
        call pgptxt(rbb(1)-0.30*(rbb(2)-rbb(1)),(rbb(3)+rbb(4))/2,
     .      90.0,0.5,"Relative Flux")
     
        if(ratio.eq.1.0)then
            width=(ph2-ph1)*period*24.0d0/2.0d0
        else
            width=1.0
        endif
        call pgwindow(real(-width),real(width),rbb(3),rbb(4))        
        if(ii.eq.nplanet)then
            call pgbox('BCNTS1',0.0,0,'BCNTSV1',0.0,0)
        else
            call pgbox('BCTS1',0.0,0,'BCNTSV1',0.0,0)
        endif
        call pgwindow(real(-width),real(width),rbb(3),rbb(4)) 
        call pgsch(1.0)
        
c        read(5,*)
        
 12   continue
      
      call pgclos()
      
      goto 999
 901  write(0,*) "Usage: transitplot5 <photfile> <fitpars>"
      goto 999
 902  write(0,*) "Error opening ",inputsol
      goto 999
 903  write(0,*) "Error opening ",obsfile
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine binp(npt,phase,mag,merr,bins,flag)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,bins,N,i,j
      parameter(N=600000)
      integer flag(n)
      double precision phase(N),tt(N),mag(N),tn(N),merr(N),te(n),
     .  tp(N)
      
c      write(0,*) "Bins:",bins

      do 5 i=1,bins
         tn(i)=0.0
         tp(i)=0.0
         tt(i)=0.0
         te(i)=0.0
         flag(i)=0
 5    continue

      do 10 i=1,npt
         j=int(real(bins)*phase(i))+1
         if(j.le.bins)then
            tn(j)=tn(j)+1.0/merr(i)
            tp(j)=tp(j)+phase(i)/merr(i)
            tt(j)=tt(j)+mag(i)/merr(i)
            te(j)=te(j)+1.0
         endif
 10   continue

      do 20 i=1,bins
c         phase(i)=(dble(i)-0.5)/dble(bins)
         phase(i)=tp(i)/tn(i)
         mag(i)=tt(i)/tn(i)
         merr(i)=(te(i)**0.5)/tn(i)
         if(tn(i).eq.0.0) flag(i)=1
c         write(0,*) flag(i)
c         read(5,*)
c         write(6,*) phase(i),mag(i),merr(i)
 20   continue
      npt=bins

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine sort2(n,x,y)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer n,i,nmax
      parameter(nmax=600000)
      integer p(nmax)
      real x(n),y(n),tx(nmax),ty(nmax)

      call rqsort(n,x,p)

      do 10 i=1,n
         tx(i)=x(p(i))
         ty(i)=y(p(i))
 10   continue
      
      do 15 i=1,n
         x(i)=tx(i)
         y(i)=ty(i)
 15   continue

      return
      end

c**********************************************************************
      subroutine rqsort(n,a,p)
c======================================================================
c     Return integer array p which indexes array a in increasing order.
c     Array a is not disturbed.  The Quicksort algorithm is used.
c
c     B. G. Knapp, 86/12/23
c
c     Reference: N. Wirth, Algorithms and Data Structures,
c     Prentice-Hall, 1986
c======================================================================
      implicit none

c     Input:
      integer   n
      real      a(n)

c     Output:
      integer   p(n)

c     Constants
      integer   LGN, Q
      parameter (LGN=32, Q=11)
c        (LGN = log base 2 of maximum n;
c         Q = smallest subfile to use quicksort on)

c     Local:
      real      x
      integer   stackl(LGN),stackr(LGN),s,t,l,m,r,i,j

c     Initialize the stack
      stackl(1)=1
      stackr(1)=n
      s=1

c     Initialize the pointer array
      do 1 i=1,n
         p(i)=i
    1 continue

    2 if (s.gt.0) then
         l=stackl(s)
         r=stackr(s)
         s=s-1

    3    if ((r-l).lt.Q) then

c           Use straight insertion
            do 6 i=l+1,r
               t = p(i)
               x = a(t)
               do 4 j=i-1,l,-1
                  if (a(p(j)).le.x) goto 5
                  p(j+1) = p(j)
    4          continue
               j=l-1
    5          p(j+1) = t
    6       continue
         else

c           Use quicksort, with pivot as median of a(l), a(m), a(r)
            m=(l+r)/2
            t=p(m)
            if (a(t).lt.a(p(l))) then
               p(m)=p(l)
               p(l)=t
               t=p(m)
            endif
            if (a(t).gt.a(p(r))) then
               p(m)=p(r)
               p(r)=t
               t=p(m)
               if (a(t).lt.a(p(l))) then
                  p(m)=p(l)
                  p(l)=t
                  t=p(m)
               endif
            endif

c           Partition
            x=a(t)
            i=l+1
            j=r-1
    7       if (i.le.j) then
    8          if (a(p(i)).lt.x) then
                  i=i+1
                  goto 8
               endif
    9          if (x.lt.a(p(j))) then
                  j=j-1
                  goto 9
               endif
               if (i.le.j) then
                  t=p(i)
                  p(i)=p(j)
                  p(j)=t
                  i=i+1
                  j=j-1
               endif
               goto 7
            endif

c           Stack the larger subfile
            s=s+1
            if ((j-l).gt.(r-i)) then
               stackl(s)=l
               stackr(s)=j
               l=i
            else
               stackl(s)=i
               stackr(s)=r
               r=j
            endif
            goto 3
         endif
         goto 2
      endif
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function stdev(npt,pts,mean)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Calculates standard deviation of data set given the mean.
      implicit none

      integer npt,i
      double precision pts(npt),mean,s,ep,adev,p,var,sdev,mean2

      s=0.
      do 11 i=1,npt
         s=s+pts(i)
 11   continue
      mean=s/npt

      ep=0.
      var=0.
      do 10 i=1,npt
         s=pts(i)-mean
         ep=ep+s
         p=s*s
         var=var+p
 10   continue
      var=(var-ep**2/npt)/(npt-1)
      sdev=sqrt(var)

      stdev=sdev

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine phasept(npt,time,phase,period,toff,ratio)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     npt - number of data points
C     time - times of observations
C     mag - the observations
C     phase - returned phase of observation
C     period - fixed period for data
      implicit none
      integer npt
      double precision time(npt),phase(npt),period,toff,ratio

      integer i
      double precision temp

      do 10 i=1,npt
         temp=time(i)
C        Get the phase
         phase(i)=temp/period-int(temp/period)
C        apply optional phase offset to make plot pretty
         phase(i)=phase(i)+toff
C        make sure phase is between 0 and 1
         if(phase(i).lt.0.0) phase(i)=phase(i)+1.0
         if(phase(i).gt.1.0) phase(i)=phase(i)-1.0
         phase(i)=(phase(i)-0.5)*ratio+0.5
c         write(6,*) ratio,phase(i)
 10   continue
c      read(5,*)
      return
      end
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function transitdepth(np,nfit,sol,sol2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer np,nfit,col,i
      double precision sol(nfit),itime(1),dtype(1),sol2(nfit),tmodel(1),
     .  time(1)
      
      itime(1)=1765.5/86400.0d0
      dtype(1)=0
      
      col=10*(np-1)
      
      do 10 i=1,8
        sol2(i)=sol(i)
 10   continue
      do 11 i=9,18
        sol2(i)=sol(i+col)
c        write(6,*) i,sol2(i)
 11   continue       
      
      time(1)=sol2(9)
      
      call transitmodel(nfit,1,sol2,1,time,itime,tmodel,dtype)
      transitdepth=(1.0d0-tmodel(1)-sol(8))
      
      return
      end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function transitdur(np,nfit,sol)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer np,nfit
      double precision sol(nfit),b,Psec,G,aConst,Pi,adrs,cincl,temp(4),
     .  bb,rdr
      
      Pi=acos(-1.d0)   !Pi
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)
      
      bb=sol(11+10*(np-1))
      b=sqrt(bb)
      Psec=sol(10+10*(np-1))*8.64d4 !sec ; period of planet
      adrs=1000.0*sol(1)*G*(Psec)**2/(3.0d0*Pi)
      adrs=adrs**(1.0d0/3.0d0) !a/R*
      cincl=b/adrs !cos(i)
      rdr=sol(12+10*(np-1))
        
      temp(1)=Psec/Pi
      temp(2)=1.0d0/adrs
      temp(3)=(1+rdr)**2.0-bb
      temp(4)=1-cincl*cincl
C     Transit duration in days
      transitdur=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readkeplc(nunit,nmax,npt,dtime,flux,ferr,itime,
     .  Keplertime)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,nmax,npt,i
      double precision dtime(nmax),flux(nmax),ferr(nmax),itime(nmax),
     .  Keplertime,sec2day,mintime
            
C     if time=0, then time=MOSTtime (that does not mean drink up)
      Keplertime=54900.0
c      Keplertime=0.5d0

C     number of seconds in a day
      sec2day=86400.0d0
      
      i=1
      
      mintime=99.9d30
  9   continue
 10   read(nunit,*,err=9,end=20) dtime(i),flux(i),ferr(i)
cc  reset times form arbitrary start of 0.0 to HJD of field center
cc  add 53.038152 to MJD center 1st exp
c        dtime(i)=dtime(i)+53.038152
cc  add HJD 5/2/09 adjustment for field center
c        dtime(i)=dtime(i)+0.50020-Keplertime
c  add time dependent change for field center
c        dtime(i)=dtime(i)+4.1d-5*(dtime(i)-53.0)
        dtime(i)=dtime(i)-Keplertime
c        dtime(i)=dtime(i)+sin(2.0d0*pi*dtime(i)/372.5d0-
c     .      1.1208d0)*0.00280758d0
        dtime(i)=dtime(i)+0.5d0 !MJD half day offset
        mintime=min(mintime,dtime(i))
        flux(i)=flux(i)+1.0!-2.5*log10(mag(i)+1.0d0)
c        ferr(i)=0.00005
        itime(i)=1765.5/sec2day
        i=i+1
      goto 10
 20   continue
        
      npt=i-1
      write(0,*)   "-------------------------"
      write(0,500) "Observations read: ",npt
      write(0,*) "Mintime: ",mintime
      write(0,*)   "-------------------------"
 500  format(1X,A19,I6)
 
c      Keplertime=Keplertime+mintime  !correct time=0 time
c      do 30 i=1,npt
c         dtime(i)=dtime(i)-mintime
c 30   continue
 
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getfitpars(nunit,nfit,nplanet,sol,serr,err)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,nunit,i,nplanet,np,j
      double precision sol(nfit),serr(nfit,2),err(nfit)
      character*160 command
      
c      write(0,*) "----------------------------"
c      write(0,*) "Reading in starting solution"
c      write(0,*) "----------------------------"
      
      nplanet=0 !initialize number of planets
      
      i=0 !initialize line counter
C     Start of loop to read in file
  10  read(nunit,500,end=11,err=901) command
 500  format(A160) !command much be contained within first 80 characters
        i=i+1 !increase line counter
        if(command(1:1).eq."#") then
c            write(0,*) "Comment line on line ",i
        elseif(command(1:3).eq."RHO") then
            read(command(5:160),*) sol(1),serr(1,1),serr(1,2),err(1)
c            write(0,*) "WTF!",sol(1),serr(1,1)
        elseif(command(1:3).eq."NL1") then
            read(command(5:160),*) sol(2),serr(2,1),serr(2,2),err(2)
        elseif(command(1:3).eq."NL2") then
            read(command(5:160),*) sol(3),serr(3,1),serr(3,2),err(3)
        elseif(command(1:3).eq."NL3") then
            read(command(5:160),*) sol(4),serr(4,1),serr(4,2),err(4)
        elseif(command(1:3).eq."NL4") then
            read(command(5:160),*) sol(5),serr(5,1),serr(5,2),err(5)
        elseif(command(1:3).eq."DIL") then
            read(command(5:160),*) sol(6),serr(6,1),serr(6,2),err(6)
        elseif(command(1:3).eq."VOF") then
            read(command(5:160),*) sol(7),serr(7,1),serr(7,2),err(7)
        elseif(command(1:3).eq."ZPT") then
            read(command(5:160),*) sol(8),serr(8,1),serr(8,2),err(8)
        elseif(command(1:2).eq."EP") then
            read(command(3:3),*,err=901) np
            if(np.gt.nplanet)nplanet=np
C           j=planet parameters*(np-1)+8 initial parameters
            j=10*(np-1)+8+1
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."PE") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+2
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."BB") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+3
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."RD") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+4
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."EC") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+5
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."ES") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+6
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."KR") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+7
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."TE") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+8
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."EL") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+9
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        elseif(command(1:2).eq."AL") then
            read(command(3:3),*,err=901) np
            j=10*(np-1)+8+10
            read(command(5:160),*) sol(j),serr(j,1),serr(j,2),err(j)
        endif
 501    format(A5,5(1X,1PE17.10))
c        write(0,*) command
      goto 10 !loop back to read statement
 11   continue !end loop when EOF is reached
c      write(0,*) "----------------------------"

      goto 999
 901  write(0,*) "Error on line ",i+1
      pause
      goto 999       
 999  return
      end 
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine transitmodel(nfit,nplanet,sol,npt,time,itime,tmodel,
     .  dtype)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,npt,i,j,nintg,dtype(npt),ii,nplanet
      parameter(nintg=11)
      double precision sol(nfit),per,epoch,b,RpRs,tmodel(npt),
     .  time(npt),phi,adrs,bt(nintg),bs2,Pi,tpi,c1,c2,c3,c4,
     .  itime(npt),t,tflux(nintg),dnintg,dnintgm1,pid2,zpt,eccn,w,Eanom,
     .  Tanom,trueanomaly,phi0,vt(nintg),K,voff,drs,distance,dil,G,
     .  esw,ecw,Manom,Ted,ell,Ag,Cs,fDB,tide(nintg),alb(nintg),
     .  albedomod,phase,ratio,ab,tm,mu(nintg)
      
      Pi=acos(-1.d0)!define Pi and 2*Pi
      tPi=2.0d0*Pi 
      pid2=Pi/2.0d0
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      Cs=2.99792458e8 !Speed of light
      fDB=0.0!1.896 !Doppler Boosting factor
      
      c1=sol(2)      !non-linear limb-darkening
      c2=sol(3)
      c3=sol(4)
      c4=sol(5)
      dil=sol(6)     !dilution parameter (model scaling)
      voff=sol(7)    !velocity zero point
      zpt=sol(8)     !flux zero point.
      
      do 17 i=1,npt
        tmodel(i)=0.0d0
 17   continue
      
      do 16 ii=1,nplanet
      
        per=sol(10*(ii-1)+8+2)     !Period (days)
        bs2=abs(sol(10*(ii-1)+8+3))
        b=sqrt(bs2)       !impact parameter
        RpRs=sol(10*(ii-1)+8+4)    !Rp/R*
        ecw=sol(10*(ii-1)+8+5)
        esw=sol(10*(ii-1)+8+6)
        eccn=sqrt(ecw*ecw+esw*esw) !eccentricity
        if(eccn.ge.1.0) eccn=0.99
        if(eccn.eq.0.0d0)then
            w=0.0d0
        else
            w=atan(esw/ecw)
            if((ecw.gt.0.0d0).and.(esw.lt.0.0d0))then
                w=tPi+w
            elseif((ecw.lt.0.0d0).and.(esw.gt.0.0d0))then 
                w=Pi+w
            elseif((ecw.lt.0.0d0).and.(esw.lt.0.0d0))then
                w=Pi+w
            endif
        endif
c        write(0,*) sol(7),sol(8),w
c        write(0,*) "w:",acos(sol(7)/eccn),asin(sol(8)/eccn)
c        read(5,*)

C       a/R*
c        adrs=sol(5)*per/tpi*sqrt(1-sol(3))*(1+sol(8))/sqrt(1-eccn*eccn)
        adrs=1000.0*sol(1)*G*(Per*86400.0d0)**2/(3.0d0*Pi)
        adrs=adrs**(1.0d0/3.0d0)
c        write(0,*) "a/R*:",adrs

        K=sol(10*(ii-1)+8+7)
      
        ted=sol(10*(ii-1)+8+8)/1.0d6 !Occultation Depth
        ell=sol(10*(ii-1)+8+9)/1.0d6 !Ellipsoidal variations
        ag=sol(10*(ii-1)+8+10)/1.0d6 !Phase changes
      
        dnintg=dble(nintg) !convert integer to double
        dnintgm1=2.0*dnintg-2.0
      
C     Find phase at centre of transit
        Eanom=w
        epoch=sol(10*(ii-1)+8+1)   !center of transit time (days)
        Manom=w !mean anomaly
        call kepler(Manom,Eanom,eccn) !eccentric anomaly
        phi0=trueanomaly(eccn,Eanom)
      
        do 10 i=1,npt
            do 11 j=1,nintg
                tflux(j)=0.0 !initialize model
C               sample over integration time
                t=time(i)+itime(i)*(2.0*dble(j)-dnintg-1.0)/dnintgm1-
     .              epoch
C               get orbital position (mean anomaly)
                phi=t/per-floor(t/per)
                phi=phi*tPi
                Manom=phi+w
                if(Manom.gt.tPi) Manom=Manom-tPi
                if(Manom.lt.0.0d0) Manom=Manom+tPi
                call kepler(Manom,Eanom,eccn)
                Tanom=trueanomaly(eccn,Eanom)
                if(phi.gt.Pi) phi=phi-tPi            
                drs=distance(adrs,eccn,Tanom)
                bt(j)=sqrt(bs2+(drs*sin(Tanom-phi0))**2)
C               Correct for light-travel time!
c            if((abs(bt(j))-RpRs.le.1.0d0).and.(abs(phi).gt.Pid2))then
c              t=time(i)-ltt+itime(i)*(2.0*dble(j)-dnintg-1.0)/dnintgm1
c     .              -epoch
c              phi=(t-ltt)/per-floor((t-ltt)/per)
c              phi=phi*tPi
c              Manom=phi+w
c              if(Manom.gt.tPi) Manom=Manom-tPi
c              if(Manom.lt.0.0d0) Manom=Manom+tPi
c              call kepler(Manom,Eanom,eccn)
c              Tanom=trueanomaly(eccn,Eanom)
c                if(phi.gt.Pi) phi=phi-tPi            
c                drs=distance(adrs,eccn,Tanom)
c                bt(j)=sqrt(bs2+(drs*sin(Tanom-phi0))**2)
c            endif
                vt(j)=K*(cos(Pid2+Tanom-phi0)+eccn*cos(w))
                tide(j)=ell*cos(2.0d0*(Pid2+Tanom-phi0))
                alb(j)=albedomod(Pi,ag,Tanom-phi0)
            
                if(j.eq.nintg/2+1)then
                    phase=Tanom-phi0!phi(nintg/2+1)
                    if(phase.gt.Pi) phase=phase-tPi
                    if(phase.lt.-Pi) phase=phase+tPi
                endif

 11         continue
            if(dtype(i).eq.0)then
                if(abs(phase).lt.Pid2)then
C       If we have a transit
                    if((c3.eq.0.0).and.(c4.eq.0.0))then
                       call occultquad(bt,c1,c2,RpRs,tflux,mu,nintg)
                    else
                       call occultsmall(RpRs,c1,c2,c3,c4,nintg,bt,tflux)
                    endif
                    tm=0.0d0
                    do 12 j=1,nintg
                        if(RpRs.le.0.0)tflux(j)=1.0d0
C                   model=transit+doppler+ellipsodial 
                        tm=tm+tflux(j)-fDB*vt(j)/Cs+tide(j)+alb(j)
 12                 continue
                    tm=tm/dnintg
                else
C       We have an eclipse
                    tm=0.0d0
                    do 14 j=1,nintg
                        ratio=1.0d0
                        ab=dabs(bt(j))
                        if((ab.ge.1.0d0).and.(ab-RpRs.le.1.0d0))then
                            ratio=(1.0d0+RpRs-ab)/(2.0d0*RpRs)
                        elseif((ab.lt.1.0d0).and.(ab+RpRs.ge.1.0d0))then
                            ratio=(RpRs+1.0d0-ab)/(2.0d0*RpRs)
                        elseif(ab-RpRs.gt.1.0d0)then
                            ratio=0.0d0
                        endif
                        if(RpRs.le.0.0d0) ratio=0.0d0
c                    write(0,*) ab,RpRs,ratio
c                    read(5,*) 
                        tm=tm+(1.0d0-ted*ratio)
     .                      -fDB*vt(j)/Cs+tide(j)+alb(j)
 14                 continue
                    tm=tm/dnintg
                endif
                tm=tm+(1.0d0-tm)*dil-1.0d0!add dilution
            else
                tm=0.0d0
                do 13 j=1,nintg
                    tm=tm+vt(j)
 13             continue
                tm=tm/dnintg
c            write(0,*) "rv:",tmodel(i)
c            read(5,*)
            endif
            tmodel(i)=tmodel(i)+tm
 10     continue
 
C     Need to add zero points (voff, zpt)
      
c        do 9 i=1,npt
c            write(6,*) time(i),tmodel(i)
c 9      continue
  
 16   continue  
 
      do 15 i=1,npt
        if(dtype(i).eq.0)then
            tmodel(i)=tmodel(i)+zpt+1.0d0
        else
            tmodel(i)=tmodel(i)+voff
        endif
 15   continue
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function albedomod(Pi,ag,phi)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision Pi,phi,alpha,phase,ag

      phi=phi+Pi
      if(phi.gt.2.0*Pi) phi=phi-2.0*Pi


      alpha=abs(phi)      
c      alpha=2.0*Pi*t/Per+phi
      alpha=alpha-2.0*Pi*int(alpha/(2.0*Pi))
      if(alpha.gt.Pi) alpha=abs(alpha-2.0*pi)
c      write(6,*) t,alpha
c      phase=(1.0d0+cos(alpha))/2.0d0
      phase=(sin(alpha)+(Pi-alpha)*cos(alpha))/Pi  !Lambertian Sphere
      
      albedomod=ag*phase
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function mandelagol(nintg,R1,R2,x1,x2,y1,y2,c,
     .  b0,mu,mulimb0,mulimbf,dist)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Adapted from Mandel and Agol, 2002, ApJ 580, L171
      implicit none
      integer i,nintg,sflag
      double precision R1,R2,x1,x2(nintg),y1,y2(nintg),c(4),
     .  c1,c2,c3,c4,mulimb0(nintg),mulimbf(nintg,5),rl,b0(nintg),
     .  mu(nintg),dist(nintg)
      
      mu=0

      c1=c(1)
      c2=c(2)
      c3=c(3)
      c4=c(4)
      rl=R2/R1
      sflag=0
      do 10 i=1,nintg
       dist(i)=Sqrt((x2(i)-x1)*(x2(i)-x1)+(y2(i)-y1)*(y2(i)-y1))/(R1+R2)
c        if(dist(i).ge.1.0d0)then
c            sflag=sflag+1
c            b0(i)=2.0
c        else
            b0(i)=(R1+R2)*dist(i)/R1
c        endif
 10   continue
      
c      write(6,500) "hello",(b0(i),i=1,nintg)
 500  format(A5,11(1X,F5.3))
      if(sflag.lt.nintg) then
c          call occultnl(rl,c1,c2,c3,c4,b0,mulimb0,mulimbf,nintg)
          call occultsmall(rl,c1,c2,c3,c4,nintg,b0,mulimb0)
      endif
c      if(b0(1).le.1.0) write(6,*)b0(1),mu

C      mandelagol=mulimb0(1)
      mandelagol=0.0
      do 11 i=1,nintg
c        mandelagol=mandelagol+mu(i)
c        if(dist(i).ge.1.0d0) mulimb0(i)=1.0d0
        mandelagol=mandelagol+mulimb0(i)
 11   continue
      mandelagol=mandelagol/dble(nintg)
c      write(6,*) "hello2",mandelagol
   
      return
      end

      subroutine occultsmall(p,c1,c2,c3,c4,nz,z,mu)
      implicit none
      integer i,nz
c      parameter (nz=201)
      real*8 p,c1,c2,c3,c4,z(nz),mu(nz),i1,norm,
     &       x,tmp,iofr,pi
C This routine approximates the lightcurve for a small 
C planet. (See section 5 of Mandel & Agol (2002) for
C details):
C Input:
C  p      ratio of planet radius to stellar radius
C  c1-c4  non-linear limb-darkening coefficients
C  z      impact parameters (positive number normalized to stellar 
C        radius)- this is an array which MUST be input to the routine
C  NOTE:  nz must match the size of z & mu in calling routine
C Output:
C  mu     flux relative to unobscured source for each z
C
      pi=acos(-1.d0)
      norm=pi*(1.d0-c1/5.d0-c2/3.d0-3.d0*c3/7.d0-c4/2.d0)
      i1=1.d0-c1-c2-c3-c4
      do i=1,nz
        mu(i)=1.d0
        if(z(i).gt.1.d0-p.and.z(i).lt.1.d0+p) then
          x=1.d0-(z(i)-p)**2
          tmp=(1.d0-c1*(1.d0-0.8d0*x**0.25)
     &             -c2*(1.d0-2.d0/3.d0*x**0.5)
     &             -c3*(1.d0-4.d0/7.d0*x**0.75)
     &             -c4*(1.d0-0.5d0*x))
          mu(i)=1.d0-tmp*(p**2*acos((z(i)-1.d0)/p)
     &        -(z(i)-1.d0)*sqrt(p**2-(z(i)-1.d0)**2))/norm
        endif
        if(z(i).le.1.d0-p.and.z(i).ne.0.d0) then
          mu(i)=1.d0-pi*p**2*iofr(c1,c2,c3,c4,z(i),p)/norm
        endif
        if(z(i).eq.0.d0) then
          mu(i)=1.d0-pi*p**2/norm
        endif
      enddo
      return
      end

      function iofr(c1,c2,c3,c4,r,p)
      implicit none
      real*8 r,p,c1,c2,c3,c4,sig1,sig2,iofr
      sig1=sqrt(sqrt(1.d0-(r-p)**2))
      sig2=sqrt(sqrt(1.d0-(r+p)**2))
      iofr=1.d0-c1*(1.d0+(sig2**5-sig1**5)/5.d0/p/r)
     &         -c2*(1.d0+(sig2**6-sig1**6)/6.d0/p/r)
     &         -c3*(1.d0+(sig2**7-sig1**7)/7.d0/p/r)
     &         -c4*(p**2+r**2)
      return
      end

      subroutine occultnl(rl,c1,c2,c3,c4,b0,mulimb0,mulimbf,nb)
c; Please cite Mandel & Agol (2002) if making use of this routine.
      implicit none
      integer i,j,nb,nr,i1,i2,nmax
      parameter (nmax=2**16)
      real*8 mulimbf(nb,5),pi,c1,c2,c3,c4,rl,bt0(nb),b0(nb),mulimb0(nb),
     &       mulimb(nb),mulimbp(nb),dt,t(nmax+1),th(nmax+1),r(nmax+1),
     &       sig,mulimb1(nb),mulimbhalf(nb),mulimb3half(nb),mulimb2(nb),
     &       sig1,sig2,omega,dmumax,fac,mu(nb),f1,f2
      pi=acos(-1.d0)
C  This routine uses the results for a uniform source to
C  compute the lightcurve for a limb-darkened source
C  (5-1-02 notes)
C Input:
C   rl        radius of the lens   in units of the source radius
C   c1-c4     limb-darkening coefficients
C   b0        impact parameter normalized to source radius
C Output:
C  mulimb0 limb-darkened magnification
C  mulimbf lightcurves for each component
C  
C  First, make grid in radius:
C  Call magnification of uniform source:
      call occultuniform(b0,rl,mulimb0,nb)
      i1=nb
      i2=1
      fac=0.d0
      do i=1,nb
        bt0(i)=b0(i)
        mulimbf(i,1)=1.d0
        mulimbf(i,2)=0.8d0
        mulimbf(i,3)=2.d0/3.d0
        mulimbf(i,4)=4.d0/7.d0
        mulimbf(i,5)=0.5d0
        mulimb(i)=mulimb0(i)
        if(mulimb0(i).ne.1.d0) then
          i1=min(i1,i)
          i2=max(i2,i)
        endif
        fac=max(fac,abs(mulimb0(i)-1.d0))
      enddo
C print,rl
      omega=4.*((1.d0-c1-c2-c3-c4)/4.+c1/5.+c2/6.+c3/7.+c4/8.)
      nr=2
      dmumax=1.d0
c      write(6,*) 'i1,i2 ',i1,i2
      do while (dmumax.gt.fac*1.d-3)
        do i=i1,i2
          mulimbp(i)=mulimb(i)
        enddo
        nr=nr*2
c        write(6,*) 'nr ',nr
        dt=0.5d0*pi/dble(nr)
        if(nr+1.gt.nmax) write(0,*) "M&A Seg: ",nr+1,b0(1)
        do j=1,nr+1
          t(j) =dt*dble(j-1)
          th(j)=t(j)+0.5d0*dt
          r(j)=sin(t(j))
        enddo
        sig=sqrt(cos(th(nr)))
        do i=i1,i2
          mulimbhalf(i) =sig**3*mulimb0(i)/(1.d0-r(nr))
          mulimb1(i)    =sig**4*mulimb0(i)/(1.d0-r(nr))
          mulimb3half(i)=sig**5*mulimb0(i)/(1.d0-r(nr))
          mulimb2(i)    =sig**6*mulimb0(i)/(1.d0-r(nr))
        enddo
        do j=2,nr
          do i=1,nb
            b0(i)=bt0(i)/r(j)
          enddo
C  Calculate uniform magnification at intermediate radii:
          call occultuniform(b0,rl/r(j),mu,nb)
C  Equation (29):
          sig1=sqrt(cos(th(j-1)))
          sig2=sqrt(cos(th(j)))
          dmumax=0.d0
          do i=i1,i2
            f1=r(j)*r(j)*mu(i)/(r(j)-r(j-1))
            f2=r(j)*r(j)*mu(i)/(r(j+1)-r(j))
            mulimbhalf(i) =mulimbhalf(i) +f1*sig1**3-f2*sig2**3
            mulimb1(i)    =mulimb1(i)    +f1*sig1**4-f2*sig2**4
            mulimb3half(i)=mulimb3half(i)+f1*sig1**5-f2*sig2**5
            mulimb2(i)    =mulimb2(i)    +f1*sig1**6-f2*sig2**6
            mulimb(i)=((1.d0-c1-c2-c3-c4)*mulimb0(i)+c1*mulimbhalf(i)*dt
     &        +c2*mulimb1(i)*dt+c3*mulimb3half(i)*dt+c4*mulimb2(i)*dt)
     &        /omega
            if(mulimb(i)+mulimbp(i).ne.0.d0) then 
              dmumax=max(dmumax,abs(mulimb(i)-mulimbp(i))/(mulimb(i)+
     &               mulimbp(i)))
            endif
          enddo
        enddo
      enddo
      do i=i1,i2
        mulimbf(i,1)=mulimb0(i)
        mulimbf(i,2)=mulimbhalf(i)*dt
        mulimbf(i,3)=mulimb1(i)*dt
        mulimbf(i,4)=mulimb3half(i)*dt
        mulimbf(i,5)=mulimb2(i)*dt
        mulimb0(i)=mulimb(i)
      enddo
      do i=1,nb
        b0(i)=bt0(i)
      enddo
      return
      end
      
      subroutine occultuniform(b0,w,muo1,nb)
      implicit none
      integer i,nb
      real*8 muo1(nb),w,b0(nb),z,pi,lambdae,kap0,kap1
      if(abs(w-0.5d0).lt.1.d-3) w=0.5d0
      pi=acos(-1.d0)
C  This routine computes the lightcurve for occultation
C  of a uniform source without microlensing  (Mandel & Agol 2002).
C Input:
C 
C  rs   radius of the source (set to unity)
C  b0   impact parameter in units of rs
C  w    occulting star size in units of rs
C 
C Output:
C  muo1 fraction of flux at each b0 for a uniform source
C 
C  Now, compute pure occultation curve:
      do i=1,nb
C  substitute z=b0(i) to shorten expressions
        z=b0(i)
C  the source is unocculted:
C  Table 3, I.
        if(z.ge.1.d0+w) then
          muo1(i)=1.d0
          goto 1
        endif
C  the  source is completely occulted:
C  Table 3, II.
        if(w.ge.1.d0.and.z.le.w-1.d0) then
          muo1(i)=0.d0
          goto 1
        endif
C  the source is partly occulted and the occulting object crosses the limb:
C  Equation (26):
        if(z.ge.abs(1.d0-w).and.z.le.1.d0+w) then
          kap1=acos(min((1.d0-w*w+z*z)/2.d0/z,1.d0))
          kap0=acos(min((w*w+z*z-1.d0)/2.d0/w/z,1.d0))
          lambdae=w*w*kap0+kap1
          lambdae=(lambdae-0.5d0*sqrt(max(4.d0*z*z-(1.d0+z*z-w*w)**2,
     &            0.d0)))/pi
          muo1(i)=1.d0-lambdae
        endif
C  the occulting object transits the source star (but doesn't
C  completely cover it):
        if(z.le.1.d0-w) muo1(i)=1.d0-w*w
 1      continue
      enddo
C muo1=1.d0-lambdae
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function distance(asep,eccn,Tanom)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision asep,eccn,Tanom
      
      distance=asep*(1.0d0-eccn*eccn)/(1+eccn*cos(Tanom))
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function trueanomaly(eccn,Eanom)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision eccn,Eanom,temp(2)
      
      temp(1)=sqrt((1.0d0+eccn)/(1.0d0-eccn))
      temp(2)=tan(Eanom/2.)
      trueanomaly=2.0d0*atan(temp(1)*temp(2))
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine kepler(Manom,Eanom,eccn)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer i,itmax
      parameter(itmax=100)
      double precision Manom,Eanom,Eold,eccn,diff,thres
      thres=1.0d-8

      Eold=Eanom
      Eanom=Manom+eccn*sin(Eanom)
      diff=abs(1.0d0-Eanom/Eold)
      Eold=Eanom
      i=0
      do while ((diff.gt.thres).and.(i.lt.itmax))
        Eanom=Manom+eccn*sin(Eanom)
        diff=abs(1.0d0-Eanom/Eold)
        Eold=Eanom
        i=i+1
      enddo

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine invkepler(Eanom,Manom,eccn)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer i,itmax
      parameter(itmax=100)
      double precision Manom,Eanom,Mold,eccn,diff,thres
      thres=1.0d-6

      Mold=Manom
      Manom=Eanom-eccn*sin(Manom)
      diff=abs(1.0d0-Manom/Mold)
      Mold=Manom
      i=0
      do while ((diff.gt.thres).and.(i.lt.itmax))
        Manom=Eanom-eccn*sin(Manom)
        diff=abs(1.0d0-Manom/Mold)
        Mold=Manom
        i=i+1
      enddo
c      if(i.ge.itmax) write(0,*) "invkepler itmax"

      return
      end

      subroutine occultquad(z0,u1,u2,p,muo1,mu0,nz)
C  This routine computes the lightcurve for occultation
C  of a quadratically limb-darkened source without microlensing.
C  Please cite Mandel & Agol (2002) if you make use of this routine
C  in your research.  Please report errors or bugs to agol@tapir.caltech.edu
      implicit none
      integer i,nz
      double precision z0(nz),u1,u2,p,muo1(nz),mu0(nz),
     &       mu(nz),lambdad(nz),etad(nz),lambdae(nz),lam,
     &       pi,x1,x2,x3,z,omega,kap0,kap1,q,Kk,Ek,Pk,n,ellec,ellk,rj
      if(abs(p-0.5d0).lt.1.d-3) p=0.5d0
C
C Input:
C
C rs   radius of the source (set to unity)
C z0   impact parameter in units of rs
C p    occulting star size in units of rs
C u1   linear    limb-darkening coefficient (gamma_1 in paper)
C u2   quadratic limb-darkening coefficient (gamma_2 in paper)
C
C Output:
C
C muo1 fraction of flux at each z0 for a limb-darkened source
C mu0  fraction of flux at each z0 for a uniform source
C
C Limb darkening has the form:
C  I(r)=[1-u1*(1-sqrt(1-(r/rs)^2))-u2*(1-sqrt(1-(r/rs)^2))^2]/(1-u1/3-u2/6)/pi
C 
C To use this routine
C
C Now, compute pure occultation curve:
      omega=1.d0-u1/3.d0-u2/6.d0
      pi=acos(-1.d0)
C Loop over each impact parameter:
      do i=1,nz
C substitute z=z0(i) to shorten expressions
        z=z0(i)
        x1=(p-z)**2
        x2=(p+z)**2
        x3=p**2-z**2
C the source is unocculted:
C Table 3, I.
        if(z.ge.1.d0+p) then
          lambdad(i)=0.d0
          etad(i)=0.d0
          lambdae(i)=0.d0
          goto 10
        endif
C the  source is completely occulted:
C Table 3, II.
        if(p.ge.1.d0.and.z.le.p-1.d0) then
          lambdad(i)=1.d0
          etad(i)=1.d0
          lambdae(i)=1.d0
          goto 10
        endif
C the source is partly occulted and the occulting object crosses the limb:
C Equation (26):
        if(z.ge.abs(1.d0-p).and.z.le.1.d0+p) then
          kap1=acos(min((1.d0-p*p+z*z)/2.d0/z,1.d0))
          kap0=acos(min((p*p+z*z-1.d0)/2.d0/p/z,1.d0))
          lambdae(i)=p*p*kap0+kap1
          lambdae(i)=(lambdae(i)-0.5d0*sqrt(max(4.d0*z*z-
     &               (1.d0+z*z-p*p)**2,0.d0)))/pi
        endif
C the occulting object transits the source star (but doesn't
C completely cover it):
        if(z.le.1.d0-p) lambdae(i)=p*p
C the edge of the occulting star lies at the origin- special 
C expressions in this case:
        if(abs(z-p).lt.1.d-4*(z+p)) then
C Table 3, Case V.:
          if(z.ge.0.5d0) then
            lam=0.5d0*pi
            q=0.5d0/p
            Kk=ellk(q)
            Ek=ellec(q)
C Equation 34: lambda_3
            lambdad(i)=1.d0/3.d0+16.d0*p/9.d0/pi*(2.d0*p*p-1.d0)*Ek-
     &                 (32.d0*p**4-20.d0*p*p+3.d0)/9.d0/pi/p*Kk
C Equation 34: eta_1
            etad(i)=1.d0/2.d0/pi*(kap1+p*p*(p*p+2.d0*z*z)*kap0-
     &              (1.d0+5.d0*p*p+z*z)/4.d0*sqrt((1.d0-x1)*(x2-1.d0)))
            if(p.eq.0.5d0) then
C Case VIII: p=1/2, z=1/2
              lambdad(i)=1.d0/3.d0-4.d0/pi/9.d0
              etad(i)=3.d0/32.d0
            endif
            goto 10
          else
C Table 3, Case VI.:
            lam=0.5d0*pi
            q=2.d0*p
            Kk=ellk(q)
            Ek=ellec(q)
C Equation 34: lambda_4
            lambdad(i)=1.d0/3.d0+2.d0/9.d0/pi*(4.d0*(2.d0*p*p-1.d0)*Ek+
     &                 (1.d0-4.d0*p*p)*Kk)
C Equation 34: eta_2
            etad(i)=p*p/2.d0*(p*p+2.d0*z*z)
            goto 10
          endif
        endif
C the occulting star partly occults the source and crosses the limb:
C Table 3, Case III:
        if((z.gt.0.5d0+abs(p-0.5d0).and.z.lt.1.d0+p).or.(p.gt.0.5d0.
     &      and.z.gt.abs(1.d0-p)*1.0001d0.and.z.lt.p)) then
          lam=0.5d0*pi
          q=sqrt((1.d0-(p-z)**2)/4.d0/z/p)
          Kk=ellk(q)
          Ek=ellec(q)
          n=1.d0/x1-1.d0
          Pk=Kk-n/3.d0*rj(0.d0,1.d0-q*q,1.d0,1.d0+n)
C Equation 34, lambda_1:
          lambdad(i)=1.d0/9.d0/pi/sqrt(p*z)*(((1.d0-x2)*(2.d0*x2+
     &        x1-3.d0)-3.d0*x3*(x2-2.d0))*Kk+4.d0*p*z*(z*z+
     &        7.d0*p*p-4.d0)*Ek-3.d0*x3/x1*Pk)
          if(z.lt.p) lambdad(i)=lambdad(i)+2.d0/3.d0
C Equation 34, eta_1:
          etad(i)=1.d0/2.d0/pi*(kap1+p*p*(p*p+2.d0*z*z)*kap0-
     &          (1.d0+5.d0*p*p+z*z)/4.d0*sqrt((1.d0-x1)*(x2-1.d0)))
          goto 10
        endif
C the occulting star transits the source:
C Table 3, Case IV.:
        if(p.le.1.d0.and.z.le.(1.d0-p)*1.0001d0) then
          lam=0.5d0*pi
          q=sqrt((x2-x1)/(1.d0-x1))
          Kk=ellk(q)
          Ek=ellec(q)
          n=x2/x1-1.d0
          Pk=Kk-n/3.d0*rj(0.d0,1.d0-q*q,1.d0,1.d0+n)
C Equation 34, lambda_2:
          lambdad(i)=2.d0/9.d0/pi/sqrt(1.d0-x1)*((1.d0-5.d0*z*z+p*p+
     &         x3*x3)*Kk+(1.d0-x1)*(z*z+7.d0*p*p-4.d0)*Ek-3.d0*x3/x1*Pk)
          if(z.lt.p) lambdad(i)=lambdad(i)+2.d0/3.d0
          if(abs(p+z-1.d0).le.1.d-4) then
            lambdad(i)=2/3.d0/pi*acos(1.d0-2.d0*p)-4.d0/9.d0/pi*
     &            sqrt(p*(1.d0-p))*(3.d0+2.d0*p-8.d0*p*p)
          endif
C Equation 34, eta_2:
          etad(i)=p*p/2.d0*(p*p+2.d0*z*z)
        endif
 10     continue
C Now, using equation (33):
        muo1(i)=1.d0-((1.d0-u1-2.d0*u2)*lambdae(i)+(u1+2.d0*u2)*
     &      lambdad(i)+u2*etad(i))/omega
C Equation 25:
        mu0(i)=1.d0-lambdae(i)
      enddo
      return
      end

      FUNCTION rc(x,y)
      REAL*8 rc,x,y,ERRTOL,TINY,SQRTNY,BIG,TNBG,COMP1,COMP2,THIRD,C1,C2,
     *C3,C4
      PARAMETER (ERRTOL=.04d0,TINY=1.69d-38,SQRTNY=1.3d-19,BIG=3.d37,
     *TNBG=TINY*BIG,COMP1=2.236d0/SQRTNY,COMP2=TNBG*TNBG/25.d0,
     *THIRD=1.d0/3.d0,C1=.3d0,C2=1.d0/7.d0,C3=.375d0,C4=9.d0/22.d0)
      REAL*8 alamb,ave,s,w,xt,yt
      if(x.lt.0..or.y.eq.0..or.(x+abs(y)).lt.TINY.or.(x+
     *abs(y)).gt.BIG.or.(y.lt.-COMP1.and.x.gt.0..and.x.lt.COMP2))pause 
     *'invalid arguments in rc'
      if(y.gt.0.d0)then
        xt=x
        yt=y
        w=1.
      else
        xt=x-y
        yt=-y
        w=sqrt(x)/sqrt(xt)
      endif
1     continue
        alamb=2.d0*sqrt(xt)*sqrt(yt)+yt
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        ave=THIRD*(xt+yt+yt)
        s=(yt-ave)/ave
      if(abs(s).gt.ERRTOL)goto 1
      rc=w*(1.d0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software

      FUNCTION rj(x,y,z,p)
      REAL*8 rj,p,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6,C7,C8
      PARAMETER (ERRTOL=.05d0,TINY=2.5d-13,BIG=9.d11,C1=3.d0/14.d0,
     *C2=1.d0/3.d0,C3=3.d0/22.d0,C4=3.d0/26.d0,C5=.75d0*C3,
     *C6=1.5d0*C4,C7=.5d0*C2,C8=C3+C3)
CU    USES rc,rf
      REAL*8 a,alamb,alpha,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,ed,
     *ee,fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt,rc,rf
      if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z,abs(p)).lt.TINY.or.max(x,y,
     *z,abs(p)).gt.BIG)pause 'invalid arguments in rj'
      sum=0.d0
      fac=1.d0
      if(p.gt.0.d0)then
        xt=x
        yt=y
        zt=z
        pt=p
      else
        xt=min(x,y,z)
        zt=max(x,y,z)
        yt=x+y+z-xt-zt
        a=1.d0/(yt-p)
        b=a*(zt-yt)*(yt-xt)
        pt=yt+b
        rho=xt*zt/yt
        tau=p*pt/yt
        rcx=rc(rho,tau)
      endif
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
        beta=pt*(pt+alamb)**2
        sum=sum+fac*rc(alpha,beta)
        fac=.25d0*fac
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        zt=.25d0*(zt+alamb)
        pt=.25d0*(pt+alamb)
        ave=.2d0*(xt+yt+zt+pt+pt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        delp=(ave-pt)/ave
      if(max(abs(delx),abs(dely),abs(delz),abs(delp)).gt.ERRTOL)goto 1
      ea=delx*(dely+delz)+dely*delz
      eb=delx*dely*delz
      ec=delp**2
      ed=ea-3.d0*ec
      ee=eb+2.d0*delp*(ea-ec)
      rj=3.d0*sum+fac*(1.d0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*
     *(-C8+delp*C4))+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave))
      if (p.le.0.d0) rj=a*(b*rj+3.d0*(rcx-rf(xt,yt,zt)))
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software

      function ellec(k)
      implicit none
      double precision k,m1,a1,a2,a3,a4,b1,b2,b3,b4,ee1,ee2,ellec
C Computes polynomial approximation for the complete elliptic
C integral of the second kind (Hasting's approximation):
      m1=1.d0-k*k
      a1=0.44325141463d0
      a2=0.06260601220d0
      a3=0.04757383546d0
      a4=0.01736506451d0
      b1=0.24998368310d0
      b2=0.09200180037d0
      b3=0.04069697526d0
      b4=0.00526449639d0
      ee1=1.d0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
      ee2=m1*(b1+m1*(b2+m1*(b3+m1*b4)))*log(1.d0/m1)
      ellec=ee1+ee2
      return
      end

      function ellk(k)
      implicit none
      double precision a0,a1,a2,a3,a4,b0,b1,b2,b3,b4,ellk,
     &       ek1,ek2,k,m1
C Computes polynomial approximation for the complete elliptic
C integral of the first kind (Hasting's approximation):
      m1=1.d0-k*k
      a0=1.38629436112d0
      a1=0.09666344259d0
      a2=0.03590092383d0
      a3=0.03742563713d0
      a4=0.01451196212d0
      b0=0.5d0
      b1=0.12498593597d0
      b2=0.06880248576d0
      b3=0.03328355346d0
      b4=0.00441787012d0
      ek1=a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
      ek2=(b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*log(m1)
      ellk=ek1-ek2
      return
      end

      FUNCTION rf(x,y,z)
      REAL*8 rf,x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
      PARAMETER (ERRTOL=.08d0,TINY=1.5d-38,BIG=3.d37,THIRD=1.d0/3.d0,
     *C1=1.d0/24.d0,C2=.1d0,C3=3.d0/44.d0,C4=1.d0/14.d0)
      REAL*8 alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
      if(min(x,y,z).lt.0.d0.or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,
     *z).gt.BIG)pause 'invalid arguments in rf'
      xt=x
      yt=y
      zt=z
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        zt=.25d0*(zt+alamb)
        ave=THIRD*(xt+yt+zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      e2=delx*dely-delz**2
      e3=delx*dely*delz
      rf=(1.d0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0NL&WR2.
