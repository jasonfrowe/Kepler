      program starplot5
      implicit none
      integer nmax,nfit,npt,iargc,nunit,nplanet,i,nplanet2,nplanetmax,j,
     .  col,k,order(10),ii,nplot,npt2,bins,iii,Teff
      parameter (nmax=2000000,nfit=108,nplanetmax=10,nplot=1000)
      integer dtype(nmax),dtype2(nmax),flag(nmax),ntype(nplanetmax)
      real px(nmax),py(nmax),rbb(4),pz(nmax),srad,per(nplanetmax)
      double precision sol(nfit),time(nmax),tmodel(nmax),
     .  flux(nmax),ferr(nmax),itime(nmax),Keplertime,serr(nfit,2),
     .  err(nfit),sol2(nfit),tdur(nplanetmax),transitdur,tdurmax,
     .  transitdepth,tdepth(nplanetmax),tdepthmax,toff,ph1,ph2,
     .  phase(nmax),period,mmax,mmin,width,time2(nmax),itime2(nmax),t1,
     .  t2,ratio,bx(2),phase2(nmax),flux2(nmax),ferr2(nmax),logg,
     .  b(nplanetmax)
C     TT variations
      integer ntt(nplanetmax),ntt2(nplanetmax)
      double precision tobs(nplanetmax,nmax),omc(nplanetmax,nmax),
     .  tobs2(nplanetmax,nmax),omc2(nplanetmax,nmax),tcor(nmax),ttcor
      character*80 inputsol,obsfile,ttfile,cline
c      data order /4,2,1,3,5,6,7,8,9,10/
      data order /1,2,3,4,5,6,7,8,9,10/
      
      if(iargc().lt.5) goto 901
      
C     Parse the name of the observations data file from the commandline
      call getarg(1,obsfile)
      nunit=10
      open(unit=nunit,file=obsfile,status='old',err=903)
c      call readdata(nunit,nmax,npt,time,flux,ferr,itime,Keplertime)
      call readkeplc(nunit,nmax,npt,time,flux,ferr,itime,Keplertime)
      close(nunit)!release unit number as we are done with the file.
      
      mmin=flux(1)
      mmax=flux(1)
      do 15 i=2,npt
        mmin=min(mmin,flux(i))
        mmax=max(mmax,flux(i))
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

      call getarg(3,cline)
      read(cline,*) srad
      if(srad.le.0.0) goto 901
      call getarg(4,cline)
      read(cline,*) Teff
      if(Teff.le.0) goto 901
      call getarg(5,cline)
      read(cline,*) logg
      if(logg.le.0.0d0) goto 901

      do 22 i=1,nplanet
        if(iargc().ge.5+i)then
            call getarg(5+i,ttfile)

            if(ttfile.eq.'null')then
                ntt(i)=0
            else
                nunit=10
                open(unit=nunit,file=ttfile,status='old',err=905)
                call readttfile(nunit,nplanetmax,nmax,i,ntt,tobs,omc)
                close(nunit)
            endif
            
        else
            ntt(i)=0
        endif
c        write(0,*) "ntt",i,ntt(i)
 22   continue

      read(5,*) (ntype(i),i=1,nplanetmax)
c      write(0,*) (ntype(i),i=1,nplanetmax)
      read(5,*) (b(i),i=1,nplanetmax)
c      write(0,*) "b: ",(b(i),i=1,nplanetmax)

      call pgopen('?') !open PGPlot device
      call pgask(.true.) !don't ask for new page.. just do it.
c      call PGPAP (10.0/real(nplanet+1) ,real(nplanet+1)/2.0) !paper size
      call PGPAP(8.0/real(nplanet+2)*2.0 ,real(nplanet+2)/2.0)
c      call pgpap(8.0,0.7)
      call pgslw(3)
      call pgsubp(1,nplanet+2)  !break up plot into grid

      do 26 i=1,nplanet
         col=10*(i-1)
         per(i)=real(sol(10+col))
 26   continue
      call rqsort(nplanet,per,order)

CCCCCCCCC
C     Do star plot
CCCCCCCCC
      call pgpage()
      call pgvport(0.25,0.95,0.0,1.0)
      call starplot(nplanet,nfit,sol,srad,Teff,logg,ntype,b,order)

 
      do 10 i=1,nplanet
        tdur(i)=1.0d0*transitdur(i,nfit,sol)/8.64d4
        tdepth(i)=transitdepth(i,nfit,sol,sol2,ntt,tobs,omc)
c        write(6,*) "tdur:",tdur(i),tdepth(i)
 10   continue
      
      tdurmax=tdur(1)
      tdepthmax=tdepth(1)
      do 11 i=2,nplanet
        tdurmax=max(tdur(i),tdurmax)
        tdepthmax=max(tdepth(i),tdepthmax)
 11   continue
      tdurmax=max(tdurmax,2.0/24.0)

      rbb(3)=real(mmin)
      rbb(4)=real(mmax)

      do 12 ii=1,nplanet
        i=order(ii)
        col=10*(i-1)

        do 24 j=1,npt
            call lininterp(tobs,omc,nplanetmax,nmax,i,ntt,time(j),
     .          ttcor)
            tcor(j)=time(j)-ttcor
 24     continue

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
        call transitmodel(nfit,nplanet,nplanetmax,sol2,nmax,npt,time,
     .      itime,ntt,tobs,omc,tmodel,dtype)
        
        toff=0.5-(sol(9+col)/sol(10+col)-int(sol(9+col)/sol(10+col)))
        if(toff.lt.0.0)toff=toff+1.0
        ph1=0.5-1.0d0*tdurmax/sol(10+col)
        if(ph1.lt.0.25)ph1=0.25
        ph2=0.5+1.0d0*tdurmax/sol(10+col)
        if(ph2.gt.0.75)ph2=0.75
        
        rbb(1)=real(ph1)
        rbb(2)=real(ph2)
        
        call phasept(npt,tcor,phase,period,toff,ratio)
        
        k=0
        do 14 j=1,npt
            if((phase(j).gt.ph1).and.(phase(j).lt.ph2))then
                k=k+1
                px(k)=real(phase(j))
                py(k)=real(flux(j)-tmodel(j))+1.0
c                write(0,*) px(k),py(k)
            endif
 14     continue
        call pgpage()
        call pgvport(0.25,0.95,0.0,1.0)
        call pgwindow(rbb(1),rbb(2),rbb(3),rbb(4))
c        write(6,*) rbb(1),rbb(2),rbb(3),rbb(4)
        call pgpt(k,px,py,17)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Plot the binned observations
        k=0
        do 6 j=1,npt
                if((real(flux(j)).gt.rbb(3)).and.(real(flux(j)).lt.
     .            rbb(4)))then
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
        if(ntype(i).eq.2)then
         call pgsci(4)
        else
         call pgsci(8)
        endif
c        call pgsch(1.5)
        call pgpt(npt2,px,py,17)
c        call pgsch(0.8)
        call pgerrb(6,npt2,px,py,pz,1.0)
        call pgsci(1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        do 16 j=1,8
            sol2(j)=sol(j)
  16    continue
        do 18 j=9,18          
            sol2(j)=sol(j+col)
  18    continue
        
        ntt2(1)=ntt(i)
        do 23 j=1,ntt(i)
            tobs2(1,j)=tobs(i,j)
            omc2(1,j)=omc(i,j)
 23     continue
        
        call lininterp(tobs2,omc2,nplanetmax,nmax,1,ntt2,sol2(9+col),
     .        ttcor)
C       We reverse correct transit-time to match model
        t1=sol2(9+col)-1.5*tdurmax+ttcor
        t2=sol2(9+col)+1.5*tdurmax+ttcor
        do 19 k=1,nplot
            time2(k)=t1+(t2-t1)/dble(nplot)*dble(k-1)
            itime2(k)=itime(1)
            dtype2(k)=0
 19     continue
c        call transitmodel(nfit,1,sol2,nplot,time2,itime2,tmodel,dtype2)
        call transitmodel(nfit,1,nplanetmax,sol2,nmax,nplot,time2,
     .      itime2,ntt2,tobs2,omc2,tmodel,dtype2)

        do 25 j=1,nplot!npt
            ttcor=0.0
            call lininterp(tobs2,omc2,nplanetmax,nmax,1,ntt2,time2(j),
     .          ttcor)
            tcor(j)=time2(j)-ttcor
 25     continue     
     
        call phasept(nplot,tcor,phase,period,toff,ratio)
        do 20 k=1,nplot
            px(k)=real(phase(k))
            py(k)=real(tmodel(k))
c            write(0,*) px(k),py(k)
 20     continue
c        read(5,*)
        call sort2(nplot,px,py)
        call pgsci(i+1)
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
        call pgsch(1.0)
        
c        read(5,*)
        
 12   continue
      
      call pgclos()
      
      goto 999
 901  write(0,*) "Usage: transitplot5 <photfile> <fitpars> <Rstar> <Teff
     .> <logg>"
      goto 999
 902  write(0,*) "Error opening ",inputsol
      goto 999
 903  write(0,*) "Error opening ",obsfile
      goto 999
 905  write(0,*) "Error opening ",ttfile
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
      parameter(nmax=2000000)
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
      double precision function transitdepth(np,nfit,sol,sol2,ntt,tobs,
     .  omc)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer np,nfit,col,i,nplanetmax,nmax
      parameter(nplanetmax=10,nmax=2000000)
      integer ntt(nplanetmax),ntt2(nplanetmax)
      double precision sol(nfit),itime(1),dtype(1),sol2(nfit),tmodel(1),
     .  time(1),tobs(nplanetmax,nmax),omc(nplanetmax,nmax),
     .  tobs2(nplanetmax,nmax),omc2(nplanetmax,nmax),ttcor
      
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

      ntt2(1)=ntt(np)
      do 23 i=1,ntt(np)
        tobs2(1,i)=tobs(np,i)
        omc2(1,i)=omc(np,i)
 23   continue
      
      time(1)=sol2(9)
      call lininterp(tobs2,omc2,nplanetmax,nmax,1,ntt2,time(1),
     .  ttcor)
      time(1)=time(1)+ttcor !we reverse correct the OC to get center of
                            !transit time
      
c      call transitmodel(nfit,1,sol2,1,time,itime,tmodel,dtype)
      call transitmodel(nfit,1,nplanetmax,sol2,nmax,1,time,
     .      itime,ntt2,tobs2,omc2,tmodel,dtype)
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
c      write(0,*) bb,adrs,rdr
        
      temp(1)=Psec/Pi
      temp(2)=1.0d0/adrs
      temp(3)=(1+rdr)**2.0-bb
      temp(4)=1-cincl*cincl
C     Transit duration in days
      transitdur=temp(1)*asin(temp(2)*sqrt(temp(3)/temp(4)))
      
      return
      end
