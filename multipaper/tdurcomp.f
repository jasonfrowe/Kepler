      program tdurcomp
      implicit none
      integer nmax,n,i,j,ndata,nunit,k,np,seed,nhavetransit,ieccn,imass,
     .   mnr,ii,ii2,nplot,flag,niter,nc,nitermax,nfit,nbuffer,nmov,nup,
     .   npars,nupdate,nupcor,nb,nacor,nacorsub,naprob,naprobsub,nsel,
     .   nsel2,nmc,nsum,nbin,skip,nbstar
      parameter(nmax=5000,ndata=26,nitermax=1000000,nfit=5,nbuffer=100)
      integer ic(nmax),isym(nmax),now(3),dumi(4),ia(nfit),ngcor(nfit),
     .   ngcorsub(nfit),ngprob(nfit),ngprobsub(nfit),nmp(6),nmp2(6),
     .   ngs(6),nas
      real data(nmax,ndata),rhoavg,x(nmax),y(nmax),z(nmax),miny,maxy,
     .   rhostar,rhostarerr,Per,adrs,G,Pi,ran2,eccn,w,Eanom,Tanom,
     .   trueanomaly,distance,ddrs,TransitProb,rTransit,x2(nmax),
     .   y2(nmax),z2(nmax),rmax,eccnsig,Schwarz,r1,r2,dumr(9),Rsun,
     .   mass(nmax),Msun,mmax,m1,m2,kroupa,mbinary,rbinary,rhobin,
     .   mrrel(nmax,4),massradfn,binfrac,rbin,minmass,dumr2,gasdev,
     .   rhoavgmeas(nmax),bx2(nmax),bx(nmax),comppop,rchi,bdx1(nmax),
     .   bdy1(nmax),bdx2(nmax),bdy2(nmax),chi1,chi2,fratio,u,alpha,
     .   esig,bsig,dsig,dbinf,eccnsign,binfracn,rpick,mc(nitermax,nfit),
     .   b(nitermax),eave,evar,estd,bave,bvar,bstd,ans(nfit),ansn(nfit),
     .   asig(nfit),buffer(nfit,nbuffer),gscale(nfit),mcmctype,corscale,
     .   gcorsub,corsub,mp(nmax,6),bx3(nmax),mp2(nmax,6),datamin,
     .   datamax,bmax,bdy3(nmax),sig
      character*80 multiprops,dumc,spars,massrad,cout

      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      Pi=acos(-1.d0)!define Pi and 2*Pi
      Rsun=6.9599E8
      Msun=1.989E30

      nbin=30
      datamin= 0.0
      datamax= 1.2
      nsum=100

      do 18 i=1,nmax
         mass(i)=0.0
         rhoavgmeas(i)=0.0
 18   continue

      do 54 i=1,6
         nmp(i)=0
         nmp2(i)=0
 54   continue

      multiprops="multiprops.dat" !file that contains multi-analysis

      nunit=10
C     open up multiprops file that has the multi-planet analysis
      open(unit=nunit,file=multiprops,status='old',err=901)
      read(nunit,*) dumc !first line is a header
      i=1
 14   read(nunit,*,end=15) (data(i,j),j=1,ndata)
c         write(0,*) i,koi(i)  !this line was for debugging
         i=i+1
      goto 14
 15   continue
      close(nunit)
      n=i-1

      spars="spars.dat" !file that contains stellar parameters
      open(unit=nunit,file=spars,status='old',err=902)
      read(nunit,*) dumc !file line is a header
 16   read(nunit,*,end=17) (dumi(i),i=1,4),(dumr(j),j=1,9)
         mass(dumi(1))=10.0**(dumr(1)-2.0)*(dumr(5)*Rsun)**2.0/G/Msun
         if(mass(dumi(1)).le.0.0)then
            mass(dumi(1))=1.0  !replace unclassified with solar
         endif
c         write(0,*) dumr(1),dumr(5),mass(dumi(1))
c         read(5,*)
      goto 16
 17   continue
      close(nunit)

C     Read in one of the Dartmouth tracks as a M-R relationship
      massrad="massrad.txt"
      open(unit=nunit,file=massrad,status='old',err=903)
      read(nunit,*) dumc !first line is a header
      mnr=1
 19   read(nunit,*,end=34) dumi(1),(mrrel(mnr,i),i=1,4)
c         write(0,*) "x",(mrrel(mnr,i),i=1,4)
         mnr=mnr+1
      goto 19
 34   continue
      mnr=mnr-1
      close(nunit)

c      call pgopen('?')
c      call pgopen('/xserve')
      call pgopen('/null')
      call pgask(.false.)
      call pgpage()
      call PGPAP ( 8.0 ,1.0)
      call pgsubp(2,3)
      call pgpanl(1,1)
      call pgvport(0.15,0.85,0.2,0.9)
      call pgwindow(-2.0,log10(20.0),-0.1,1.2)
      call pgslw(3)
      call pgsch(1.5)
      CALL PGBOX('BCLNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("<\(2143)> (e=0)","|(\(2143)\di\u-\(2143)\dj\u)|/(\(2
     .143)\di\u+\(2143)\dj\u) (e=0)","")

      ii=0
      miny= 99.9e30
      maxy=-99.9e30
      do 10 i=1,nmax
         rhoavg=0.0
         j=0
         do 12 k=1,n
            if((i.eq.int(data(k,1))).and.(data(k,21).ge.7.1).and.
     .       (data(k,24).gt.1))then
               j=j+1
c               rhoavg=rhoavg+data(k,10)
               x(j)=data(k,10) !<rhostar>
               y(j)=data(k,10) !<rhostar>
               z(j)=data(k,11) !<rhostar>_sig
            endif
 12      continue
         np=j
         if(np.gt.1)then
c            rhoavg=rhoavg/real(np)
c            rhoavgmeas(i)=rhoavg
            do 11 j=1,np
               do 53 k=j+1,np
                  rhoavg=(y(j)+y(k))/2
                  call pgsci(np)
                  ii=ii+1
                  bx(ii)=(y(j)-y(k))/(y(j)+y(k))
                  bx(ii)=abs(bx(ii))
                  call pgpt1(log10(rhoavg),bx(ii),4)
                  nmp(np)=nmp(np)+1
                  mp(nmp(np),np)=bx(ii) !store histograms for 2,3,4..
 53            continue
 11         continue
c            write(6,*) i,rhoavg,x(2),y(2),np
         endif
 10   continue
      call pgsci(1)
C     bin up the data
C     bdx1 and byd1 contain the observables that we will match to the
C     models..
      call bindata(nbin,ii,bx,bdx1,bdy1,datamin,datamax,bmax)
c      goto 900

      call itime(now)
      seed=abs(now(3)+now(1)*now(2)+now(1)*now(3)+now(2)*now(3)*100)
      dumr(1)=ran2(-seed)

      goto 400 !just do fitting


C     We now need to simulate the eccentric population..
C     Use the stellar properties for <rhostar> and same Periods..
C     Then

      nplot=0
      call pgpage()
      call pgwindow(-0.1,1.1,0.0,150.0)
      call pgslw(3)
c      call pgsch(1.5)
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("|(\(2143)\di\u-\(2143)\dj\u)/(\(2143)\di\u+\(2143)\d
     .j\u)| (e=0)","N","")
      call bindata(nbin,ii,bx,bdx1,bdy1,datamin,datamax,bmax)
      call pgpt(nbin,bdx1,bdy1,17)
      call pgslw(5)
      call pgline(nbin,bdx1,bdy1)
      do 93 i=1,nbin
         call pgerr1(6,bdx1(i),bdy1(i),sqrt(bdy1(i)),1.0)
 93   continue
      call pgslw(3)

      ans(1)=0.0
      ans(2)=0.867
      ans(3)=3.03
      ans(4)=1.0
      ans(5)=-1.0
      do 62 i=1,nbin
         bdy2(i)=0.0
 62   continue
      do 59 i=1,nsum
         call makepop(nfit,ans,seed,nmax,ndata,data,n,mass,
     .      x,y,z,ii2,bx2,dumr,mnr,mrrel,nplot,nmp2,mp2,nbstar)
         call bindata(nbin,ii2,bx2,bdx1,bdy1,datamin,datamax,bmax)
         do 60 j=1,nbin
            bdy2(j)=bdy2(j)+bdy1(j)
 60      continue
 59   continue
      do 61 i=1,nbin
         bdy2(i)=bdy2(i)/real(nsum)
 61   continue
      call pgsci(2)
      call pgpt(nbin,bdx1,bdy2,17)
      call pgline(nbin,bdx1,bdy2)

      ans(1)=0.5
      ans(2)=0.867
      ans(3)=3.03
      ans(4)=1.0
      ans(5)=-1.0
      do 63 i=1,nbin
         bdy2(i)=0.0
 63   continue
      do 64 i=1,nsum
         call makepop(nfit,ans,seed,nmax,ndata,data,n,mass,
     .      x,y,z,ii2,bx2,dumr,mnr,mrrel,nplot,nmp2,mp2,nbstar)
         call bindata(nbin,ii2,bx2,bdx1,bdy1,datamin,datamax,bmax)
         do 65 j=1,nbin
            bdy2(j)=bdy2(j)+bdy1(j)
 65      continue
 64   continue
      do 66 i=1,nbin
         bdy2(i)=bdy2(i)/real(nsum)
 66   continue
      call pgsci(3)
      call pgpt(nbin,bdx1,bdy2,17)
      call pgline(nbin,bdx1,bdy2)

      ans(1)=0.5
      ans(2)=0.867
      ans(3)=3.03
      ans(4)=1.0
      ans(5)=0.0
      do 67 i=1,nbin
         bdy2(i)=0.0
 67   continue
      do 68 i=1,nsum
         call makepop(nfit,ans,seed,nmax,ndata,data,n,mass,
     .      x,y,z,ii2,bx2,dumr,mnr,mrrel,nplot,nmp2,mp2,nbstar)
         call bindata(nbin,ii2,bx2,bdx1,bdy1,datamin,datamax,bmax)
         do 69 j=1,nbin
            bdy2(j)=bdy2(j)+bdy1(j)
 69      continue
 68   continue
      do 70 i=1,nbin
         bdy2(i)=bdy2(i)/real(nsum)
 70   continue
      call pgsci(4)
      call pgpt(nbin,bdx1,bdy2,17)
      call pgline(nbin,bdx1,bdy2)

      ans(1)=0.0
      ans(2)=0.0
      ans(3)=0.0
      ans(4)=1.0
      ans(5)=-1.0
      do 71 i=1,nbin
         bdy2(i)=0.0
 71   continue
      do 72 i=1,nsum
         call makepop(nfit,ans,seed,nmax,ndata,data,n,mass,
     .      x,y,z,ii2,bx2,dumr,mnr,mrrel,nplot,nmp2,mp2,nbstar)
         call bindata(nbin,ii2,bx2,bdx1,bdy1,datamin,datamax,bmax)
         do 73 j=1,nbin
            bdy2(j)=bdy2(j)+bdy1(j)
 73      continue
 72   continue
      do 74 i=1,nbin
         bdy2(i)=bdy2(i)/real(nsum)
 74   continue
      call pgsci(5)
      call pgpt(nbin,bdx1,bdy2,17)
      call pgline(nbin,bdx1,bdy2)

      skip=1
      if(skip.eq.0)then
      ans(1)=0.013
      ans(2)=0.156
      ans(3)=3.097
      ans(4)=1.857
      ans(5)=-1.0
      do 75 i=1,nbin
         bdy2(i)=0.0
 75   continue
      do 76 i=1,nsum
         call makepop(nfit,ans,seed,nmax,ndata,data,n,mass,
     .      x,y,z,ii2,bx2,dumr,mnr,mrrel,nplot,nmp2,mp2,nbstar)
         call bindata(nbin,ii2,bx2,bdx1,bdy1,datamin,datamax,bmax)
         do 77 j=1,nbin
            bdy2(j)=bdy2(j)+bdy1(j)
 77      continue
 76   continue
      do 78 i=1,nbin
         bdy2(i)=bdy2(i)/real(nsum)
c         write(0,*) i, bdy2(i)
 78   continue
      call pgsci(6)
      call pgpt(nbin,bdx1,bdy2,17)
      call pgline(nbin,bdx1,bdy2)
      endif

      goto 900

 400  continue !jump to just do fitting

      call pgpage()
      call pgwindow(-2.0,log10(20.0),-0.1,1.2)
      call pgslw(3)
c      call pgsch(1.5)
      CALL PGBOX('BCLNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("<\(2143)> (e=0)","|(\(2143)\di\u-\(2143)\dj\u)/(\(21
     .43)\di\u+\(2143)\dj\u)| (e=0)","")

c      eccnsig=0.3
c      esig=0.03
c      binfrac=0.05
c      bsig=0.03

      ans(1)=0.1!0.097
      asig(1)=0.03
c      ans(2)=0.809
      ans(2)=0.867
      asig(2)=0.044
c      ans(3)=1.932
      ans(3)=3.03
      asig(3)=0.01
      ans(4)=1.00
      asig(4)=0.01
      ans(5)=-1.0
      asig(5)=0.05

      do 49 i=1,nfit
         ansn(i)=ans(i)
 49   continue

      nsum=100
      nplot=0
      do 83 i=1,nbin
         bdy3(i)=0.0
 83   continue
      do 79 i=1,nsum
         call makepop(nfit,ans,seed,nmax,ndata,data,n,mass,
     .      x,y,z,ii2,bx2,dumr,mnr,mrrel,nplot,nmp2,mp2,nbstar)
         call bindata(nbin,ii2,bx2,bdx2,bdy2,datamin,datamax,bmax)
         do 80 j=1,nbin
            bdy3(j)=bdy3(j)+bdy2(j)
 80      continue
 79   continue
      do 81 i=1,nbin
         bdy3(i)=bdy3(i)/real(nsum)
 81   continue
      chi1=0.0
      do 82 i=1,nbin
         sig=sqrt(bdy1(i)+bdy3(i))
         if(sig.gt.0) then
            chi1=chi1+(bdy1(i)-bdy3(i))*(bdy1(i)-bdy3(i))/(sig*sig)
         endif
 82   continue

c      do 55 i=2,6
c         do 56 j=1,nmp(i)
c            bx3(j)=mp(j,i)
c 56      continue
c         call pgsci(i-1)
c         if(i.eq.2)then
c            call pghist(nmp(i),bx3,-1.2,1.2,30,0)
c         else
c            call pghist(nmp(i),bx3,-1.2,1.2,30,1)
c         endif
c 55   continue
c      call pgsci(1)

c      do 57 i=2,6
c         do 58 j=1,nmp2(i)
c            bx3(j)=mp2(j,i)
c 58      continue
c         call pgsci(i-1)
c         if(i.eq.2)then
c            call pghist(nmp2(i),bx3,-1.2,1.2,30,0)
c         else
c            call pghist(nmp2(i),bx3,-1.2,1.2,30,1)
c         endif
c 57   continue
c      call pgsci(1)

c      goto 900

C     Calculate Chi-square distribution for exoplanet e dist.
      call bindata(nbin,ii,bx,bdx1,bdy1,datamin,datamax,bmax)
      nplot=0
      do 51 i=1,0
         ans(1)=real(i-1)/100.0
         chi1=0.0
         do 84 j=1,nbin
            bdy3(j)=0.0
 84      continue
         do 52 j=1,nsum
            call makepop(nfit,ans,seed,nmax,ndata,data,n,mass,
     .         x,y,z,ii2,bx2,dumr,mnr,mrrel,nplot,nmp2,mp2,nbstar)
            call bindata(nbin,ii2,bx2,bdx2,bdy2,datamin,datamax,bmax)
            do 85 k=1,nbin
               bdy3(k)=bdy3(k)+bdy2(k)
 85         continue
 52      continue
         do 86 j=1,nbin
            bdy3(j)=bdy3(j)/real(nsum)
 86      continue

         chi1=0.0
         do 87 j=1,nbin
            sig=sqrt(bdy1(j)+bdy3(j))
            if(sig.gt.0) then
               chi1=chi1+(bdy1(j)-bdy3(j))*(bdy1(j)-bdy3(j))/(sig*sig)
            endif
 87      continue
         write(6,*) ans(1),chi1
 51   continue
c      goto 900

      call pgpage()
      call pgwindow(0.0,10.0,0.0,1.0)
      call pgslw(3)
c      call pgsch(1.5)
      CALL PGBOX('BCNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("a","Hierarchical Blends","")

      nplot=0
      nsum=10
      niter=1000000

C     Initialize Gibbs scaler
      do 44 i=1,nfit
        gscale(i)=1.0e0
        ngcor(i)=0
        ngcorsub(i)=0
        ngprob(i)=0
        ngprobsub(i)=0
        ngs(i)=0
 44   continue
      corscale=1.0e0
      nacor=0
      nacorsub=0
      naprob=0
      naprobsub=0
      nb=0
      nas=0

      nmov=0
      nup=0
      npars=1
      nupdate=0

      ia(1)=0 !controls which variables are fit.
      ia(2)=0
      ia(3)=0
      ia(4)=0
      ia(5)=0 !slope of binary mass ratio distribution

      call bindata(nbin,ii,bx,bdx1,bdy1,datamin,datamax,bmax)
      do 37 i=1,niter
         nupcor=0

         mcmctype=ran2(seed)

         if((nmov.lt.nbuffer).or.(mcmctype.le.0.5))then
 48         nc=int(ran2(seed)*real(nfit)+1.0e0)
            if(ia(nc).gt.0) goto 48 !only fit fitted variables
            ansn(nc)=ans(nc)+asig(nc)*gasdev(seed)*gscale(nc)
            nmc=0
         else
            nsel=int(ran2(seed)*real(nbuffer-1)+1.0d0)
            nsel2=int(ran2(seed)*real(nbuffer-1)+1.0d0)
            do 46 j=1,nfit
               ansn(j)=ans(j)+(buffer(j,nsel2)-buffer(j,nsel))*corscale
 46         continue
            nmc=1
         endif

         if((ansn(1).lt.0).or.(ansn(1).ge.0.5))then
            flag=1
            goto 47
         endif
         do 39 j=2,3
            if((ansn(j).le.0).or.(ansn(j).gt.10.0))then
               flag=1
               goto 47
            endif
 39      continue
         if(ansn(4).le.0)then
            flag=1
            goto 47
         endif
         if((ansn(5).lt.-2.0).or.(ansn(5).gt.2.0))then
            flag=1
            goto 47
         endif

c         chi2=0.0
c         do 50 j=1,nsum
c            call makepop(nfit,ansn,seed,nmax,ndata,data,n,mass,
c     .         x,y,z,ii2,bx2,dumr,mnr,mrrel,nplot,nmp2,mp2,nbstar)
c            chi2=chi2+comppop(nmax,ii,bx,ii2,bx2,bdx1,bdy1,bdx2,bdy2,
c     .         30,-2.0,5.0)
c 50      continue
c         chi2=chi2/real(nsum)

         do 88 j=1,nbin
            bdy3(j)=0.0
 88      continue
         do 89 j=1,nsum
            call makepop(nfit,ansn,seed,nmax,ndata,data,n,mass,
     .         x,y,z,ii2,bx2,dumr,mnr,mrrel,nplot,nmp2,mp2,nbstar)
            call bindata(nbin,ii2,bx2,bdx2,bdy2,datamin,datamax,bmax)
            do 90 k=1,nbin
               bdy3(k)=bdy3(k)+bdy2(k)
 90         continue
 89      continue
         do 91 j=1,nbin
            bdy3(j)=bdy3(j)/real(nsum)
 91      continue

         chi2=0.0
         do 92 j=1,nbin
            sig=sqrt(bdy1(j)+bdy3(j))
            if(sig.gt.0) then
               chi2=chi2+(bdy1(j)-bdy3(j))*(bdy1(j)-bdy3(j))/(sig*sig)
            endif
 92      continue

c         write(0,*) chi1,chi2,nmc
         fratio=exp(0.5d0*(chi1-chi2))
         u=ran2(seed)
         alpha=min(fratio,1.0d0)
         if(u.le.alpha)then
            flag=0  !accept chain
            chi1=chi2
            do 40 j=1,nfit
               ans(j)=ansn(j)
 40         continue
         else
            flag=1  !reject chain
         endif

 47      continue !line we jump to with bad parameters

C     Fill Buffer (initial fill)
C     Filling buffer
      if((flag.eq.0).and.(nmov.lt.nbuffer))then
        nmov=nmov+1
        do 42 j=1,nfit
            buffer(j,nmov)=ans(j)
 42     continue
        write(cout,502) "nmov:",nmov
c         write(0,502) "nmov:",nmov
 502    format(A5,I4)
        call ovrwrt(cout,2)
      endif

C     Update buffer if accepted
      if((flag.eq.0).and.(nup.ge.npars).and.(nmov.eq.nbuffer).and.
     .   (mcmctype.le.0.5))then
        nup=0
        nb=nb+1
        if(nb.eq.nbuffer)then
            nupcor=1
            nupdate=nupdate+1
c            npars=npars+2  !keep a larger history..
        endif
        if(nb.gt.nbuffer) nb=1
        write(cout,503) "nb: ",nb,nupdate
c         write(0,503) "nb: ",nb,nupdate
 503    format(A4,I4,1X,I4)
        call ovrwrt(cout,2)
        do 43 j=1,nfit
            buffer(j,nb)=ans(j)
 43     continue
      elseif((flag.eq.0).and.(nup.lt.npars).and.(nmov.eq.nbuffer).and.
     .   (mcmctype.le.0.5))then
        nup=nup+1
      endif
cccccccccccccccc
      if((nmov.lt.nbuffer).or.(mcmctype.le.0.5))then
        if(flag.eq.0)then
            ngcor(nc)=ngcor(nc)+1
            ngcorsub(nc)=ngcorsub(nc)+1
        endif
        ngprob(nc)=ngprob(nc)+1
        ngprobsub(nc)=ngprobsub(nc)+1
      else
        if(flag.eq.0)then
            nacor=nacor+1
            nacorsub=nacorsub+1
        endif
        naprob=naprob+1
        naprobsub=naprobsub+1
      endif

      if(nupcor.eq.1)then
        do 45 j=1,nfit
            if((ngprob(j).ne.ngprobsub(j)).and.(ngs(j).eq.0))then
            write(0,*) "ngs:",j,real(ngcorsub(j))/real(ngprobsub(j))
            if((real(ngcorsub(j))/real(ngprobsub(j)).lt.0.2).or.
     .       (real(ngcorsub(j))/real(ngprobsub(j)).gt.0.3))then
                gcorsub=real(ngcor(j)-ngcorsub(j))/
     .              real(ngprob(j)-ngprobsub(j))
                gscale(j)=gscale(j)*(0.75*(gcorsub+0.01)/(0.25*
     .              (1.0e0-gcorsub+0.01)))**0.25e0
            else
               ngs(j)=1 !stop updating scaling factor
            endif
            endif
            ngcorsub(j)=0 !reset counter
            ngprobsub(j)=0 !reset counter
 45     continue
        if((naprob.ne.naprobsub).and.(nas.eq.0))then
        write(0,*) "nas: ",real(nacorsub)/real(naprobsub)
        if((real(nacorsub)/real(naprobsub).lt.0.2).or.(real(nacorsub)/
     .   real(naprobsub).gt.0.3))then
            corsub=real(nacor-nacorsub)/
     .          real(naprob-naprobsub)
            corscale=corscale*(0.75*(corsub+0.01)/(0.25*
     .          (1.0d0-corsub+0.01)))**0.25d0
c            write(0,*) naprob,naprobsub,corsub
c            write(0,*)  nacor,nacorsub
        else
            nas=1 !stop updating scale factor
        endif
        endif
        nacorsub=0 !reset counter
        naprobsub=0 !reset counter
        write(0,*)
        write(0,*) "corscale"
        write(0,504) corscale
        write(0,*) "gscale"
        write(0,504) (gscale(j),j=1,nfit)
        write(0,*)
 504    format(20(F5.1,1X))
      endif

  38     write(6,500) flag,nc,nmc,chi1,(ans(j),j=1,nfit),nbstar
 500     format(I1,2(1X,I1),6(1X,1PE17.10),1X,I5)

         do 41 j=1,nfit
            mc(i,j)=ans(j)
 41      continue
         if(mod(i,10).eq.0)then
            call pgpt1(ans(2),ans(1),17)
c            call avevar(e,i,eave,evar)
c            call avevar(b,i,bave,bvar)
c            write(0,501) eave,sqrt(evar),bave,sqrt(bvar)
 501        format(4(F6.3,1X))
         endif

 37   continue

c      call pgpage()

      goto 900

      call pgpage()
      call pgwindow(-2.0,log10(20.0),-15.0,15.0)
      call pgslw(3)
c      call pgsch(1.5)
      CALL PGBOX('BCLNTS1',0.0,0,'BCNTS1',0.0,0)
      call pglabel("<rhostar> (YY)","d<rhostar>/sig","")

      miny= 99.9e30
      maxy=-99.9e30
      do 20 i=1,nmax
         rhoavg=0.0
         j=0
         do 22 k=1,n
            if((i.eq.int(data(k,1))).and.(data(k,21).ge.7.1))then
               j=j+1
               rhostar=data(k,12)
               rhostarerr=data(k,13)
               x(j)=real(i)    ! KOI
               y(j)=data(k,10) !<rhostar>
               z(j)=data(k,11) !<rhostar>_sig
               if(data(k,24).gt.1)then
                  ic(j)=1
                  isym(j)=4
               else
                  ic(j)=2
                  isym(j)=2
               endif
            endif
 22      continue
         np=j
         if(np.gt.0)then
            do 21 j=1,np
c               write(0,*) y(j),rhoavg
               y(j)=(y(j)-rhostar)/sqrt(z(j)*z(j)+rhostarerr*rhostarerr)
               miny=min(miny,y(j))
               maxy=max(maxy,y(j))
 21         continue
            do 23 j=1,np
               call pgsci(ic(j))
               call pgpt1(log10(rhostar),y(j),isym(j))
 23         continue
c            write(6,*) i,rhoavg,x(2),y(2),np
         endif
 20   continue


 900  continue
      call pgclos()

c      write(6,*) "min/max y: ",miny,maxy

      goto 999
 901  write(0,*) "Cannot open ",multiprops
      goto 999
 902  write(0,*) "Cannot open ",spars
      goto 999
 903  write(0,*) "Cannot open ",massrad
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE avevar(data,n,ave,var)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n,m
      REAL ave,var,data(n)
C Given array data(1:n), returns its mean as ave and its variance as var.
      INTEGER j
      REAL s,ep
      ave=0.0
      m=min(1000,n)
      do 11 j=1,m
         ave=ave+data(j)
 11   continue
      ave=ave/m
      var=0.0
      ep=0.0
      do 12 j=1,n
         s=data(j)-ave
         ep=ep+s
         var=var+s*s
 12   continue
      var=(var-ep**2/n)/(n-1) !Corrected two-pass formula (14.1.8).
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE avevar2(data,data2,n,ave,var)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n,m,n2
      REAL ave,var,data(n),std,data2(n),i
C Given array data(1:n), returns its mean as ave and its variance as var.
      INTEGER j
      REAL s,ep
c      ave=0.0
c      m=min(1000,n)
c      do 11 j=1,m
c         ave=ave+data(j)
c 11   continue
c      ave=ave/m
      var=0.0
      ep=0.0
      do 12 j=1,n
         s=data(j)-ave
         ep=ep+s
         var=var+s*s
 12   continue
      var=(var-ep**2/n)/(n-1) !Corrected two-pass formula (14.1.8).
      std=sqrt(var)

      do 15 i=1,3

      n2=0
      do 13 j=1,n
         if(abs(data(j)-ave).lt.3.0*sqrt(var))then
            n2=n2+1
            data2(n2)=data(j)
         endif
 13   continue

      var=0.0
      ep=0.0
      do 14 j=1,n2
         s=data2(j)-ave
         ep=ep+s
         var=var+s*s
 14   continue
      var=(var-ep**2/n2)/(n2-1) !Corrected two-pass formula (14.1.8).

 15   continue

      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function comppop(nmax,ii,bx,ii2,bx2,bdx1,bdy1,bdx2,bdy2,
     .   nbin,datamin,datamax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer ii,ii2,nbin,nmax,i
      real bx(nmax),bx2(nmax),chi,bmax,bdx1(nmax),bdy1(nmax),bdx2(nmax),
     .   bdy2(nmax),datamin,datamax,sig

      call bindata(nbin,ii,bx,bdx1,bdy1,datamin,datamax,bmax)
      call bindata(nbin,ii2,bx2,bdx2,bdy2,datamin,datamax,bmax)

c      do 10 i=1,nbin
c         write(6,*) i,bdx1(i),bdx2(i),bdy1(i),bdy2(i)
c 10   continue

      chi=0.0
      do 11 i=1,nbin
         sig=sqrt(bdy1(i)+bdy2(i))
         if(sig.gt.0) then
            chi=chi+(bdy1(i)-bdy2(i))*(bdy1(i)-bdy2(i))/(sig*sig)
         endif
 11   continue
      comppop=chi

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bindata(nbin,npt,pdata,bdatax,bdatay,datamin,datamax,
     .   bmax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nbin,npt,i,bin
      real pdata(npt),bdatax(nbin),bdatay(nbin),datamin,datamax,
     .   binsize,bmax,tsum

C     get size of bin.
      binsize=(datamax-datamin)/real(nbin-1)

C     initalize bdata
      do 20 i=1,nbin
         bdatay(i)=0.
 20   continue

C     loop over data and place data into proper bin.
      do 10 i=1,npt
C        calculate bin number
         bin=int((pdata(i)-datamin)/binsize)+1
         if((bin.gt.0).and.(bin.le.nbin)) bdatay(bin)=bdatay(bin)+1.0
 10   continue

C     get the max value in a bin and assign bin value.
      tsum=0.
      bmax=0.
      do 30 i=1,nbin
         bmax=max(bdatay(i),bmax)
         bdatax(i)=datamin+real(i-1)*binsize!+binsize/2.0 !added -1
         tsum=tsum+bdatay(i)
 30   continue
c      write(6,*) "Sum:",tsum

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine makepop(nfit,ans,seed,nmax,ndata,data,n,mass,
     .   x,y,z,ii2,bx2,dumr,mnr,mrrel,nplot,nmp,mp,nbstar)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer ii2,i,j,k,nmax,imass,seed,mnr,n,ndata,nhavetransit,ieccn,
     .   nplot,nfit,np,nmp(6),jj,usebin(nmax),nbstar,nbflag,nbt
      real rhoavg,mass(nmax),minmass,kroupa,mbinary,ran2,mmax,m1,m2,
     .   rbinary,massradfn,mrrel(nmax,4),rhobin,Msun,Rsun,G,Pi,
     .   data(nmax,ndata),rbin,dumr2,rhostar,x(nmax),y(nmax),z(nmax),
     .   Per,adrs,rmax,schwarz,eccn,r1,r2,w,Eanom,trueanomaly,Tanom,
     .   ddrs,distance,TransitProb,rTransit,bx2(nmax),binfrac,eccnsig,
     .   dumr(nmax),lummass,ans(nfit),a,b,Pbeta,gasdev,dy,mp(nmax,6),
     .   escale,qmin,q,qpow,fratio,lumin,la,lb

      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      Pi=acos(-1.d0)!define Pi and 2*Pi
      Rsun=6.9599E8
      Msun=1.989E30

      nbstar=0
      nbt=0
      fratio=0.02

      binfrac=ans(1)
      a=ans(2)
      b=ans(3)
      escale=ans(4)
      qpow=ans(5)

      eccnsig=1.0
      if((a.le.0.0).and.(b.le.0.0)) eccnsig=0.0

      do 10 i=1,6
         nmp(i)=0
 10   continue


      ii2=0
      do 30 i=1,nmax
         nbflag=0 !count number of blends per star
         j=0 !reset counter
         rhoavg=0.0
C        Generate a potential binary-star companion
         if(mass(i).gt.0.0)then
            minmass=0.08  !replace with luminosity limit
            jj=0
            do 36 k=1,n
               if((i.eq.int(data(k,1))).and.(data(k,21).ge.7.1).and.
     .          (data(k,24).gt.1))then
                  jj=jj+1
                  rbin=ran2(seed)
                  if(rbin.gt.binfrac)then
                     usebin(jj)=0
                  else
                     usebin(jj)=1
                     minmass=max(minmass,
     .                lummass(nmax,mnr,mrrel,mass(i),data(i,14)))
                  endif
               endif
 36         continue

c            imass=0
c            mmax=kroupa(minmass)
c            do while(imass.eq.0)
c               mbinary=ran2(seed)*(mass(i)-minmass)+minmass
c               m1=mmax*ran2(seed)
c               m2=kroupa(mbinary)
c               if(m1.le.m2) imass=1
c            enddo

C           binary star mass ratio N(q)~q-1 where q=M2/M1
            imass=0
            mmax=(minmass/mass(i))**qpow!mass(i)/minmass
            qmin=minmass/mass(i)
            do while(imass.eq.0)
               q=ran2(seed)*(1.0-qmin)+qmin
               m1=mmax*ran2(seed)
               m2=q**qpow!1.0/q
               if(m1.le.m2) then
                  imass=1
                  mbinary=mass(i)*q
                endif
            enddo
c            write(0,*) qmin,q,mbinary
c            read(5,*)

c            mbinary=ran2(seed)*(mass(i)-minmass)+minmass
c               write(0,*) "Mbinary: ",mbinary,mmax,m1,m2
            rbinary=massradfn(nmax,mnr,mrrel,mbinary)
            rhobin=mbinary*Msun/(4.0/3.0*pi*(rbinary*Rsun)**3.0)/1000.0
         endif
c         write(0,*) "jj:",jj!(usebin(k),k=1,jj)
         do 31 k=1,n
            if((i.eq.int(data(k,1))).and.(data(k,21).ge.7.1).and.
     .        (data(k,24).gt.1))then
               j=j+1
               rbin=ran2(seed)
c               dumr2=data(k,12)
c               if(rbin.gt.binfrac)then
               if(usebin(j).eq.0)then
                  rhostar=data(k,12) !use stellar properties <rhostar>
c                  rhostar=rhoavgmeas(i)
               else
                  rhostar=rhobin !adapt binary star companion
                  la=lumin(nmax,mnr,mrrel,mass(i))
                  lb=lumin(nmax,mnr,mrrel,mbinary)
                  if((lb/la.gt.fratio).and.(nbflag.eq.0)) then
                     nbstar=nbstar+1
                     nbflag=1
c                     write(0,*) "nbstar: ",nbstar,lb,la
                  endif
               endif
               nbt=nbt+1
               if(rhostar.le.0) rhostar=1.0 !avoid non-sense
               dumr(j)=rhostar
c               write(0,*) "x:",data(k,12),rhobin,rhostar
c               x(j)=rhostar
               z(j)=data(k,11)/data(k,10) !use transit derived error
               Per=data(k,4)
               adrs=1000.0*rhostar*G*(Per*86400.0d0)**2/(3.0d0*Pi)
               adrs=adrs**(1.0d0/3.0d0)
               nhavetransit=0
               do while (nhavetransit.eq.0)
c                  rmax=schwarz(eccnsig,eccnsig)
                  rmax=Pbeta(0.001,a,b)
                  if(eccnsig.gt.0)then
                     ieccn=0
                     do while (ieccn.eq.0)
                        eccn=ran2(seed)
                        r1=rmax*ran2(seed)
c                        r2=schwarz(eccn,eccnsig)
                        r2=Pbeta(eccn,a,b)
                        if(r1.le.r2) ieccn=1
                     enddo
c                     write(0,*) eccn
                  else
                     eccn=0.0
                  endif
c                  eccn=ran2(seed)
                  w=2.0*pi*ran2(seed)
                  call kepler(w,Eanom,eccn)
                  Tanom=trueanomaly(eccn,Eanom)
                  ddrs=distance(adrs,eccn,Tanom)
                  TransitProb=1.0/ddrs  !probability of transit
                  rTransit=ran2(seed) !see if we have a transit
                  if(rTransit.le.TransitProb)then
                     nhavetransit=1
                  endif
               enddo
               y(j)=ddrs**3.0/(1000.0*G*(Per*86400.0)**2)*(3.0*Pi)
 34            dy=z(j)*y(j)*gasdev(seed)*escale
               if(y(j)+dy.gt.0.0)then
                  y(j)=y(j)+dy
               else
c                  write(0,*) i,y(j),dy
c                  read(5,*)
                  goto 34
               endif
c               rhoavg=rhoavg+y(j)
c               y(j)=(y(j)-rhostar)/z(j)
c               write(6,*) i,eccn,w
            endif
 31      continue
         np=j
c         rhoavg=rhoavg/real(j)
         if(np.gt.1)then
            do 33 j=1,np
               do 35 k=j+1,np
                  rhoavg=(y(j)+y(k))/2
                  ii2=ii2+1
                  bx2(ii2)=(y(j)-y(k))/(y(j)+y(k))
                  bx2(ii2)=abs(bx2(ii2))
                  if(nplot.eq.1)then
                     call pgsci(np)
                     call pgpt1(log10(rhoavg),bx2(ii2),4)
                  endif
                  nmp(np)=nmp(np)+1
                  mp(nmp(np),np)=bx2(ii2) !store histograms for 2,3,4..
c            write(6,500) dumr2,rhobin,dumr(k),rhoavg,y(k),z(k)
 500           format(6(F6.2,1X))
 35            continue
 33         continue
         endif
 30   continue
      if(nplot.eq.1) call pgsci(1)

c      write(0,*) nbt,nbstar,real(nbstar)/real(nbt)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function massradfn(nmax,mnr,mrrel,mass)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nmax,mnr,i
      real mrrel(nmax,4),mass,rad,G,Rsun,Msun,lg

      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      Rsun=6.9599E8
      Msun=1.989E30

      if(mass.lt.mrrel(1,1))then
         rad=sqrt(Msun*mrrel(1,1)*G/(10**(mrrel(1,3)-2.0)))/Rsun
      elseif(mass.gt.mrrel(mnr,1))then
         rad=sqrt(Msun*mrrel(mnr,1)*G/(10**(mrrel(mnr,3)-2.0)))/Rsun
      else
         do 10 i=1,mnr
            if((mass.ge.mrrel(i,1)).and.(mass.lt.mrrel(i+1,1)))then
               lg=(mass-mrrel(i,1))/(mrrel(i+1,1)-mrrel(i,1))*
     .          (mrrel(i+1,3)-mrrel(i,3))+mrrel(i,3)
               rad=sqrt(Msun*mass*G/(10**(lg-2.0)))/Rsun
            endif
 10      continue
      endif

c      write(0,*) "Rad:", rad,mass
      massradfn=rad

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function lumin(nmax,mnr,mrrel,mass)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer nmax,mnr,i
      real mrrel(nmax,4),mass,logl

      if(mass.lt.mrrel(1,1))then
         logl=mrrel(1,4)
      elseif(mass.gt.mrrel(mnr,1))then
         logl=mrrel(mnr,4)
      else
         do 10 i=1,mnr
            if((mass.ge.mrrel(i,1)).and.(mass.lt.mrrel(i+1,1)))then
               logl=(mass-mrrel(i,1))/(mrrel(i+1,1)-mrrel(i,1))*
     .          (mrrel(i+1,4)-mrrel(i,4))+mrrel(i,4)
            endif
 10      continue
      endif

      lumin=10.0**logl

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function lummass(nmax,mnr,mrrel,mass,tdepth)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nmax,mnr,i
      real mrrel(nmax,4),mass,rad,G,Rsun,Msun,lum,f1,f2,massb,
     .   tdepth

      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      Rsun=6.9599E8
      Msun=1.989E30

      if(mass.lt.mrrel(1,1))then
         f1=10**mrrel(1,4)
      elseif(mass.gt.mrrel(mnr,1))then
         f1=10**mrrel(mnr,4)
      else
         do 11 i=1,mnr
            if((mass.ge.mrrel(i,1)).and.(mass.lt.mrrel(i+1,1)))then
               f1=(mass-mrrel(i,1))/(mrrel(i+1,1)-mrrel(i,1))*
     .          (mrrel(i+1,4)-mrrel(i,4))+mrrel(i,4)
               f1=10**f1
            endif
 11      continue
      endif

      f2=f1*2.0*tdepth/1.0e6
c      write(0,*) f1,f2

      lum=log10(f2) !find lowest mass based on this luminosity

      if(lum.lt.mrrel(1,4))then
         massb=mrrel(1,1)
      elseif(lum.gt.mrrel(mnr,4))then
         massb=mrrel(mnr,1)
      else
         do 10 i=1,mnr
            if((lum.ge.mrrel(i,4)).and.(lum.lt.mrrel(i+1,4)))then
               massb=(lum-mrrel(i,4))/(mrrel(i+1,4)-mrrel(i,4))*
     .          (mrrel(i+1,1)-mrrel(i,1))+mrrel(i,1)
            endif
 10      continue
      endif

c      write(0,*) "minmass: ",f1,tdepth,10**lum,massb
c      read(5,*)
      lummass=massb

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function kroupa(mass)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      real mass,a1,a2,a3

      a1=0.3
      a2=1.3
      a3=2.3

      if(mass.lt.0.8)then
         kroupa=mass**(-1.0*a1)
      elseif((mass.ge.0.8).and.(mass.lt.0.5))then
         kroupa=mass**(-1.0*a2)
      else
         kroupa=mass**(-1.0*a3)
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function Pbeta(x,a,b)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      real x,a,b,gammln,beta

      beta=exp(gammln(a)+gammln(b)-gammln(a+b))

      Pbeta=(1.0/beta)*x**(a-1.0)*(1-x)**(b-1.0)
c      write(6,*) Pbeta,a,b,gammln(a+b)
c      read(5,*)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real FUNCTION gammln(xx)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
C     From Numerical Recipes
      real xx
C     Returns the value ln[Γ(xx)] for xx > 0.
      INTEGER j
      double precision ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     * 24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     * -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
         y=y+1.d0
         ser=ser+cof(j)/y
 11   continue
      gammln=tmp+log(stp*ser/x)
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function schwarz(x,sx)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      real x,sx

      schwarz=x/(sx*sx)*exp(-x*x/(2.0*sx*sx))

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real FUNCTION ran2(idum)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     * IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     * IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 11 j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
 11      continue
         iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real FUNCTION gasdev(idum)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER idum
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran2
      SAVE iset,gset
      DATA iset/0/
      if (idum.lt.0) iset=0
      if (iset.eq.0) then
 1       v1=2.0d0*ran2(idum)-1.0d0
         v2=2.0d0*ran2(idum)-1.0d0
         rsq=v1**2+v2**2
         if(rsq.ge.1.0d0.or.rsq.eq.0.0d0)goto 1
         fac=sqrt(-2.0d0*log(rsq)/rsq)
         gset=v1*fac
         gasdev=v2*fac
         iset=1
      else
         gasdev=gset
         iset=0
      endif
      return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function distance(asep,eccn,Tanom)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      real asep,eccn,Tanom

      distance=asep*(1.0d0-eccn*eccn)/(1+eccn*cos(Tanom))

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real function trueanomaly(eccn,Eanom)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      real eccn,Eanom,temp(2)

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
      real Manom,Eanom,Eold,eccn,diff,thres
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
      real Manom,Eanom,Mold,eccn,diff,thres
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ovrwrt (line, iwhich)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Taken from the DAOPhot code by Stetson.
C     This is a cool routine for writing lots of info to the screen
C     and keeping everything on the same line.
C
      character*(*) line
      character*79 output
      integer len
      if (iwhich .eq. 1) then
         write (0,1) line
    1    format (a)
      else if (iwhich .eq. 2) then
         if (len(line) .lt. 79) then
            output = ' '
            output = line
            write (0,2) output, char(13), char(13)
            write (0,2) output, char(13), char(13)
            write (0,2) output, char(13), char(13)
    2       format (a, 2a1, $)
         else
            write (0,2) line, char(13), char(13)
         end if
      else if (iwhich .eq. 3) then
         write (0,3) line
    3    format (a)
      else
         write (0,4) line, char(13), char(13)
    4    format (/a, 2a1, $)
         write (0,2) line, char(13), char(13)
      end if
      return
      end
