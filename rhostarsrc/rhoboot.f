      program mcmchisto
      implicit none
      integer nmax,npt,col,nbin,i,j,nfit,nunit,nbinmax,k,npt2
      parameter(nmax=1000000,nfit=7,nbinmax=500)
      integer nd(nmax)
      real dd(nmax),datamin,datamax,bdatax(nmax),
     .   bdatay(nmax),bmax,fc,ave,var,pi,ans(2,nmax),
     .   errs(nfit),sol(nfit),serr(nfit,2),err(nfit),pmin,pmax
      real G,Msun,Mearth,Mjup,
     .  Rsun,aConst,out1(nbinmax,nfit),AU,
     .  out2(nbinmax,nfit),bsize,Rjup,med(nfit),mode(nfit),
     .  momax,pid2,tpi,outradmass(4)
      real per(nmax),bb(nmax),zdr(nmax),ecw(nmax),esw(nmax),eccn,w,adrs
      character*3 titles2(nfit)
      character*80 filename,inputsol,cline
      character*80 titles(nfit),name,names(nfit)
      titles(1)='Stellar Mass (M\d\(2281)\u)'
      titles(2)='Age (Gyr)'
      titles(3)='Z'
      titles(4)="Stellar Radius (R\d\(2281)\u)"
      titles(5)="Mean Stellar Density (g/cm\u3\d)"
      titles(6)='Temperature (K)'
      titles(7)='Stellar Luminosity (L\d\(2281)\u)'
      
      names(1)='M* (Msun)'
      names(2)='Age (Gyr)'
      names(3)='Z'
      names(4)='R* (Rsun)'
      names(5)='rhostar (g/cm^3)'
      names(6)='Teff (K)'
      names(7)='L* (Lsun)'
            
      Pi=acos(-1.d0)   !Pi
      pid2=Pi/2.0 !Pi/2
      tPi=2.0*Pi !2Pi
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      Msun=1.9891d30 !kg  mass of Sun
      Mearth=5.974d24 !kg mass of Earth
      Mjup=317.833d0*Mearth !kg  mass of Jupiter  
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      Rjup=142980.0d0*1000.0d0/2.0d0 !m  radius of Jupiter 
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)
      AU=1.49598e11
      
      nbin=30 !default binning for histograms
      if(iargc().lt.1) goto 903 !need one commandline arguement
           
      call getarg(1,filename)
      open(unit=10,file=filename,status='old',err=901)
      if(iargc().ge.2)then
        call getarg(2,cline)
        read(cline,*) nbin
      endif
      if(nbin.lt.2) goto 902

      if(nbin.gt.nbinmax) goto 904
            
      call pgopen('?') !open PGPlot device
c      call pgopen('/null')
      call pgask(.false.) !don't ask for new page.. just do it.
      call PGPAP ( 6.0 ,1.0)
      call pgsubp(4,4)
c      call pgpage() !start a multi-page document

      write(6,503) "Parameter    ","Median    ", "Stdev     ", 
     .  "+1 sig     ", "-1 sig     ", "+2 sig     ", "-2 sig     ", 
     .  "+3 sig     ", "-3 sig     "
      write(6,504) "----------------------------------------------------
     .------------------------------------------------------------------
     .---"
 503  format(A18,8(A13))   
 504  format(A122)  

C     Now we loop over all fitted parameters 
      do 10 i=1,nfit
c         if(i.eq.2) goto 10
c         if(i.eq.5) goto 10
         col=i
         npt=nmax
         call getdata(col,npt,dd) !read column of bootstrapped values
         npt2=npt

C     These lines will sigclip at 4 sigma if required.  (probably not)
         call avevar(dd,npt,ave,var)
         k=0
         do 12 j=1,npt
            if(abs(dd(j)-ave).lt.4.0*sqrt(var))then
               k=k+1
               dd(k)=dd(j)
            endif
 12      continue
         npt=k
C     End of sigmaclipping

         call avevar(dd,npt,ave,var)
C     If we have a phase parameter, then we need to check for cyclic 
C     numbers around Pi

         ans(1,i)=ave
         ans(2,i)=sqrt(var)

C        Find median          
         call rqsort(npt,dd,nd) !changed npt to k
         med(i)=dd(nd(npt/2))
      
         call findrange(npt,dd,datamin,datamax)
         if(datamin.eq.datamax)goto 11
         
c         pmin=dd(1)
c         pmax=dd(2)
c         do 28 j=2,npt
c            pmin=min(dd(j),pmin)
c            pmax=max(dd(j),pmax)
c 28      continue
c         write(6,500) (pmax+pmin)/2.0,(pmax-pmin)/2.0
         
         
c         write(6,*) "datamin,datamax",datamin,datamax
         call bindata(nbin,npt,dd,bdatax,bdatay,datamin,datamax,bmax)

         bsize=(bdatax(2)-bdatax(1))/2.0
         do 19 j=1,nbin
            out1(j,i)=bdatax(j)+bsize
            out2(j,i)=bdatay(j)/bmax
 19      continue
 
C        find mode
         mode(i)=bdatax(1)+bsize
         momax=bdatay(1)
         do 27 j=2,nbin
            if(bdatay(j).gt.momax)then
                momax=bdatay(j)
                mode(i)=bdatax(j)+bsize
            endif
 27      continue           
      
         call pgpage() !fresh plotting surface
         call pgslw(1)
         call pgsch(2.0)
c         call windowsetup(xb1,xb2,yb1,yb2) !make a square plotting surface
c         call pgvport(xb1,xb2,yb1,yb2)
         call pgvport(0.2,1.0,0.2,0.9)
         fc=2.*(datamax-datamin)/real(nbin) !center histograms
         call pgwindow(datamin,datamax,0.,bmax+0.1*bmax) !set size
         
C        Add axis labels
         call pglabel(titles(i),"Relative Probability","")
         call pgbin(nbin,bdatax,bdatay,.false.) !plot the histogram
         
         
c         if(i.eq.6)then
c            call pgwindow(datamin*1.0e6,datamax*1.0e6,0.,1.1)
c         elseif(i.eq.7)then
c            call pgwindow(datamin-60.0,datamax-60.0,0.,1.1)
c         else
            call pgwindow(datamin,datamax,0.,1.1) !set size
c         endif
         call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0) !add boarders
          
c         if(i.eq.6) sol(i)=ans(1,i)
c         call errorest(npt,dd,nbin,bdatax,bdatay,bmax,mode(i),errs)
         call errorest2(npt,dd,nbin,bdatax,bdatay,bmax,med(i),errs)
         
         serr(i,1)=med(i)
         serr(i,2)=ans(2,i)

c         write(6,500) med(i),ans(2,i),(errs(j),j=1,6)
         name=names(i)
         call writetable(med(i),ans(2,i),errs,name)
         if(i.eq.1)then
            outradmass(1)=med(i)
            outradmass(2)=ans(2,i)
         elseif(i.eq.4)then
            outradmass(3)=med(i)
            outradmass(4)=ans(2,i)
         endif

c 500     format(9(1PE16.9,1X))
c 500     format(9(F11.7,1X))
 500     format(9(F12.6,1X))
c 500     format(9(F15.12,1X))
         
C        this line is heavy on disk-io but easy on memory.
 11      rewind(10) ! rewind the bootstrap file so we can read it again.
c         read(10,*) dumi !not needed, but no harm. Just skips first line
 10   continue
      call pgslw(1)
      call pgclos() !close plotting device

      write(6,505) (outradmass(i),i=1,4)
 505  format(4(F6.4,1X))
      
c      call exportfit(nfit,sol,serr,err,titles2)
      
c      do 12 i=0,nfreq-1
c         write(6,500) (ans(1,j),j=3*i+1,3*i+3),(ans(2,j),j=3*i+1,3*i+3)
c 12   continue
c 500  format(200(E14.7,1X))

c      open(unit=13,file="probs.dat")
c      write(13,502) "M*","M*_P","Mp","Mp_P","R*","R*_P","Rp","Rp_P",
c     .  "Per","Per_P","i","i_P","E","E_p","Ecl","Ecl_P","a/R*","a/R*_P"
c 502  format(18(A9,3X))
c      do 20 j=1,nbin
c     write(13,501) out1(j,1),out2(j,1),out1(j,2),out2(j,2),
c     .    out1(j,3),out2(j,3),out1(j,4),out2(j,4),out1(j,5),out2(j,5),
c     .    out1(j,6),out2(j,6),out1(j,7),out2(j,7),out1(j,16),out2(j,16),
c     .    out1(j,18),out2(j,18)
cc        write(6,501) out1(j,1),out2(j,1)
c 20   continue
c      close(13)
c 501  format(18(F11.7,1X))
      
      
      close(10)
      goto 999
 901  write(6,*) "Cannot open ",filename
      goto 999
 902  write(6,*) "Error: Nbin must be at least 2"
      goto 999
 903  write(6,*) "Usage: boothist <bootstats> [nbin]"
      write(6,*) "<bootstats>: Bootstrap statistics from rhostar"
      write(6,*) "[nbin] (optional) : binned parameter for histrograms, 
     .default=30"
      goto 999
 904  write(6,*) "Increase nbinmax to at least",nbin
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine writetable(ave,std,errs,name)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer i
      real ave,std,errs(6),minerr
      character*80 name
      
c      digit=log10(abs(errs(1)))
c      digit=min(log10(abs(errs(2))),digit)
c      ndigit=int(digit)
c      ndigit=ndigit-2
c      if(ndigit.lt.1) ndigit=ndigit-2
c      if(ndigit.ge.1) ndigit=ndigit-1
c      write(0,*) errs(1),errs(2)
c      write(0,*) ndigit,digit,log10(abs(errs(1))),log10(abs(errs(2)))
c      write(0,*) "ndigit:",ndigit

      minerr=min(errs(1),-errs(2))
      
      if(minerr.le.1.0e-5)then
        write(6,500) name,ave,std,(errs(i),i=1,6)
 500    format(A18,8(F12.7,1X))
      elseif(minerr.le.1.0e-4)then
        write(6,501) name,ave,std,(errs(i),i=1,6)
 501    format(A18,8(F11.6,2X))
      elseif(minerr.le.1.0e-3)then
        write(6,502) name,ave,std,(errs(i),i=1,6)
 502    format(A18,8(F10.5,3X))
      elseif(minerr.le.1.0e-2)then
        write(6,503) name,ave,std,(errs(i),i=1,6)
 503    format(A18,8(F9.4,4X))
      elseif(minerr.le.1.0e-1)then
        write(6,504) name,ave,std,(errs(i),i=1,6)
 504    format(A18,8(F8.3,5X))
      elseif(minerr.le.1.0)then
        write(6,505) name,ave,std,(errs(i),i=1,6)
 505    format(A18,8(F7.2,6X))
      elseif(minerr.le.10.0)then
        write(6,506) name,ave,std,(errs(i),i=1,6)
 506    format(A18,8(F6.1,7X))
      else
        write(6,507) name,ave,std,(errs(i),i=1,6)
 507    format(A18,8(F5.0,8X))
      endif
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine errorest2(npt,dd,nbin,bdatax,bdatay,bmax,frsol,errs)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nbin,i,j,mdpnt,k,nbinmax
      parameter(nbinmax=500)
      integer bsy(nbinmax)
      real dd(npt),bdatax(nbin),bdatay(nbin),bmax,frsol,errs(6),shgt,
     .  intl,inth,per,perold,e1,e2,inthold,intlold
      logical loop
      
      call rqsort(nbin,bdatay,bsy)
      
C     find out which bin contains the middle.      
      do 5 i=1,nbin-1
c        write(0,*) frsol,bdatax(i),bdatax(i+1)
        if((bdatax(i).le.frsol).and.(bdatax(i+1).gt.frsol))then
            mdpnt=i
        endif
 5    continue
C     case when "bestfit" is outside histogram (i.e. eccentricity)
      if(frsol.le.bdatax(1)) mdpnt=1
      if(frsol.ge.bdatax(nbin)) mdpnt=nbin 
c      write(0,*) "mdpnt:",mdpnt
      
      
      perold=0 !initalization of variables
      intlold=0
      inthold=0
      do 14 i=1,6  !initialize variable to zero.
        errs(i)=0.0
 14   continue
C     bmax is the value of the largest histogram bin
c      do 13 k=1,bmax-1
c        shgt=real(bmax-k)
      do 13 k=nbin,1,-1
        shgt=real(bdatay(bsy(k)))
        
C       first we move in the forward direction      
        loop=.true. !loop until loop is false
        i=mdpnt  !removed -1
        do 10 while(loop)
            if((bdatay(i).ge.shgt).and.(bdatay(i+1).lt.shgt))then
                loop=.false.
c               inth=bdatax(i)
                inth=(bdatay(i+1)-shgt)/(bdatay(i+1)-bdatay(i))
     .              *(bdatax(i+1)-bdatax(i))+
     .              bdatax(i)
            endif
            if(i.ge.nbin-1) then
                loop=.false.
                inth=bdatax(nbin)
            endif
            i=i+1
c           write(6,*) "i:",i
 10     enddo
       
C       now we move in the reverse direction
        loop=.true. !loop until loop is false
        i=mdpnt+1
        do 11 while(loop)
            if((bdatay(i).ge.shgt).and.(bdatay(i-1).lt.shgt))then
                loop=.false.
c               intl=bdatax(i)
                intl=bdatax(i)-
     .              (bdatay(i)-shgt)/(bdatay(i)-bdatay(i-1))
     .              *(bdatax(i)-bdatax(i-1))
            endif
            if(i.le.2) then
                loop=.false.
                intl=bdatax(1)
            endif
            i=i-1
 11     enddo

        j=0
        do 12 i=1,npt
            if((dd(i).gt.intl).and.(dd(i).lt.inth))then
                j=j+1    
            endif
 12     continue
 
        per=real(j)/real(npt)
        if((per.ge.0.683).and.(perold.lt.0.683))then
            e1=inth-(per-0.683)/(per-perold)*(inth-inthold)
            e2=intl+(per-0.683)/(per-perold)*(intlold-intl)
            errs(1)=e1-frsol
            errs(2)=e2-frsol
c           write(6,*) "1 sigma:",frsol,errs(1),errs(2)
        endif
        if((per.ge.0.954).and.(perold.lt.0.954))then
            e1=inth-(per-0.954)/(per-perold)*(inth-inthold)
            e2=intl+(per-0.954)/(per-perold)*(intlold-intl)
            errs(3)=e1-frsol
            errs(4)=e2-frsol
c           write(6,*) "2 sigma:",frsol,errs(3),errs(4)
        endif 
        if((per.ge.0.9973).and.(perold.lt.0.9973))then
            e1=inth-(per-0.9973)/(per-perold)*(inth-inthold)
            e2=intl+(per-0.9973)/(per-perold)*(intlold-intl)
            errs(5)=e1-frsol
            errs(6)=e2-frsol
c           write(6,*) "3 sigma:",frsol,errs(5),errs(6)
        endif       
        perold=per
        intlold=intl
        inthold=inth
c       write(6,*) per,intl,inth,shgt
 13   continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine errorest(npt,dd,nbin,bdatax,bdatay,bmax,frsol,errs)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nbin,i,j,mdpnt,k
      real dd(npt),bdatax(nbin),bdatay(nbin),bmax,frsol,errs(6),shgt,
     .  intl,inth,per,perold,e1,e2,inthold,intlold
      logical loop
      
C     find out which bin contains the middle.      
      do 5 i=1,nbin-1
c        write(0,*) frsol,bdatax(i),bdatax(i+1)
        if((bdatax(i).le.frsol).and.(bdatax(i+1).gt.frsol))then
            mdpnt=i
        endif
 5    continue
C     case when "bestfit" is outside histogram (i.e. eccentricity)
      if(frsol.le.bdatax(1)) mdpnt=1
      if(frsol.ge.bdatax(nbin)) mdpnt=nbin 
c      write(0,*) "mdpnt:",mdpnt
      
      perold=0 !initalization of variables
      intlold=0
      inthold=0
      do 14 i=1,6  !initialize variable to zero.
        errs(i)=0.0
 14   continue
C     bmax is the value of the largest histogram bin
      do 13 k=1,bmax-1
        shgt=real(bmax-k)
C       first we move in the forward direction      
        loop=.true. !loop until loop is false
        i=mdpnt  !removed -1
        do 10 while(loop)
            if((bdatay(i).ge.shgt).and.(bdatay(i+1).lt.shgt))then
                loop=.false.
c               inth=bdatax(i)
                inth=(bdatay(i+1)-shgt)/(bdatay(i+1)-bdatay(i))
     .              *(bdatax(i+1)-bdatax(i))+
     .              bdatax(i)
            endif
            if(i.ge.nbin-1) then
                loop=.false.
                inth=bdatax(nbin)
            endif
            i=i+1
c           write(6,*) "i:",i
 10     enddo
       
C       now we move in the reverse direction
        loop=.true. !loop until loop is false
        i=mdpnt+1
        do 11 while(loop)
            if((bdatay(i).ge.shgt).and.(bdatay(i-1).lt.shgt))then
                loop=.false.
c               intl=bdatax(i)
                intl=bdatax(i)-
     .              (bdatay(i)-shgt)/(bdatay(i)-bdatay(i-1))
     .              *(bdatax(i)-bdatax(i-1))
            endif
            if(i.le.2) then
                loop=.false.
                intl=bdatax(1)
            endif
            i=i-1
 11     enddo

        j=0
        do 12 i=1,npt
            if((dd(i).gt.intl).and.(dd(i).lt.inth))then
                j=j+1    
            endif
 12     continue
 
        per=real(j)/real(npt)
        if((per.ge.0.683).and.(perold.lt.0.683))then
            e1=inth-(per-0.683)/(per-perold)*(inth-inthold)
            e2=intl+(per-0.683)/(per-perold)*(intlold-intl)
            errs(1)=e1-frsol
            errs(2)=e2-frsol
c           write(6,*) "1 sigma:",frsol,errs(1),errs(2)
        endif
        if((per.ge.0.954).and.(perold.lt.0.954))then
            e1=inth-(per-0.954)/(per-perold)*(inth-inthold)
            e2=intl+(per-0.954)/(per-perold)*(intlold-intl)
            errs(3)=e1-frsol
            errs(4)=e2-frsol
c           write(6,*) "2 sigma:",frsol,errs(3),errs(4)
        endif 
        if((per.ge.0.9973).and.(perold.lt.0.9973))then
            e1=inth-(per-0.9973)/(per-perold)*(inth-inthold)
            e2=intl+(per-0.9973)/(per-perold)*(intlold-intl)
            errs(5)=e1-frsol
            errs(6)=e2-frsol
c           write(6,*) "3 sigma:",frsol,errs(5),errs(6)
        endif       
        perold=per
        intlold=intl
        inthold=inth
c       write(6,*) per,intl,inth,shgt
 13   continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine readfreqs(nfreq,nmax,frsol)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfreq,nmax,nvar,i,j
      real frsol(nmax)
      double precision ztime
      
      read(11,*) ztime
      read(11,*) nvar,nfreq
      
      do 10 i=1,nfreq
        j=3*(i-1)+1 !read all variables in order.
        read(11,*) frsol(j),frsol(j+1),frsol(j+2)
 10   continue
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getfitpars(nunit,nfit,sol,serr,err)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nfit,nunit,i
      real sol(nfit),serr(nfit,2),err(nfit)
      character*160 command
      
c      write(0,*) "----------------------------"
c      write(0,*) "Reading in starting solution"
c      write(0,*) "----------------------------"
      
      i=0 !initialize line counter
C     Start of loop to read in file
  10  read(nunit,500,end=11) command
 500  format(A160) !command much be contained within first 80 characters
        i=i+1 !increase line counter
        if(command(1:1).eq."#") then
            write(0,*) "Comment line on line ",i
        elseif(command(1:3).eq."EPO") then
            read(command(5:160),*) sol(1),serr(1,1),serr(1,2),err(1)
        elseif(command(1:3).eq."PER") then
            read(command(5:160),*) sol(2),serr(2,1),serr(2,2),err(2)
        elseif(command(1:3).eq."BBB") then
            read(command(5:160),*) sol(3),serr(3,1),serr(3,2),err(3)
        elseif(command(1:3).eq."RDR") then
            read(command(5:160),*) sol(4),serr(4,1),serr(4,2),err(4)
        elseif(command(1:3).eq."ZDR") then
            read(command(5:160),*) sol(5),serr(5,1),serr(5,2),err(5)
        elseif(command(1:3).eq."ZPT") then
            read(command(5:160),*) sol(6),serr(6,1),serr(6,2),err(6)
        elseif(command(1:3).eq."ECW") then
            read(command(5:160),*) sol(7),serr(7,1),serr(7,2),err(7)
        elseif(command(1:3).eq."ESW") then
            read(command(5:160),*) sol(8),serr(8,1),serr(8,2),err(8)
        elseif(command(1:3).eq."KRV") then
            read(command(5:160),*) sol(9),serr(9,1),serr(9,2),err(9)
        elseif(command(1:3).eq."VOF") then
            read(command(5:160),*) sol(10),serr(10,1),serr(10,2),err(10)
        elseif(command(1:3).eq."NL1") then
            read(command(5:160),*) sol(11),serr(11,1),serr(11,2),err(11)
        elseif(command(1:3).eq."NL2") then
            read(command(5:160),*) sol(12),serr(12,1),serr(12,2),err(12)
        elseif(command(1:3).eq."NL3") then
            read(command(5:160),*) sol(13),serr(13,1),serr(13,2),err(13)
        elseif(command(1:3).eq."NL4") then
            read(command(5:160),*) sol(14),serr(14,1),serr(14,2),err(14)
        endif
 501    format(A5,5(1X,1PE17.10))
      goto 10 !loop back to read statement
 11   continue !end loop when EOF is reached
c      write(0,*) "----------------------------"
             
      return
      end
 
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
         bdatax(i)=datamin+real(i-1)*binsize !added -1
         tsum=tsum+bdatay(i)
 30   continue
c      write(6,*) "Sum:",tsum
 
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine findrange(npt,pdata,datamin,datamax)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i
      real pdata(npt),datamin,datamax
      
C     initialize datamin and datamax
      datamin=pdata(1)
      datamax=pdata(1)
      
      do 10 i=2,npt
         datamin=min(pdata(i),datamin)
         datamax=max(pdata(i),datamax)
 10   continue
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine windowsetup(xb1,xb2,yb1,yb2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      setup plotting area
      real xb1,xb2,yb1,yb2,xs1,xs2,
     .   ys1,ys2,pratio
      
C     set text size      
      call pgsch(1.0)
C     inquire window size
      call pgqvsz(2,xs1,xs2,ys1,ys2)
      pratio=(xs2-xs1)/(ys2-ys1)
C     if we have a landscape surface
      if(xs2-xs1.gt.ys2-ys1)then
C        start 10% from the edge
         yb1=(ys2-ys1)*0.10+ys1
C        use 60% of surface available
         yb2=(ys2-ys1)*0.70+yb1
C        now normalize to 1
         yb1=yb1/(ys2-ys1)         
         yb2=yb2/(ys2-ys1)
C        now use same size in other dimension      
         xb1=yb1
         xb2=yb1+(yb2-yb1)/pratio
      else
C        start 10% from the edge
         xb1=(xs2-xs1)*0.10+xs1
C        use 60% of surface available
         xb2=(xs2-xs1)*0.60+xb1
C        now normalize to 1
         xb1=xb1/(xs2-xs1)         
         xb2=xb2/(xs2-xs1)
C        now use same size in other dimension      
         yb1=xb1
         yb2=xb1+(xb2-xb1)*pratio
      endif
      
      return
      end 
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getdata(col,npt,dd)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer col,npt,i,j,k
      real dd(npt),dumr
      
      k=0
      i=1
  9   continue
 10   read(10,*,end=11) (dumr,j=1,col-1),dd(i)
c          if(col.eq.6) dd(i)=log10(dd(i))
c         write(6,*) dd(i)
c         if((col.eq.6).and.(dd(i).gt.89.95)) goto 10
c         if(col.eq.8) dd(i)=dd(i)*1.0e6
          if((col.eq.2).and.(dd(i).gt.14.0))goto 10
          if(col.eq.7) dd(i)=10**dd(i)
          if(col.eq.5) dd(i)=dd(i)/1.0d3

c         k=k+1
c         if(k.lt.10) goto 10
c         k=0

         i=i+1
         goto 10
 11   continue
      
      npt=i-1
      
 500  format(200(E14.7,1X))
      goto 999
 901  write(6,*) "Error on line: ",i
      goto 999
      
 999  return
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

