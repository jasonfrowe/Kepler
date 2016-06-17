      program rhoboot2
      implicit none
      integer iargc,nbin,nbinmax,nmax,nunit,np,i
      parameter(nmax=5000000,nbinmax=500)
      integer nd(nmax)
      real rp(nmax),bdatax(nbinmax),bdatay(nbinmax),errs(6),ave,std,
     .  output(16)
      double precision mstar(nmax),age(nmax),z(nmax),rstar(nmax),
     .  rhostar(nmax),temp(nmax),lum(nmax),work(nmax),dd(nmax)
      character*80 rhofile,cline,title,name

      if(iargc().lt.1) goto 901 !check number of input parameters
      
      nbin=30 !change to optional input parameter and check nbinmax
      if(iargc().ge.2)then
        call getarg(2,cline)
        read(cline,*) nbin
      endif
      if(nbin.lt.2) goto 905
      if(nbin.gt.nbinmax) goto 906
      
      call getarg(1,rhofile) !get mcmc parameters from rhostar analysis.
      nunit=10
      open(unit=nunit,file=rhofile,status='old',err=902)
      call getrhostar(nunit,np,nmax,mstar,age,z,rstar,rhostar,temp,lum)
      close(nunit)

c      call pgopen('/null') !open PGPlot device
      call pgopen('?')
      call pgask(.true.) !don't ask for new page.. just do it.
      call PGPAP ( 8.0 ,1.0) !paper size
      call pgsubp(4,4)  !break up plot into grid
      
      write(6,501) "Parameter    ","Median    ", "Stdev     ", 
     .  "+1 sig     ", "-1 sig     ", "+2 sig     ", "-2 sig     ", 
     .  "+3 sig     ", "-3 sig     "
      write(6,502) "----------------------------------------------------
     .------------------------------------------------------------------
     .---"
 501  format(A18,8(A13))   
 502  format(A122)  

C     Now lets make histograms of the stellar-parameters.
      title='Stellar Mass (M\d\(2281)\u)'
      call histogram(np,rp,mstar,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs)
c      write(6,500) ave,std,(errs(i),i=1,6)
      name="M* (Msun)"
      call writetable(ave,std,errs,name)
 500  format(9(F11.6,1X))
      output(1)=ave
      output(2)=std
        
      title='Stellar Radius (R\d\(2281)\u)'
      call histogram(np,rp,rstar,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs)
c      write(6,500) ave,std,(errs(i),i=1,6)
      name="R* (Rsun)"
      call writetable(ave,std,errs,name)
      output(3)=ave
      output(4)=std
     
      title='Z'
      call histogram(np,rp,z,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs)
c      write(6,500) ave,std,(errs(i),i=1,6)
      name="Z"
      call writetable(ave,std,errs,name)
      output(5)=ave
      output(6)=std
     
      title='T\deff\u (K)'
      call histogram(np,rp,temp,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs)
c      write(6,500) ave,std,(errs(i),i=1,6)
      name="Teff (K)"
      call writetable(ave,std,errs,name)
      output(7)=ave
      output(8)=std
            
      title='Stellar Luminosity (L\d\(2281)\u)'
      call histogram(np,rp,lum,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs)
c      write(6,500) ave,std,(errs(i),i=1,6)
      name="L (Lsun)"
      call writetable(ave,std,errs,name)
      output(9)=ave
      output(10)=std
      
      title='Age (Gyr)'
      call histogram(np,rp,age,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs)
c      write(6,500) ave,std,(errs(i),i=1,6)
      name="Age (Gyr)"
      call writetable(ave,std,errs,name)
      output(11)=ave
      output(12)=std
  
C     logg here
      call getlogg(np,dd,mstar,rstar)
      title='log(g)'
      call histogram(np,rp,dd,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs)
c      write(6,500) ave,std,(errs(i),i=1,6)
      name="log(g) (cgs)"
      call writetable(ave,std,errs,name)
      output(13)=ave
      output(14)=std

C     rhostar here
      title='Stellar Density (g/cm\u3\d)'
      call histogram(np,rp,rhostar,work,nd,nbin,nbinmax,bdatax,bdatay,
     .  title,ave,std,errs)
c      write(6,500) ave,std,(errs(i),i=1,6) 
      name="rhostar (g/cm^3)"
      call writetable(ave,std,errs,name)
      output(15)=ave
      output(16)=std
      
      call pgclos()

      write(6,503) (output(i),i=1,16)
 503  format(16(F9.4,1X))
c      write(6,504) int(output(7)),'|',int(output(8)),'|',output(5),'|',
c     .   output(6),'|',output(15),'|',output(16),'|',output(3),'|',
c     .   output(4),'|',output(1),'|',output(2)
 504  format(2(I4,A1),2(F5.3,A1),2(F8.6,A1),4(F5.3,A1))

      goto 999
 901  write(0,*) "Usage: rhoboot2 rhoboot.file <bins>"
      goto 999
 902  write(0,*) "Cannot open ",rhofile
      goto 999
 905  write(0,*) "Error: Nbin must be at least 2"
      goto 999
 906  write(0,*) "Increase nbinmax to at least",nbin
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
        if(ave.lt.10000.0)then
            write(6,507) name,ave,std,(errs(i),i=1,6)
 507        format(A18,8(F5.0,8X))
        else
            write(6,508) name,ave,std,(errs(i),i=1,6)
 508        format(A18,8(F7.0,6X))
        endif
      endif
      
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getrhostar(nunit,np,nmax,mstar,age,z,rstar,rhostar,
     .  temp,lum)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer np,nmax,nunit,i
      double precision mstar(nmax),age(nmax),z(nmax),rstar(nmax),
     .  rhostar(nmax),temp(nmax),lum(nmax)
     
      i=1
 10   read(nunit,*,end=11) mstar(i),age(i),z(i),rstar(i),rhostar(i),
     .  temp(i),lum(i)
        lum(i)=10.0**lum(i)
        rhostar(i)=rhostar(i)/1.0d3
        i=i+1
      goto 10
 11   continue
     
      np=i-1
c      write(0,*) "Stellar parameters read:",np
     
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getlogg(np,logg,mstar,rstar)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer np,i
      double precision logg(np),mstar(np),rstar(np),Msun,Rsun,G
      
      Msun=1.9891d30 !kg  mass of Sun
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      G=6.674d-11 !N m^2 kg^-2  Gravitation constant
      
      do 10 i=1,np
        logg(i)=log10(1.0e2*G*Msun*mstar(i)/
     .      (Rsun*Rsun*rstar(i)*rstar(i)))
  10  continue
      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine histogram(npt,rp,dp,dp2,np,nbin,nbinmax,bdatax,bdatay,
     .  title,rmed,std,errs)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,i,j,npt2,nbin,nbinmax
      integer np(npt)
      real rp(npt),datamin,datamax,bdatax(nbinmax),bdatay(nbinmax),bmax,
     .  rave,errs(6),std,rmed,bsize,mode,momax
      double precision dp(npt),dp2(npt),ave,var,med
      character*80 title
      
C     First we remove the average (helps a lot with double -> real)
      call avevar(dp,npt,ave,var)
      rave=real(ave) !convert to real*4
      std=real(sqrt(var))
      
C     Now we convert dp to rp
      j=0 !because we have sigma-clipping, we need a counter.
      do 10 i=1,npt
        if(abs(dp(i)-ave).lt.4.0*sqrt(var))then
            j=j+1
            dp2(j)=dp(i)-ave
            rp(j)=real(dp2(j))
        endif
 10   continue
      npt2=j

C     Find median          
      call rqsort(npt2,dp2,np) !changed npt to k
      i=npt2/2
      if(i.le.0) i=1
      med=dp2(np(i))
      rmed=real(med)

C     Find datarange
      datamin=rp(1)
      datamax=rp(1)
      do 12 i=2,npt2
        datamin=min(rp(i),datamin)
        datamax=max(rp(i),datamax)
 12   continue
 
      call bindata(nbin,npt2,rp,bdatax,bdatay,datamin,datamax,bmax)

      bsize=(bdatax(2)-bdatax(1))/2.0
C     find mode
      mode=bdatax(1)+bsize
      momax=bdatay(1)
      do 27 j=2,nbin
        if(bdatay(j).gt.momax)then
            momax=bdatay(j)
            mode=bdatax(j)+bsize
        endif
 27   continue 
c      rmed=mode
      
      call pgpage() !fresh plotting surface
      call pgslw(1)
      call pgsch(2.0)
c         call windowsetup(xb1,xb2,yb1,yb2) !make a square plotting surface
c         call pgvport(xb1,xb2,yb1,yb2)
      call pgvport(0.2,1.0,0.2,0.9)
c      fc=2.*(datamax-datamin)/real(nbin) !center histograms
      call pgwindow(datamin,datamax,0.,bmax+0.1*bmax) !set size
         
C        Add axis labels
      call pglabel(title,"Relative Probability","")
      call pgbin(nbin,bdatax,bdatay,.false.) !plot the histogram
      
C     Shift axis scale to account for average removal
      call pgwindow(datamin+rave,datamax+rave,0.,1.0+0.1*1.0)
      call pgbox('BCNTS1',0.0,0,'BCNTS1',0.0,0) !add boarders
      
      call errorest(npt2,rp,nbin,bdatax,bdatay,bmax,rmed,errs)
c      call errorest2(nbin,bdatax,bdatay,bmax,rmed,errs)
      rmed=rmed+rave !correct for average removal
C     Need to recalulate standard deviation at this point!      

      return
      end      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine errorest2(nbin,bdatax,bdatay,bmax,frsol,errs)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nbin,i,j,k,nbinmax,mdpnt
      parameter(nbinmax=500)
      integer bsy(nbinmax)
      real bdatax(nbin),bdatay(nbin),bmax,frsol,errs(6),sumby,shgt,
     .  inth,intl,sumval,per,perold,inthold,intlold,e1,e2,bx
      logical loop
      
      sumby=0
      do 6 i=1,nbin
        sumby=sumby+bdatay(i)
 6    continue
c      write(0,*) "sumby:",sumby
      
      call rqsortr(nbin,bdatay,bsy)
      bx=bdatax(bsy(nbin))
      
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
c      write(0,*) "mdpnt:",mdpnt,bdatax(mdpnt)

      perold=0 !initalization of variables
      intlold=0
      inthold=0
      do 14 i=1,6  !initialize variable to zero.
        errs(i)=0.0
 14   continue
      do 13 k=nbin,1,-1
        shgt=real(bdatay(bsy(k)))
c        write(0,*) shgt

        loop=.true. !loop until loop is false
        i=mdpnt  !removed -1
        do 10 while(loop)
            if((bdatay(i).ge.shgt).and.(bdatay(i+1).lt.shgt))then
                loop=.false.
                inth=real(i)!bdatax(i)
c                inth=(bdatay(i+1)-shgt)/(bdatay(i+1)-bdatay(i))
c     .              *(bdatax(i+1)-bdatax(i))+
c     .              bdatax(i)
            endif
            if(i.ge.nbin-1) then
                loop=.false.
                inth=real(nbin)!bdatax(nbin)
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
                intl=real(i)!bdatax(i)
c                intl=bdatax(i)-
c     .              (bdatay(i)-shgt)/(bdatay(i)-bdatay(i-1))
c     .              *(bdatax(i)-bdatax(i-1))
            endif
            if(i.le.2) then
                loop=.false.
                intl=real(1)!bdatax(1)
            endif
            i=i-1
 11     enddo

        sumval=0.0
        do 15 i=int(intl),int(inth)
            sumval=sumval+bdatay(i)
 15     continue
        per=sumval/sumby
        
        if((per.ge.0.683).and.(perold.lt.0.683))then
            e1=bdatax(inth)-(per-0.683)/(per-perold)*
     .          (bdatax(inth)-bdatax(inthold))
            e2=bdatax(intl)-(per-0.683)/(per-perold)*
     .          (bdatax(intl)-bdatax(intlold))
            errs(1)=e1-bx!-frsol
            errs(2)=e2-bx!-frsol
        endif
        if((per.ge.0.954).and.(perold.lt.0.954))then
            e1=bdatax(inth)-(per-0.954)/(per-perold)*
     .          (bdatax(inth)-bdatax(inthold))
            e2=bdatax(intl)-(per-0.954)/(per-perold)*
     .          (bdatax(intl)-bdatax(intlold))
            errs(3)=e1-bx!-frsol
            errs(4)=e2-bx!-frsol
        endif
        if((per.ge.0.9973).and.(perold.lt.0.9973))then
            e1=bdatax(inth)-(per-0.9973)/(per-perold)*
     .          (bdatax(inth)-bdatax(inthold))
            e2=bdatax(intl)-(per-0.9973)/(per-perold)*
     .          (bdatax(intl)-bdatax(intlold))
            errs(5)=e1-bx!-frsol
            errs(6)=e2-bx!-frsol
        endif
        
        perold=per
        intlold=intl
        inthold=inth

c      write(0,*) shgt,intl,inth,per
c      write(0,*) e1,e2
        
 13   continue
      
c      write(0,*) frsol
c      read(5,*)      
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine errorest(npt,dd,nbin,bdatax,bdatay,bmax,frsol,errs)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer npt,nbin,i,j,mdpnt,k,nbinmax
      parameter(nbinmax=500)
      integer bsy(nbinmax)
      real dd(npt),bdatax(nbin),bdatay(nbin),bmax,frsol,errs(6),shgt,
     .  intl,inth,per,perold,e1,e2,inthold,intlold
      logical loop
      
      call rqsortr(nbin,bdatay,bsy)
      
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

c**********************************************************************
      subroutine rqsortr(n,a,p)
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
      SUBROUTINE avevar(data,n,ave,var)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      INTEGER n,m
      REAL*8 ave,var,data(n)
C Given array data(1:n), returns its mean as ave and its variance as var.
      INTEGER j
      REAL*8 s,ep
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
      real*8      a(n)

c     Output:
      integer   p(n)

c     Constants
      integer   LGN, Q
      parameter (LGN=32, Q=11)
c        (LGN = log base 2 of maximum n;
c         Q = smallest subfile to use quicksort on)

c     Local:
      real*8      x
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
