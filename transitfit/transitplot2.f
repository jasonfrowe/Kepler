      program transitplot2
      implicit none
      integer iargc,nunit,ntype,kID,Teff,bins,nmax,nfit,i,j,ci1,
     .  ci2,icol,npt,inl(3),ncolour,kIDold,flag,nc,nplanet,np,nsc
      parameter(nmax=650000,nfit=18,ncolour=37)
      integer colourtemp(ncolour),sR(ncolour),sG(ncolour),sB(ncolour)
      real bounds(4),srad,sscale,r,nl(4),mu,Intensity,rp,In(3),br,x,y,
     .  plotsize,ymax,sradold,shiftx,shifty,ymin,th,dth
      double precision Kmag,logg,rad,eoff,escale,xwidth,sol(nfit),
     .  serr(nfit,2),Dpvary(nfit),err(nfit),doe,toff,Rsun,Rjup,Pi,G,M1,
     .  M2,Psec,asemi,aConst,Msun,Mearth,Mjup,time(nmax),mag(nmax),
     .  merr(nmax),itime(nmax),Zerotime,koi,T0,dumr,Per,adrs,rprs,tdur,
     .  Rstar,Teq,Tdepth,b,dx,dy
      double precision nl1(76,11,8),nl2(76,11,8),nl3(76,11,8),
     .  nl4(76,11,8)
      character*80 parsfile,filename,fitparsfile,dumc
      
      open(unit=10,
     .  file="/Users/rowe/Documents/transitfit/starcoloursRGB.txt",
     .  status='old',err=905)
      do 15 i=1,ncolour
        read(10,*) dumc,colourtemp(i),sR(i),sG(i),sB(i)
 15   continue
      close(10)
      
      call limbdarkensetup(nl1,nl2,nl3,nl4)

      Pi=acos(-1.d0)   !Pi
      G=6.674d-11 !N m^2 kg^-2  Gravitation constamt
      Rsun=696265.0d0*1000.0d0 !m  radius of Sun
      Rjup=142980.0d0*1000.0d0/2.0d0 !m  radius of Jupiter
      Msun=1.9891d30 !kg  mass of Sun
      Mearth=5.974d24 !kg mass of Earth
      Mjup=317.833d0*Mearth !kg  mass of Jupiter
      aConst=(G/(4.0*Pi*Pi))**(1.0d0/3.0d0)
      nplanet=7

      if(iargc().lt.1) goto 901 !number of required commandline pars
      call getarg(1,parsfile) !get first parameter - inputs
      
      nc=0 !counter
      nsc=0 !count number of stars
      nunit=10
      open(unit=nunit,file=parsfile,status='old',err=903)
c      call getplotpars(nunit,filename,ntype,fitparsfile,id,id2,kID,
c     .  Kmag,Teff,logg,rad,eoff,escale,xwidth,bins)
c      read(nunit,*) dumc
c      read(nunit,*) dumc

      plotsize=50.0
c      plotsize=112.0!105.0!76.5!95
      call pgopen('?')
      call pgpage()
c      call pgpap(7.0,1.0)
      call pgpap(10.25,1.0)
      sscale=1.5
      bounds(1)=0.0 !-sscale*srad
      bounds(2)=plotsize! sscale*srad
      bounds(3)=0.0 !-sscale*srad
      bounds(4)=plotsize! sscale*srad
      call pgwindow(bounds(1),bounds(2),bounds(3),bounds(4))
      
      call pgscr(14,0.0,0.0,0.0)
      call pgsci(14)
      call pgrect(0.0,plotsize,0.0,plotsize)
      
      CI1 = 16 
      CI2 = 86
      call PGSCIR(CI1,CI2)
      shiftx=0.0
      kidold=0
      
c      goto 20
      
      x=0.0
      y=20.7!5.7
      y=8.7!40.7
      ymax=0
      ymin=100.0
 14   read(nunit,*,end=13) koi,kid,Rstar,logg,Teff,b,rprs,flag
c         write(0,*) "KOI: ",koi,Rstar,rprs
!         if(rprs*Rstar*109.17.gt.30.0) goto 14
        if(rprs.gt.0.5) goto 14
        if(flag.ne.0) goto 14
        nc=nc+1
        
        srad=real(Rstar)
        sradold=0.0
        rp=real(rprs*Rstar)
        br=min(real(b),0.99)

        if(kid.eq.kidold) then
         np=np+1
        else
         np=1
         nsc=nsc+1 !count total number of stars
        endif

         i=np
         dth=(Pi-2.0*asin(br))/real(nplanet+2)
         th=asin(br)+dth*i
         dx= br*srad/tan(th)
         dy=-br*srad


        if(kid.eq.kidold) then
            shiftx=-1.0*(shiftx+0.2)
            shifty=-shifty
            goto 16 !only plot planet
        else
            shiftx=0.0
            shifty=1.0
        endif

       write(0,*) dx,dy

 
        x=x+srad*1.1
        ymax=max(srad,ymax)
        ymin=min(srad,ymin)
        if(x+srad.gt.plotsize)then
            x=srad
            y=y+ymax+1.2*ymin
            ymax=0
        endif
        sradold=srad
c        write(0,*) x,y
 
C koi,kid,Kmag,T0,dumr,Per,dumr,adrs,dumr,
C     .  rprs,dumr,b,dumr,tdur,Teff,logg,Rstar,asemi,Teq,tdepth
     
c      write(0,*) Rstar,Teff,rprs
      
        call starcolour(Teff,In,ncolour,colourtemp,sR,sG,sB)
      
        write(0,*) "KID",kID,rstar,nc
      
c      open(unit=nunit,file=filename,status='old',err=902)  
c      if (ntype.eq.0)then
c        call readdata(nunit,nmax,npt,time,mag,merr,itime,Zerotime)
c      elseif(ntype.eq.1)then
c        Zerotime=54900.0
c        call readkeplc(nunit,nmax,npt,time,mag,merr,itime,Zerotime)
c      endif
c      close(nunit)
      
c      open(unit=nunit,file=fitparsfile,status='old',err=904)
c      call getfitpars(nunit,nfit,sol,serr,Dpvary,err,doe,toff,eoff)
c      close(nunit)
            
C     Get limbdarkening
        if(Teff.le.13000.0)then
            inl(1)=(Teff-3500)/250+1
        elseif(Teff.gt.13000)then
            inl(1)=(Teff-13000)/1000+39
        endif
        inl(2)=int(10.0*logg)/5+1      
        inl(3)=6 ![Fe/H]=0
      
        nl(1)=real(nl1(inl(1),inl(2),inl(3)))
        nl(2)=real(nl2(inl(1),inl(2),inl(3)))
        nl(3)=real(nl3(inl(1),inl(2),inl(3)))
        nl(4)=real(nl4(inl(1),inl(2),inl(3)))
c        write(0,*) (nl(i),i=1,4)
            
        do 12 i=ci1,ci2
            Intensity=real(i-ci1)/real(ci2-ci1)
            call pgscr(i,In(1)*Intensity/255.0,In(2)*Intensity/255.0,
     .          In(3)*Intensity/255.0)     
 12     continue

C        Drawing the Sun      
c        if(flag.lt.0)then 
cc            x=x+0.5
c            call pgsci(2)
c            call pgcirc(x,y,srad*1.2)
c        endif
      
        do 10 i=ci2,ci1,-1
            r=real(i-ci1)/real(ci2-ci1)
            mu=sqrt(1.0-r*r)
            Intensity=1.0
            do 11 j=1,4
                Intensity=Intensity-nl(j)*(1.0-mu**(real(j)/2.0))
 11         continue
c           write(0,*) i,r,Intensity
            icol=Intensity*(ci2-ci1)+ci1
            call pgsci(icol)
            call pgcirc(x,y,r*srad)
 10     continue     

C       Plot planet
 16     call pgscr(14,0.0,0.0,0.0)
        call pgsci(14)
        if(br.lt.0.9)then
            call pgcirc(x-(1.0+shiftx)*(srad*cos(br)/2.0),
     .          y+shifty*srad*br,rp)
        else
            call pgcirc(x,y+shifty*srad*br,rp)
        endif
c        call pgcirc(x+dx,y+dy,rp)
        call pgsci(1)
  
        kidold=kid !remember previous KID
      goto 14
 13   continue
 
C     Draw the Sun
      call drawSun(ncolour,colourtemp,sR,sG,sB,nl1,nl2,nl3,nl4,ci1,ci2)

      
c      call pgpage()
 20     call drawscale(ncolour,colourtemp,sR,sG,sB,ci1,ci2)
      
      call pgclos()
      
      close(nunit)
      
      write(0,*) "nc, nsc", nc, nsc

      goto 999
 901  write(6,*) "Usage: transitplot <parsfile>"
      goto 999
 902  write(6,*) "Cannot open ",filename
      goto 999
 903  write(6,*) "Cannot open ",parsfile
      goto 999
 904  write(6,*) "Cannot open ",fitparsfile
      goto 999
 905  write(6,*) "Cannot open starcolourRGB.txt"
      goto 999
 999  end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine drawscale(ncolour,colourtemp,sR,sG,sB,ci1,ci2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer ncolour,colourtemp(ncolour),sR(ncolour),sG(ncolour),
     .  sB(ncolour),i,ci1,ci2,Teff,ia(100,1)
      real In(3)


C     Create the colour map
      do 12 i=ci1,ci2
        Teff=int(3000.0+(10000.0-3000.0)*dble(i-ci1)/dble(ci2-ci1))

        call starcolour(Teff,In,ncolour,colourtemp,sR,sG,sB)
        call pgscr(i,(In(1)/255.0),
     .               (In(2)/255.0),
     .               (In(3)/255.0))
        ia(i-ci1+1,1)=i 
c        write(0,*) "Teff:",Teff,i-ci1+1,ia(i-ci1+1,1)
 12   continue

C     Create the gradient
      call pgsvp(0.05,0.95,0.05,0.95)
      call pgwnad(0.0,1.0,0.0,1.0)
      call pgpixl(ia,100,1,1,ci2-ci1,1,1,0.0,1.0,0.0,0.05)
      
     
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine drawSun(ncolour,colourtemp,sR,sG,sB,nl1,nl2,nl3,nl4,
     .  ci1,ci2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer Teff,ncolour,colourtemp(ncolour),sR(ncolour),sG(ncolour),
     .  sB(ncolour),inl(3),i,ci1,ci2,j,icol
      real In(3),nl(4),logg,Intensity,r,mu,srad,x,y,b,rp
      double precision nl1(76,11,8),nl2(76,11,8),nl3(76,11,8),
     .  nl4(76,11,8)
      
      Teff=5781
      logg=4.438
      srad=1
      x=80.0
      y=34.0
      
      call starcolour(Teff,In,ncolour,colourtemp,sR,sG,sB)
      
      if(Teff.le.13000.0)then
        inl(1)=(Teff-3500)/250+1
      elseif(Teff.gt.13000)then
        inl(1)=(Teff-13000)/1000+39
      endif
      inl(2)=int(10.0*logg)/5+1      
      inl(3)=6 ![Fe/H]=0
      
      nl(1)=real(nl1(inl(1),inl(2),inl(3)))
      nl(2)=real(nl2(inl(1),inl(2),inl(3)))
      nl(3)=real(nl3(inl(1),inl(2),inl(3)))
      nl(4)=real(nl4(inl(1),inl(2),inl(3)))
      
      do 12 i=ci1,ci2
        Intensity=real(i-ci1)/real(ci2-ci1)
        call pgscr(i,(In(1)*Intensity/255.0),
     .               (In(2)*Intensity/255.0),
     .               (In(3)*Intensity/255.0))    
 12   continue
      
      do 10 i=ci2,ci1,-1
        r=real(i-ci1)/real(ci2-ci1)
        mu=sqrt(1.0-r*r)
        Intensity=1.0
        do 11 j=1,4
            Intensity=Intensity-nl(j)*(1.0-mu**(real(j)/2.0))
 11     continue
c       write(0,*) i,r,Intensity
        icol=Intensity*(ci2-ci1)+ci1
        call pgsci(icol)
        call pgcirc(x,y,r*srad)
 10   continue         

C     Plot Earth
      b=0.0
      rp=0.0092*srad
      call pgscr(14,0.0,0.0,0.0)
      call pgsci(14)
      if(b.lt.0.9)then
        call pgcirc(x-(srad*cos(b)/2.0),y+srad*b,rp)
      else
        call pgcirc(x,y+srad*b,rp)
c        write(0,*) x,y+shifty*srad*br
      endif
      call pgsci(1)
      
C     Plot Jupiter
      b=0.5
      rp=0.0892*srad
      call pgscr(14,0.0,0.0,0.0)
      call pgsci(14)
      if(b.lt.0.9)then
        call pgcirc(x-(srad*cos(b)/2.0),y+srad*b,rp)
      else
        call pgcirc(x,y+srad*b,rp)
c        write(0,*) x,y+shifty*srad*br
      endif
      call pgsci(1)
        
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine starcolour(Teff,In,ncolour,colourtemp,sR,sG,sB)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer Teff,ncolour,colourtemp(ncolour),sR(ncolour),sG(ncolour),
     .  sB(ncolour),i
      real In(3)
      logical loop
 
      if(Teff.ge.colourtemp(1))then
        In(1)=real(sR(1))
        In(2)=real(sG(1))
        In(3)=real(sB(1))
      elseif(Teff.le.colourtemp(ncolour))then
        In(1)=real(sR(ncolour))
        In(2)=real(sG(ncolour))
        In(3)=real(sB(ncolour))
      else
        loop=.true.
        i=1
        do while(loop)
           if((Teff.lt.colourtemp(i)).and.(Teff.ge.colourtemp(i+1)))then
                In(1)=real(sR(i))
                In(2)=real(sG(i))
                In(3)=real(sB(i))
                loop=.false.
            endif
            i=i+1
        enddo
      endif
c      if(Teff.gt.30000)then
c        In(1)=155.0
c        In(2)=176.0
c        In(3)=255.0
c      elseif(Teff.gt.9250)then
c        In(1)=170.0
c        In(2)=191.0
c        In(3)=255.0
c      elseif(Teff.gt.7200)then
c        In(1)=202.0
c        In(2)=215.0
c        In(3)=255.0
c      elseif(Teff.gt.5250)then
c        In(1)=255.0
c        In(2)=244.0
c        In(3)=234.0
c      elseif(Teff.gt.3850)then
c        In(1)=255.0
c        In(2)=210.0
c        In(3)=161.0
c      else
c        In(1)=255.0
c        In(2)=204.0
c        In(3)=111.0
c      endif
        
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine limbdarkensetup(nl1,nl2,nl3,nl4)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer tabteff,tablogg,tabfeh,i,j,k
      double precision dumc,nl(4),temp(7),nl1(76,11,8),nl2(76,11,8),
     .  nl3(76,11,8),nl4(76,11,8)
      character*80 file1,file2

c      file1="/home/rowe/p555/transitfit/kepler.list"
      file1="/Users/rowe/Documents/transitfit/models.list"
c      file2="/home/rowe/p555/transitfit/kepler.new"
      file2="/Users/rowe/Documents/transitfit/kepler.ch_25_26.ld"
c      file2="/home/rowe/p555/transitfit/table.most"
      
      open(unit=12,file=file1,status='old',err=902)
      open(unit=13,file=file2,status='old',err=903)
      do 5 i=1,7
        read(13,502) dumc
 5    continue
 502  format(A1)
      
 10   read(12,501,end=11) tabteff,tablogg,tabfeh
 501  format(I5,1X,I2,1X,I3)
      read(13,*,end=11) (temp(i),i=1,7),(nl(j),j=1,4)
C       get array position 
        if(tabteff.le.13000)then
            i=(tabteff-3500)/250+1
c            write(6,*) i,temps(i),tabteff
        elseif(tabteff.gt.13000)then
            i=(tabteff-13000)/1000+39
c            write(6,*) i,temps(i),tabteff
        endif
        j=tablogg/5+1
c        write(6,*) j,loggs(j),tablogg
        if(tabfeh.le.0)then
            k=(tabfeh+25)/5+1
c            write(6,*) k,fehs(k),tabfeh
        elseif(tabfeh.eq.2)then
            k=7
        elseif(tabfeh.eq.5)then
            k=8
        endif
        nl1(i,j,k)=nl(1)
        nl2(i,j,k)=nl(2)
        nl3(i,j,k)=nl(3)
        nl4(i,j,k)=nl(4)
c        write(6,*) nl(1),nl(2),nl(3),nl(4)
c        read(5,*)
      goto 10
 11   continue
     
      close(12)
      close(13)
      
      goto 999
 901  write(0,*) "Usage: limbdarkening teff logg"
      write(0,*) "Where Teff is in K and logg is cgs"
      goto 999
 902  write(0,*) "Cannot open ",file1
      goto 999
 903  write(0,*) "Cannot open ",file2
      goto 999
 999  return
      end
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine getplotpars(nunit,filename,ntype,fitparsfile,id,id2,
     .  kID,Kmag,Teff,logg,rad,eoff,escale,xwidth,bins)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer nunit,ntype,id,kID,Teff,id2,bins
      double precision Kmag,logg,rad,eoff,escale,xwidth
      character*80 filename,fitparsfile,cline
      
      read(nunit,*) filename  !1
      read(nunit,*) ntype     !2
      read(nunit,*) fitparsfile !3
      read(nunit,500) cline 
 500  format(A80)
      read(cline,*,err=10,end=10) id,id2 !4
      goto 11
 10     read(cline,*) id
        id2=1
 11   continue
      read(nunit,*) kID  !5
      read(nunit,*) Kmag !6
      read(nunit,*) Teff !7
      read(nunit,*) logg !8
      read(nunit,*) rad  !9
      read(nunit,*) eoff !10
      read(nunit,*) escale  !11
      read(nunit,*) xwidth !12
      bins=0
      read(nunit,*,end=12) bins !13
      write(6,*)"BINS:",bins
 12   continue
      
      return
      end     
